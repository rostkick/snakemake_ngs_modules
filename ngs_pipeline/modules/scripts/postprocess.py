import pandas as pd
import numpy as np
import sys
import os


log_file = snakemake.log[0]
input_file = snakemake.input.tsv
mart_file = snakemake.params.mart
rank_file = snakemake.params.rank
output_file = snakemake.output.tsv
bed_file = snakemake.params.get('bed_file', None)
ngs_type = snakemake.params.get('ngs_type', None)

gnomad_filter_enabled = snakemake.params.get('gnomad_filter', False)
gnomad_af_threshold = snakemake.params.get('gnomad_af_threshold', 0.01)
consequence_filter_enabled = snakemake.params.get('consequence_filter', False)

EXTREME_DP_THRESHOLD = 5
EXTREME_GQ_THRESHOLD = 10

def clean_text(text):
    if isinstance(text, str):
        return text.replace('\xa0', ' ')
    return text


def remove_refcall_variants(df):
    initial_count = len(df)
    refcall_mask = pd.Series([False] * len(df), index=df.index)
    
    for col in df.columns:
        col_str = str(col).upper()
        if any(keyword in col_str for keyword in ['ZYG', 'GT', 'GENOTYPE', 'CALL']):
            col_values = df[col].astype(str)
            refcall_mask |= (
                col_values.str.contains('^0/0$') |
                col_values.str.contains('^0\|0$') |
                col_values.str.contains('^REF$', case=False) |
                col_values.str.contains('^HOM_REF$', case=False) |
                col_values.str.contains('^0$') |
                col_values.str.contains('^\./\.$') |
                col_values.str.contains('^\.\|\.$')
            )
    
    if 'AlleleDepth' in df.columns:
        try:
            df_temp = df.copy()
            split_values = df_temp['AlleleDepth'].str.split(',', expand=True)
            if split_values.shape[1] >= 2:
                ad1 = pd.to_numeric(split_values[0], errors='coerce').fillna(0)
                ad2 = pd.to_numeric(split_values[1], errors='coerce').fillna(0)
                ref_from_ad = (ad2 == 0) & (ad1 > 0)
                refcall_mask |= ref_from_ad
        except:
            pass
    
    df = df[~refcall_mask]
    return df


def remove_extreme_low_quality_variants(df):
    initial_count = len(df)
    extreme_mask = pd.Series([False] * len(df), index=df.index)
    
    if 'CoverageDepth' in df.columns:
        df['CoverageDepth'] = pd.to_numeric(df['CoverageDepth'], errors='coerce')
        extreme_dp_mask = df['CoverageDepth'] <= EXTREME_DP_THRESHOLD
        extreme_mask |= extreme_dp_mask
    
    if 'GenotypeQual' in df.columns:
        df['GenotypeQual'] = pd.to_numeric(df['GenotypeQual'], errors='coerce')
        extreme_gq_mask = df['GenotypeQual'] <= EXTREME_GQ_THRESHOLD
        extreme_mask |= extreme_gq_mask
    
    df = df[~extreme_mask]
    return df


def apply_quality_filters(df):
    if 'GATK_FILTER' not in df.columns:
        df['GATK_FILTER'] = 'PASS'
    
    if 'CoverageDepth' in df.columns:
        df['CoverageDepth'] = pd.to_numeric(df['CoverageDepth'], errors='coerce')
        dp_mask = (df['CoverageDepth'] < 10) & (df['CoverageDepth'] > EXTREME_DP_THRESHOLD)
        df.loc[dp_mask, 'GATK_FILTER'] = 'DP'

    if 'GenotypeQual' in df.columns:
        df['GenotypeQual'] = pd.to_numeric(df['GenotypeQual'], errors='coerce')
        gq_mask = (df['GenotypeQual'] < 20) & (df['GenotypeQual'] > EXTREME_GQ_THRESHOLD) & (df['GATK_FILTER'] == 'PASS')
        df.loc[gq_mask, 'GATK_FILTER'] = 'GQ'
    
    if 'AlleleDepth' in df.columns:
        df[['AD1', 'AD2']] = df['AlleleDepth'].str.split(',', expand=True).iloc[:,:2]
        df['AD1'] = pd.to_numeric(df['AD1'], errors='coerce').fillna(0)
        df['AD2'] = pd.to_numeric(df['AD2'], errors='coerce').fillna(0)
        
        df['AlleleBalance'] = df['AD2'] / (df['AD1'] + df['AD2'])
        df['AlleleBalance'] = df['AlleleBalance'].fillna(0)

        df['GT_inferred'] = 'unknown'
        df.loc[df['AlleleBalance'] <= 0.1, 'GT_inferred'] = 'hom_ref'
        df.loc[(df['AlleleBalance'] >= 0.25) & (df['AlleleBalance'] <= 0.75), 'GT_inferred'] = 'het'
        df.loc[df['AlleleBalance'] >= 0.9, 'GT_inferred'] = 'hom_alt'

        suspicious_ab_mask = ~(
            (df['GT_inferred'] == 'hom_ref') |
            (df['GT_inferred'] == 'het') |
            (df['GT_inferred'] == 'hom_alt')
        )
        
        ab_apply_mask = suspicious_ab_mask & (df['GATK_FILTER'] == 'PASS')
        df.loc[ab_apply_mask, 'GATK_FILTER'] = 'AlleleBalance'
    
    return df


def filter_by_panel_genes(df, bed_file, ngs_type):
    """Filter dataframe to only include genes from the panel (only for panel sequencing)"""
    if ngs_type != 'panel':
        print(f"NGS type is '{ngs_type}', skipping panel gene filter")
        return df
    
    if not bed_file or not os.path.exists(bed_file):
        print("BED file not provided or does not exist, skipping panel filter")
        return df
    
    panel_genes = []
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('browser') and not line.startswith('track') and not line.startswith('#'):
                fields = line.split('\t')
                if len(fields) >= 4:
                    panel_genes.append(fields[3])
    
    panel_genes_set = set(panel_genes)
    print(f"Loaded {len(panel_genes_set)} genes from panel BED file")
    
    gene_col = None
    for col in ['RefGene', 'HGNC symbol', 'Gene', 'Gene_Symbol', 'SYMBOL']:
        if col in df.columns:
            gene_col = col
            break
    
    if not gene_col:
        print("WARNING: Could not find gene column, skipping panel filter")
        return df
    
    initial_count = len(df)
    df = df[df[gene_col].isin(panel_genes_set)]
    filtered_count = initial_count - len(df)
    
    print(f"Panel filter: removed {filtered_count} variants ({initial_count} -> {len(df)})")
    
    return df


def main():
    df = pd.read_csv(input_file, sep='\t', encoding='utf-8')
    df_mart = pd.read_csv(mart_file, sep='\t', encoding='utf-8')

    df_rank = None
    if rank_file and os.path.exists(rank_file):
        df_rank = pd.read_csv(rank_file, sep='\t', encoding='utf-8')
    
    str_columns = df.select_dtypes(include=['object']).columns
    for col in str_columns:
        df[col] = df[col].apply(clean_text)
    
    df = remove_refcall_variants(df)
    df = remove_extreme_low_quality_variants(df)
    
    if consequence_filter_enabled:
        consequences_to_keep = [
            'missense_variant',
            'inframe_deletion',
            'inframe_insertion',
            'feature_truncation',
            'feature_elongation',
            'transcript_amplification',
            'start_lost',
            'stop_lost',
            'frameshift_variant',
            'stop_gained',
            'stop_codon_variant',
            'splice_acceptor_variant',
            'transcript_ablation'
        ]
        df = df[df['Consequence'].isin(consequences_to_keep)]

    gnomad_cols = [col for col in df.columns if 'gnomAD' in col]
    if gnomad_cols:
        df[gnomad_cols] = df[gnomad_cols].replace('.', 0)
        df[gnomad_cols] = df[gnomad_cols].astype(float)

    if 'RefGene' in df.columns and 'HGNC symbol' in df_mart.columns:
        df = pd.merge(df, df_mart, left_on='RefGene', right_on='HGNC symbol', how='left')

    if df_rank is not None:
        gene_symbol_col = None
        for col in ['HGNC symbol', 'RefGene', 'Gene', 'Gene_Symbol', 'SYMBOL']:
            if col in df.columns:
                gene_symbol_col = col
                break
        
        if gene_symbol_col:
            df = pd.merge(df, df_rank, left_on=gene_symbol_col, right_on='hgnc_symbol', how='left')
            if 'hgnc_symbol' in df.columns and gene_symbol_col != 'hgnc_symbol':
                df.drop('hgnc_symbol', axis=1, inplace=True)
        else:
            df['rank_gse71613'] = np.nan
    else:
        df['rank_gse71613'] = np.nan

    df = apply_quality_filters(df)

    if gnomad_filter_enabled:
        gnomad_cols = [col for col in ['gnomAD_exome_NFE', 'gnomAD_genome_NFE'] if col in df.columns]
        if gnomad_cols:
            gnomad_mask = pd.Series([False] * len(df), index=df.index)
            for col in gnomad_cols:
                col_mask = df[col] >= gnomad_af_threshold
                gnomad_mask = gnomad_mask | col_mask
            df = df[~gnomad_mask]

    df = filter_by_panel_genes(df, bed_file, ngs_type)

    consequence_order = {
        'sequence_variant': 0,
        'intergenic_variant': 1,
        'regulatory_region_variant': 2,
        'regulatory_region_amplification': 3,
        'regulatory_region_ablation': 4,
        'TF_binding_site_variant': 5,
        'TFBS_amplification': 6,
        'TFBS_ablation': 7,
        'downstream_gene_variant': 8,
        'upstream_gene_variant': 9,
        'coding_transcript_variant': 10,
        'non_coding_transcript_variant': 11,
        'NMD_transcript_variant': 12,
        'intron_variant': 13,
        'non_coding_transcript_exon_variant': 14,
        '3_prime_UTR_variant': 15,
        '5_prime_UTR_variant': 16,
        'mature_miRNA_variant': 17,
        'coding_sequence_variant': 18,
        'synonymous_variant': 19,
        'stop_retained_variant': 20,
        'start_retained_variant': 21,
        'incomplete_terminal_codon_variant': 22,
        'splice_polypyrimidine_tract_variant': 23,
        'splice_donor_region_variant': 24,
        'splice_region_variant': 25,
        'splice_donor_5th_base_variant': 26,
        'protein_altering_variant': 27,
        'missense_variant': 28,
        'inframe_deletion': 29,
        'inframe_insertion': 30,
        'feature_truncation': 31,
        'feature_elongation': 32,
        'transcript_amplification': 33,
        'start_lost': 34,
        'stop_lost': 35,
        'frameshift_variant': 36,
        'stop_gained': 37,
        'stop_codon_variant': 38,
        'splice_acceptor_variant': 39,
        'transcript_ablation': 40
    }

    if 'Chr' in df.columns and 'Consequence' in df.columns:
        df['Consequence_severity'] = df['Consequence'].map(consequence_order)
        df = df.sort_values('Consequence_severity', ascending=False).drop_duplicates(subset='Chr', keep='first')

    cols_to_drop = ['AD1', 'AD2', 'AlleleBalance', 'GT_inferred', 'Consequence_severity']
    if 'HGNC symbol' in df.columns:
        cols_to_drop.extend(['HGNC symbol', 'NCBI gene (formerly Entrezgene) ID', 
                           'Gene name', 'BIOTYPE', 'CANONICAL', 'Gene stable ID'])
    
    cols_to_drop = [col for col in cols_to_drop if col in df.columns]
    df.drop(cols_to_drop, axis=1, inplace=True)

    if 'rsID' in df.columns:
        columns = list(df.columns)
        cols_to_move = ['CoverageDepth', 'GenotypeQual', 'AlleleDepth']
        rsid_index = columns.index('rsID') if 'rsID' in columns else -1
        
        if rsid_index != -1:
            for col in cols_to_move:
                if col in columns:
                    columns.remove(col)
            
            insert_position = rsid_index + 1
            for col in reversed(cols_to_move):
                if col in df.columns:
                    columns.insert(insert_position, col)
            
            df = df[columns]
    
    df.to_csv(output_file, index=False, sep='\t', na_rep='.')

if __name__ == "__main__":
    with open(log_file, "w") as f:
        sys.stderr = sys.stdout = f
        main()