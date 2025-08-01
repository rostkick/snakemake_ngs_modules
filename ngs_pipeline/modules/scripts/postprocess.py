import pandas as pd
import numpy as np
import sys
import os


log_file = snakemake.log[0]
input_file = snakemake.input.tsv
mart_file = snakemake.params.mart
rank_file = snakemake.params.rank
output_file = snakemake.output.tsv


def clean_text(text):
    """Clean text by replacing problematic characters"""
    if isinstance(text, str):
        return text.replace('\xa0', ' ')  # Replace non-breaking space with regular space
    return text

def main():
    # Load data
    df = pd.read_csv(input_file, sep='\t', encoding='utf-8')
    df_mart = pd.read_csv(mart_file, sep='\t', encoding='utf-8')

    # Load ranking table if provided
    df_rank = None
    if rank_file and os.path.exists(rank_file):
        df_rank = pd.read_csv(rank_file, sep='\t', encoding='utf-8')
        print(f"Loaded {len(df_rank)} genes with rankings")
    else:
        print("No ranking file provided or file not found - skipping ranking annotation")
    
    # Clean any problematic characters in string columns
    str_columns = df.select_dtypes(include=['object']).columns
    for col in str_columns:
        df[col] = df[col].apply(clean_text)

    # Initialize GATK_FILTER column if it doesn't exist
    if 'GATK_FILTER' not in df.columns:
        df['GATK_FILTER'] = 'PASS'
    
    # Get only protein coding variants
    if 'BIOTYPE' in df.columns:
        df = df[df['BIOTYPE'] == 'protein_coding']
        print(f"After protein_coding filter: {len(df)} variants")

    # Handle gnomAD columns
    gnomad_cols = [col for col in df.columns if 'gnomAD' in col]
    if gnomad_cols:
        df[gnomad_cols] = df[gnomad_cols].replace('.', 0)
        df[gnomad_cols] = df[gnomad_cols].astype(float)

    # Merge with mart data
    if 'RefGene' in df.columns and 'HGNC symbol' in df_mart.columns:
        df = pd.merge(df, df_mart, left_on='RefGene', right_on='HGNC symbol', how='left')

    # Merge with ranking data
    if df_rank is not None:
        # First, try to merge using existing gene symbol column
        gene_symbol_col = None
        for col in ['HGNC symbol', 'RefGene', 'Gene', 'Gene_Symbol', 'SYMBOL']:
            if col in df.columns:
                gene_symbol_col = col
                break
        
        if gene_symbol_col:
            print(f"Merging ranking data using column: {gene_symbol_col}")
            df = pd.merge(df, df_rank, left_on=gene_symbol_col, right_on='hgnc_symbol', how='left')
            # Clean up duplicate hgnc_symbol column from merge
            if 'hgnc_symbol' in df.columns and gene_symbol_col != 'hgnc_symbol':
                df.drop('hgnc_symbol', axis=1, inplace=True)
            
            # Count how many variants got ranking annotation
            ranked_variants = len(df[df['rank_gse71613'].notna()])
            print(f"Added ranking annotation to {ranked_variants} variants")
        else:
            print("Warning: No suitable gene symbol column found for ranking annotation")
            # Add empty ranking column
            df['rank_gse71613'] = np.nan
    else:
        # Add empty ranking column if no ranking file
        df['rank_gse71613'] = np.nan

    # Apply filters
    
    # RefCall filter - exclude reference calls (0/0 genotypes)
    refcall_mask = pd.Series([False] * len(df), index=df.index)
    
    # Check different possible genotype columns and formats
    if 'Zyg' in df.columns:
        # Standard zygosity notation
        refcall_mask |= (df['Zyg'] == '0/0') | (df['Zyg'] == '0|0') | (df['Zyg'] == 'REF')
    
    # Also check if we have GT column directly
    if 'GT' in df.columns:
        refcall_mask |= (df['GT'] == '0/0') | (df['GT'] == '0|0')
    
    # Apply RefCall filter
    if refcall_mask.any():
        df.loc[refcall_mask & (df['GATK_FILTER'] == 'PASS'), 'GATK_FILTER'] = 'RefCall'
        refcall_filtered = len(df[df['GATK_FILTER'] == 'RefCall'])
        print(f"RefCall filter: {refcall_filtered} variants filtered (reference calls)")
    
    # Coverage depth and genotype quality filters
    if 'CoverageDepth' in df.columns:
        df['CoverageDepth'] = pd.to_numeric(df['CoverageDepth'], errors='coerce')
        df.loc[df['CoverageDepth'] < 10, 'GATK_FILTER'] = 'DP'
    
    if 'GenotypeQual' in df.columns:
        df['GenotypeQual'] = pd.to_numeric(df['GenotypeQual'], errors='coerce')
        df.loc[df['GenotypeQual'] < 20, 'GATK_FILTER'] = 'GQ'

    # AlleleBalance Filter
    if 'AlleleDepth' in df.columns:
        # Parse AlleleDepth values
        df[['AD1', 'AD2']] = df['AlleleDepth'].str.split(',', expand=True).iloc[:,:2]
        df['AD1'] = pd.to_numeric(df['AD1'], errors='coerce').fillna(0)
        df['AD2'] = pd.to_numeric(df['AD2'], errors='coerce').fillna(0)
        
        # Calculate AlleleBalance (AD2 is ALT, AD1 is REF)
        df['AlleleBalance'] = df['AD2'] / (df['AD1'] + df['AD2'])
        df['AlleleBalance'] = df['AlleleBalance'].fillna(0)

        # Infer genotype based on AlleleBalance
        df['GT_inferred'] = 'unknown'
        df.loc[df['AlleleBalance'] <= 0.1, 'GT_inferred'] = 'hom_ref'  # 0/0
        df.loc[(df['AlleleBalance'] >= 0.25) & (df['AlleleBalance'] <= 0.75), 'GT_inferred'] = 'het'  # 0/1
        df.loc[df['AlleleBalance'] >= 0.9, 'GT_inferred'] = 'hom_alt'  # 1/1

        # PASS variants with expected AlleleBalance for their genotype
        pass_ab_filter = (
            (df['GT_inferred'] == 'hom_ref') |  # AB <= 0.1
            (df['GT_inferred'] == 'het') |      # AB 0.25-0.75  
            (df['GT_inferred'] == 'hom_alt')    # AB >= 0.9
        )

        # FILTER variants with suspicious AlleleBalance
        df.loc[~pass_ab_filter & (df['GATK_FILTER'] == 'PASS'), 'GATK_FILTER'] = 'AlleleBalance'

        # Set zygosity based on inferred genotype
        df.loc[df['GT_inferred'] == 'het', 'Zyg'] = 'HET'
        df.loc[df['GT_inferred'].isin(['hom_ref', 'hom_alt']), 'Zyg'] = 'HOM'
        df.loc[df['GT_inferred'] == 'unknown', 'Zyg'] = '.'

        # Additional RefCall filter based on inferred genotype
        refcall_inferred_mask = df['GT_inferred'] == 'hom_ref'
        df.loc[refcall_inferred_mask & (df['GATK_FILTER'] == 'PASS'), 'GATK_FILTER'] = 'RefCall'

        ab_filtered = len(df[df['GATK_FILTER'] == 'AlleleBalance'])
        pass_count = len(df[pass_ab_filter])
        print(f"AlleleBalance: {pass_count} variants passed, {ab_filtered} variants filtered")

    # Consequence severity ranking for multiple canonical transcripts
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

    # Remove duplicates by keeping the most severe consequence per position
    if 'Chr' in df.columns and 'Consequence' in df.columns:
        df = df.drop_duplicates(subset='Chr', keep='last')
        df['Consequence_severity'] = df['Consequence'].map(consequence_order)
        df = df.sort_values('Consequence_severity', ascending=False).drop_duplicates(subset='Chr', keep='first')

    # Print filter statistics
    print("\nFilter statistics:")
    filter_counts = df['GATK_FILTER'].value_counts()
    for filt, count in filter_counts.items():
        print(f"{filt}: {count} variants")
    
    total_passed = len(df[df['GATK_FILTER'] == 'PASS'])
    total_filtered = len(df[df['GATK_FILTER'] != 'PASS'])
    print(f"\nSummary: {total_passed} variants passed, {total_filtered} variants filtered")

    # Print ranking statistics
    if 'rank_gse71613' in df.columns:
        ranked_count = len(df[df['rank_gse71613'].notna()])
        print(f"Ranking annotation: {ranked_count} variants have ranking information")
        if ranked_count > 0:
            print(f"Ranking range: {df['rank_gse71613'].min():.1f} - {df['rank_gse71613'].max():.1f}")

    # Optional: Remove RefCall variants entirely (uncomment next line if desired)
    df = df[df['GATK_FILTER'] != 'RefCall']
    print(f"After RefCall removal: {len(df)} variants remaining")

    # Clean up temporary columns
    cols_to_drop = ['AD1', 'AD2', 'AlleleBalance', 'GT_inferred', 'Consequence_severity']
    if 'HGNC symbol' in df.columns:
        cols_to_drop.extend(['HGNC symbol', 'NCBI gene (formerly Entrezgene) ID', 
                           'Gene name', 'BIOTYPE', 'CANONICAL', 'Gene stable ID'])
    
    # Only drop columns that exist
    cols_to_drop = [col for col in cols_to_drop if col in df.columns]
    df.drop(cols_to_drop, axis=1, inplace=True)

    # Save filtered results
    df.to_csv(output_file, index=False, sep='\t', na_rep='.')
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    with open(log_file, "w") as f:
        sys.stderr = sys.stdout = f
        main()