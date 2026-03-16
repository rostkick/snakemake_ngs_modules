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
# GQ threshold removed: GLnexus DeepVariantWES caps GQ at 10 during joint calling.
# Variants with GQ<=10 that passed joint calling are legitimate — they are mostly
# common polymorphisms filtered downstream by gnomAD, not true low-quality calls.

CONSEQUENCE_ORDER = {
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
    'splice_donor_variant': 39,
    'transcript_ablation': 40,
}


def max_severity(csq_str):
    """Return the highest severity score for a consequence string.
    Handles compound consequences joined by '&'
    (e.g. 'missense_variant&splice_region_variant').
    Returns -1 for missing or unrecognized consequences.
    """
    if pd.isna(csq_str):
        return -1
    return max(
        (CONSEQUENCE_ORDER.get(c.strip(), -1) for c in str(csq_str).split('&')),
        default=-1
    )


def log_step(step_name, n_before, n_after):
    removed = n_before - n_after
    pct = f"{removed / n_before * 100:.1f}%" if n_before > 0 else "N/A"
    print(f"[FILTER] {step_name}: {n_before} -> {n_after}  (removed {removed}, {pct})")


def clean_text(text):
    if isinstance(text, str):
        return text.replace('\xa0', ' ')
    return text


def remove_refcall_variants(df):
    n = len(df)
    refcall_mask = pd.Series([False] * len(df), index=df.index)

    gt_cols_found = []
    for col in df.columns:
        col_upper = str(col).upper()
        if any(k in col_upper for k in ['ZYG', 'GT', 'GENOTYPE', 'CALL']):
            gt_cols_found.append(col)
            col_values = df[col].astype(str)
            m = (
                col_values.str.fullmatch(r'0[/|]0') |
                col_values.str.fullmatch('0') |
                col_values.str.fullmatch(r'\./\.') |
                col_values.str.fullmatch(r'\.\|\.') |
                col_values.str.fullmatch(r'\.[\/|]0') |   # ./0 — alt missing
                col_values.str.fullmatch(r'0[\/|]\.') |   # 0/. — alt missing
                col_values.str.fullmatch('REF', case=False) |
                col_values.str.fullmatch('HOM_REF', case=False)
            )
            print(f"[FILTER]   GT col '{col}': {m.sum()} refcall rows")
            refcall_mask |= m

    if not gt_cols_found:
        print("[FILTER]   WARNING: no genotype column found — refcall step skipped")

    if 'AlleleDepth' in df.columns:
        try:
            split_ad = df['AlleleDepth'].astype(str).str.split(',', expand=True)
            if split_ad.shape[1] >= 2:
                ad_ref = pd.to_numeric(split_ad[0], errors='coerce').fillna(0)
                ad_alt = pd.to_numeric(split_ad[1], errors='coerce').fillna(0)
                ad_refcall = (ad_alt == 0) & (ad_ref > 0)
                print(f"[FILTER]   AlleleDepth safeguard: {ad_refcall.sum()} additional refcall rows")
                refcall_mask |= ad_refcall
        except Exception as e:
            print(f"[FILTER]   AlleleDepth safeguard failed: {e}")

    df = df[~refcall_mask]
    log_step("remove_refcall_variants", n, len(df))
    return df


def remove_extreme_low_quality_variants(df):
    n = len(df)
    extreme_mask = pd.Series([False] * len(df), index=df.index)

    if 'CoverageDepth' in df.columns:
        df['CoverageDepth'] = pd.to_numeric(df['CoverageDepth'], errors='coerce')
        m = df['CoverageDepth'] <= EXTREME_DP_THRESHOLD
        print(f"[FILTER]   DP <= {EXTREME_DP_THRESHOLD}: {m.sum()} rows")
        extreme_mask |= m

    # Hard-remove variants with GQ < 5: even accounting for GLnexus GQ cap at 10,
    # GQ 0-4 indicates DeepVariant itself has no confidence in the genotype call.
    # GQ 5-10 are retained — these are legitimate calls capped by GLnexus DeepVariantWES.
    if 'GenotypeQual' in df.columns:
        df['GenotypeQual'] = pd.to_numeric(df['GenotypeQual'], errors='coerce')
        m = df['GenotypeQual'] < 5
        print(f"[FILTER]   GQ < 5: {m.sum()} rows")
        extreme_mask |= m

    # Hard-remove variants where ref AD is missing (.) — unreliable multiallelic
    # calls where GLnexus could not compute ref allele depth.
    if 'AlleleDepth' in df.columns:
        try:
            split_ad = df['AlleleDepth'].astype(str).str.split(',', expand=True)
            ref_missing = split_ad[0].str.strip() == '.'
            print(f"[FILTER]   AD ref missing (.): {ref_missing.sum()} rows")
            extreme_mask |= ref_missing
        except Exception as e:
            print(f"[FILTER]   AD ref missing check failed: {e}")

    # Hard-remove het calls with VAF < 0.1: fewer than 10% alt reads for a het
    # is almost certainly a miscall or alignment artifact regardless of DP.
    if 'AlleleDepth' in df.columns and 'Zyg' in df.columns:
        try:
            split_ad = df['AlleleDepth'].astype(str).str.split(',', expand=True)
            ad_ref = pd.to_numeric(split_ad[0], errors='coerce').fillna(0)
            ad_alt = pd.to_numeric(split_ad[1], errors='coerce').fillna(0)
            total = ad_ref + ad_alt
            vaf = ad_alt / total.replace(0, float('nan'))
            het_mask = df['Zyg'].astype(str).str.fullmatch(r'0[/|]1')
            low_vaf_het = het_mask & (vaf < 0.10)
            print(f"[FILTER]   het VAF < 0.10: {low_vaf_het.sum()} rows")
            extreme_mask |= low_vaf_het
        except Exception as e:
            print(f"[FILTER]   VAF filter failed: {e}")

    df = df[~extreme_mask]
    log_step("remove_extreme_low_quality", n, len(df))
    return df


def apply_quality_flags(df):
    if 'GATK_FILTER' not in df.columns:
        df['GATK_FILTER'] = 'PASS'

    # Normalize GATK_FILTER: VCF PASS is represented as '.' in bcftools output
    df['GATK_FILTER'] = df['GATK_FILTER'].astype(str).replace('.', 'PASS')

    if 'CoverageDepth' in df.columns:
        df['CoverageDepth'] = pd.to_numeric(df['CoverageDepth'], errors='coerce')
        m = (df['CoverageDepth'] > EXTREME_DP_THRESHOLD) & (df['CoverageDepth'] < 10)
        df.loc[m, 'GATK_FILTER'] = 'DP'
        print(f"[FILTER]   Flagged DP (5-10): {m.sum()} rows")

    # GQ soft flag: GQ 5-9 are retained but flagged — GLnexus caps GQ at 10,
    # so these are borderline calls worth manual review.
    if 'GenotypeQual' in df.columns:
        df['GenotypeQual'] = pd.to_numeric(df['GenotypeQual'], errors='coerce')
        m = (df['GenotypeQual'] >= 5) & (df['GenotypeQual'] < 10) & (df['GATK_FILTER'] == 'PASS')
        df.loc[m, 'GATK_FILTER'] = 'GQ'
        print(f"[FILTER]   Flagged GQ (5-9): {m.sum()} rows")

    if 'AlleleDepth' in df.columns:
        df[['AD1', 'AD2']] = df['AlleleDepth'].astype(str).str.split(',', expand=True).iloc[:, :2]
        df['AD1'] = pd.to_numeric(df['AD1'], errors='coerce').fillna(0)
        df['AD2'] = pd.to_numeric(df['AD2'], errors='coerce').fillna(0)
        total_ad = df['AD1'] + df['AD2']
        df['AlleleBalance'] = np.where(total_ad > 0, df['AD2'] / total_ad, np.nan)

        df['GT_inferred'] = 'unknown'
        df.loc[df['AlleleBalance'] <= 0.10, 'GT_inferred'] = 'hom_ref'
        df.loc[(df['AlleleBalance'] >= 0.25) & (df['AlleleBalance'] <= 0.75), 'GT_inferred'] = 'het'
        df.loc[df['AlleleBalance'] >= 0.90, 'GT_inferred'] = 'hom_alt'

        m = (df['GT_inferred'] == 'unknown') & (df['GATK_FILTER'] == 'PASS')
        df.loc[m, 'GATK_FILTER'] = 'AlleleBalance'
        print(f"[FILTER]   Flagged AlleleBalance: {m.sum()} rows")

    return df


def apply_gnomad_filter(df, threshold):
    """Filter common variants by gnomAD NFE allele frequency.

    Logic for WES:
    - Filter ONLY on gnomAD_exome_NFE >= threshold
    - Variants with exome_NFE == 0 (absent from exome DB) are KEPT —
      they may be rare variants not captured by the gnomAD exome project.
    - gnomAD_genome_NFE is NOT used for filtering: it is populated for
      nearly every genomic position and would remove ~94% of variants.
    - gnomAD_genome_NFE is retained in the output for manual review.
    """
    n = len(df)
    exome_col  = 'gnomAD_exome_NFE'
    genome_col = 'gnomAD_genome_NFE'

    if exome_col not in df.columns and genome_col not in df.columns:
        print("[FILTER]   gnomAD filter: no gnomAD columns found, skipping")
        return df

    if exome_col in df.columns:
        mask = df[exome_col] >= threshold
        print(f"[FILTER]   {exome_col} >= {threshold}: {mask.sum()} rows removed")
        print(f"[FILTER]   {exome_col} == 0 (absent from exome DB, kept): {(df[exome_col] == 0).sum()} rows")
        if genome_col in df.columns:
            print(f"[FILTER]   {genome_col}: reported but not used for filtering")
    else:
        # exome column absent — fall back to genome
        mask = df[genome_col] >= threshold
        print(f"[FILTER]   {genome_col} >= {threshold} (exome col absent): {mask.sum()} rows removed")

    df = df[~mask]
    log_step(f"gnomad_filter (exome_NFE >= {threshold})", n, len(df))
    return df


def apply_consequence_filter(df):
    n = len(df)
    # Coding consequences for wes_clinical:
    # includes all functionally relevant coding + splicing variants
    # excludes: intron, downstream, upstream, UTR, non_coding_transcript
    consequences_to_keep = {
        # protein-altering
        'missense_variant',
        'inframe_deletion',
        'inframe_insertion',
        'frameshift_variant',
        'stop_gained',
        'stop_lost',
        'stop_codon_variant',
        'start_lost',
        'start_retained_variant',
        'transcript_ablation',
        'transcript_amplification',
        'feature_truncation',
        'feature_elongation',
        'protein_altering_variant',
        # splicing
        'splice_acceptor_variant',
        'splice_donor_variant',
        'splice_region_variant',
        'splice_donor_5th_base_variant',
        'splice_donor_region_variant',
        'splice_polypyrimidine_tract_variant',
        # synonymous — keep for completeness of coding region
        'synonymous_variant',
        'stop_retained_variant',
        'coding_sequence_variant',
        'incomplete_terminal_codon_variant',
        # NMD
        'NMD_transcript_variant',
    }

    def has_kept_csq(csq):
        if pd.isna(csq):
            return False
        return any(c.strip() in consequences_to_keep for c in str(csq).split('&'))

    df = df[df['Consequence'].apply(has_kept_csq)]
    log_step("consequence_filter", n, len(df))
    return df


def dedup_per_position(df):
    """Keep the most severe consequence per CHROM:POS position.
    Handles compound consequences (e.g. missense_variant&splice_region_variant)
    by taking the max severity component.
    """
    if 'Chr' not in df.columns or 'Consequence' not in df.columns:
        return df

    n = len(df)
    df = df.copy()
    df['_sev'] = df['Consequence'].apply(max_severity)
    df = (df
          .sort_values('_sev', ascending=False)
          .drop_duplicates(subset='Chr', keep='first')
          .drop(columns=['_sev']))

    log_step("dedup_per_position (keep most severe transcript per CHROM:POS)", n, len(df))
    return df


def filter_by_panel_genes(df, bed_file, ngs_type):
    if ngs_type != 'panel':
        print(f"[FILTER] filter_by_panel_genes: skipped (ngs_type='{ngs_type}')")
        return df

    if not bed_file or not os.path.exists(bed_file):
        print("[FILTER] filter_by_panel_genes: BED not found, skipping")
        return df

    panel_genes = set()
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(('browser', 'track', '#')):
                continue
            fields = line.split('\t')
            if len(fields) >= 4:
                panel_genes.add(fields[3])

    print(f"[FILTER] filter_by_panel_genes: {len(panel_genes)} genes in BED")
    gene_col = next((c for c in ['RefGene', 'HGNC symbol', 'Gene', 'Gene_Symbol', 'SYMBOL']
                     if c in df.columns), None)
    if not gene_col:
        print("[FILTER] filter_by_panel_genes: gene column not found, skipping")
        return df

    n = len(df)
    df = df[df[gene_col].isin(panel_genes)]
    log_step("filter_by_panel_genes", n, len(df))
    return df


def main():
    df = pd.read_csv(input_file, sep='\t', encoding='utf-8', low_memory=False)
    print(f"[FILTER] ===== postprocess start =====")
    print(f"[FILTER] Input: {input_file}  ({len(df)} variants)")
    print(f"[FILTER] Settings: gnomad_filter={gnomad_filter_enabled}, "
          f"threshold={gnomad_af_threshold}, consequence_filter={consequence_filter_enabled}, "
          f"ngs_type={ngs_type}")

    df_mart = pd.read_csv(mart_file, sep='\t', encoding='utf-8')
    if 'HGNC symbol' in df_mart.columns:
        n_mart_before = len(df_mart)
        df_mart = df_mart.drop_duplicates(subset='HGNC symbol', keep='first')
        print(f"[FILTER] Mart dedup: {n_mart_before} -> {len(df_mart)} rows (1 per gene)")
    df_rank = (pd.read_csv(rank_file, sep='\t', encoding='utf-8')
               if rank_file and os.path.exists(rank_file) else None)

    for col in df.select_dtypes(include=['object']).columns:
        df[col] = df[col].apply(clean_text)

    # 1. Remove reference calls (0/0, ./. etc.)
    df = remove_refcall_variants(df)

    # 2. Hard-remove extreme low quality
    df = remove_extreme_low_quality_variants(df)

    # 3. Optional: keep only pathogenic consequence classes
    if consequence_filter_enabled:
        df = apply_consequence_filter(df)

    # 4. Parse gnomAD columns
    gnomad_cols = [c for c in df.columns if 'gnomAD' in c]
    if gnomad_cols:
        df[gnomad_cols] = df[gnomad_cols].replace('.', 0)
        df[gnomad_cols] = df[gnomad_cols].astype(float)

    # 5. Merge gene annotation (mart)
    n_pre = len(df)
    if 'RefGene' in df.columns and 'HGNC symbol' in df_mart.columns:
        df = pd.merge(df, df_mart, left_on='RefGene', right_on='HGNC symbol', how='left')
    if len(df) != n_pre:
        print(f"[FILTER] WARNING: mart merge changed row count {n_pre} -> {len(df)}")

    # 6. Merge gene rank
    if df_rank is not None:
        if 'hgnc_symbol' in df_rank.columns:
            n_rank_before = len(df_rank)
            df_rank = df_rank.drop_duplicates(subset='hgnc_symbol', keep='first')
            if len(df_rank) < n_rank_before:
                print(f"[FILTER] Rank dedup: {n_rank_before} -> {len(df_rank)} rows")
        # Prefer RefGene for join (stable, always present); fall back to others
        gene_col = next((c for c in ['RefGene', 'Gene', 'Gene_Symbol', 'SYMBOL']
                         if c in df.columns), None)
        if gene_col:
            n_pre_rank = len(df)
            df = pd.merge(df, df_rank, left_on=gene_col, right_on='hgnc_symbol', how='left')
            if 'hgnc_symbol' in df.columns:
                df.drop('hgnc_symbol', axis=1, inplace=True)
            if len(df) != n_pre_rank:
                print(f"[FILTER] WARNING: rank merge changed row count {n_pre_rank} -> {len(df)}")
        else:
            df['rank_gse71613'] = np.nan
    else:
        df['rank_gse71613'] = np.nan

    # 7. Soft quality flags (tag only, do NOT remove)
    df = apply_quality_flags(df)

    # 8. gnomAD frequency filter (wes_clinical only, enabled in snakemake rule)
    if gnomad_filter_enabled:
        df = apply_gnomad_filter(df, gnomad_af_threshold)

    # 9. Panel gene filter (panel ngs_type only)
    df = filter_by_panel_genes(df, bed_file, ngs_type)

    # 10. Dedup multi-transcript rows: keep most severe per CHROM:POS
    df = dedup_per_position(df)

    # 11. Sort by consequence severity (most severe first), then by chromosome/position
    def chrom_sort_key(chrom_pos):
        chrom = str(chrom_pos).split(':')[0].replace('chr', '').replace('Chr', '')
        order = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
        try:
            return order.get(chrom, int(chrom))
        except ValueError:
            return 99

    def pos_sort_key(chrom_pos):
        try:
            return int(str(chrom_pos).split(':')[1])
        except (IndexError, ValueError):
            return 0

    if 'Chr' in df.columns and 'Consequence' in df.columns:
        df = df.copy()
        df['_sev'] = df['Consequence'].apply(max_severity)
        df['_chrom'] = df['Chr'].apply(chrom_sort_key)
        df['_pos'] = df['Chr'].apply(pos_sort_key)
        df = df.sort_values(['_sev', '_chrom', '_pos'], ascending=[False, True, True])
        df = df.drop(columns=['_sev', '_chrom', '_pos'])

    # Cleanup temp columns
    cols_to_drop = ['AD1', 'AD2', 'AlleleBalance', 'GT_inferred']
    if 'HGNC symbol' in df.columns:
        cols_to_drop += ['HGNC symbol', 'NCBI gene (formerly Entrezgene) ID',
                         'Gene name', 'BIOTYPE', 'CANONICAL', 'Gene stable ID']
    df.drop(columns=[c for c in cols_to_drop if c in df.columns], inplace=True)

    # Move DP/GQ/AD right after rsID
    if 'rsID' in df.columns:
        cols = list(df.columns)
        to_move = [c for c in ['CoverageDepth', 'GenotypeQual', 'AlleleDepth'] if c in cols]
        for c in to_move:
            cols.remove(c)
        idx = cols.index('rsID') + 1
        for c in reversed(to_move):
            cols.insert(idx, c)
        df = df[cols]

    print(f"[FILTER] ===== Final output: {len(df)} variants =====")
    df.to_csv(output_file, index=False, sep='\t', na_rep='.')


if __name__ == "__main__":
    with open(log_file, "w") as f:
        sys.stderr = sys.stdout = f
        main()