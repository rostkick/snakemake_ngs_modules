import pandas as pd
import numpy as np
import sys
import os


log_file = snakemake.log[0]
input_file = snakemake.input.tsv
mart_file = snakemake.params.mart
rank_file = snakemake.params.rank
output_file = snakemake.output.tsv

# Optional gnomAD filtering parameters
gnomad_filter_enabled = snakemake.params.get('gnomad_filter', False)
gnomad_af_threshold = snakemake.params.get('gnomad_af_threshold', 0.01)  # Default 1%

# Optional consequence filter
consequence_filter_enabled = snakemake.params.get('consequence_filter', False)

# Extreme quality thresholds - variants below these are completely removed
EXTREME_DP_THRESHOLD = 5
EXTREME_GQ_THRESHOLD = 10

def clean_text(text):
    """Clean text by replacing problematic characters"""
    if isinstance(text, str):
        return text.replace('\xa0', ' ')  # Replace non-breaking space with regular space
    return text


def remove_refcall_variants(df):
    """
    Remove ALL reference call variants (0/0 genotypes) from the dataframe
    
    Args:
        df: DataFrame with variant data
        
    Returns:
        df: DataFrame with RefCall variants removed
    """
    initial_count = len(df)
    
    # Create mask for ALL possible reference call formats
    refcall_mask = pd.Series([False] * len(df), index=df.index)
    
    # Check all columns that might contain genotype information
    for col in df.columns:
        col_str = str(col).upper()
        
        # Check if this column might contain genotype information
        if any(keyword in col_str for keyword in ['ZYG', 'GT', 'GENOTYPE', 'CALL']):
            # Convert to string and check for reference patterns
            col_values = df[col].astype(str)
            
            # Check for all possible reference call formats
            refcall_mask |= (
                col_values.str.contains('^0/0$') |  # 0/0
                col_values.str.contains('^0\|0$') |  # 0|0
                col_values.str.contains('^REF$', case=False) |  # REF
                col_values.str.contains('^HOM_REF$', case=False) |  # HOM_REF
                col_values.str.contains('^0$') |  # 0
                col_values.str.contains('^\./\.$') |  # ./.
                col_values.str.contains('^\.\|\.$')  # .|.
            )
    
    # Also check AlleleBalance for hom_ref inference
    if 'AlleleDepth' in df.columns:
        try:
            # Parse AlleleDepth to identify reference calls
            df_temp = df.copy()
            split_values = df_temp['AlleleDepth'].str.split(',', expand=True)
            if split_values.shape[1] >= 2:
                ad1 = pd.to_numeric(split_values[0], errors='coerce').fillna(0)
                ad2 = pd.to_numeric(split_values[1], errors='coerce').fillna(0)
                
                # Variants with only reference allele (AD2 = 0) are reference calls
                ref_from_ad = (ad2 == 0) & (ad1 > 0)
                refcall_mask |= ref_from_ad
        except:
            pass
    
    # Remove ALL RefCall variants
    df = df[~refcall_mask]
    refcall_removed = initial_count - len(df)
    
    print(f"Removed {refcall_removed} reference call variants")
    print(f"Variants remaining: {len(df)}")
    
    return df


def remove_extreme_low_quality_variants(df):
    """
    Remove variants with extremely low quality (DP <= 5 or GQ <= 10)
    These variants are completely removed from results
    
    Args:
        df: DataFrame with variant data
        
    Returns:
        df: DataFrame with extreme low quality variants removed
    """
    initial_count = len(df)
    
    # Create mask for extreme low quality variants
    extreme_mask = pd.Series([False] * len(df), index=df.index)
    
    # Check for extreme low DP
    if 'CoverageDepth' in df.columns:
        df['CoverageDepth'] = pd.to_numeric(df['CoverageDepth'], errors='coerce')
        extreme_dp_mask = df['CoverageDepth'] <= EXTREME_DP_THRESHOLD
        extreme_mask |= extreme_dp_mask
        
        extreme_dp_count = extreme_dp_mask.sum()
        if extreme_dp_count > 0:
            print(f"  Found {extreme_dp_count} variants with DP <= {EXTREME_DP_THRESHOLD}")
    
    # Check for extreme low GQ
    if 'GenotypeQual' in df.columns:
        df['GenotypeQual'] = pd.to_numeric(df['GenotypeQual'], errors='coerce')
        extreme_gq_mask = df['GenotypeQual'] <= EXTREME_GQ_THRESHOLD
        extreme_mask |= extreme_gq_mask
        
        extreme_gq_count = extreme_gq_mask.sum()
        if extreme_gq_count > 0:
            print(f"  Found {extreme_gq_count} variants with GQ <= {EXTREME_GQ_THRESHOLD}")
    
    # Remove extreme low quality variants
    df = df[~extreme_mask]
    extreme_removed = initial_count - len(df)
    
    if extreme_removed > 0:
        print(f"  Removed {extreme_removed} variants with extreme low quality (DP <= {EXTREME_DP_THRESHOLD} or GQ <= {EXTREME_GQ_THRESHOLD})")
    
    print(f"  Variants remaining after extreme quality filter: {len(df)}")
    
    return df


def apply_quality_filters(df):
    """
    Apply quality filters to remaining variants
    Variants that fail these filters get flags but remain in results
    
    Priority order:
    1. DP - insufficient coverage depth (most stringent)
    2. GQ - insufficient genotype quality 
    3. AlleleBalance - suspicious allele balance (lowest priority)
    
    If multiple filters apply, only the most stringent is shown in GATK_FILTER
    
    Args:
        df: DataFrame with variant data
        
    Returns:
        df: DataFrame with GATK_FILTER column updated according to priority
    """
    
    # Initialize all variants as PASS
    if 'GATK_FILTER' not in df.columns:
        df['GATK_FILTER'] = 'PASS'
    
    print("\nApplying quality filters (variants remain with flags)...")
    
    # Priority 1: Coverage Depth filter - most stringent technical filter
    if 'CoverageDepth' in df.columns:
        print("1. Applying CoverageDepth filter (DP < 10)...")
        
        df['CoverageDepth'] = pd.to_numeric(df['CoverageDepth'], errors='coerce')
        dp_mask = (df['CoverageDepth'] < 10) & (df['CoverageDepth'] > EXTREME_DP_THRESHOLD)
        
        # Apply DP filter (most stringent)
        df.loc[dp_mask, 'GATK_FILTER'] = 'DP'
        
        dp_filtered = len(df[df['GATK_FILTER'] == 'DP'])
        print(f"   DP filter: {dp_filtered} variants filtered (DP < 10)")
    else:
        print("1. CoverageDepth filter: column not found, skipping")

    # Priority 2: Genotype Quality filter 
    if 'GenotypeQual' in df.columns:
        print("2. Applying GenotypeQual filter (GQ < 20)...")
        
        df['GenotypeQual'] = pd.to_numeric(df['GenotypeQual'], errors='coerce')
        gq_mask = (df['GenotypeQual'] < 20) & (df['GenotypeQual'] > EXTREME_GQ_THRESHOLD) & (df['GATK_FILTER'] == 'PASS')
        
        # Apply GQ filter only to PASS variants
        df.loc[gq_mask, 'GATK_FILTER'] = 'GQ'
        
        gq_filtered = len(df[df['GATK_FILTER'] == 'GQ'])
        print(f"   GQ filter: {gq_filtered} variants filtered (GQ < 20)")
    else:
        print("2. GenotypeQual filter: column not found, skipping")
    
    # Priority 3: AlleleBalance Filter - lowest priority
    if 'AlleleDepth' in df.columns:
        print("3. Applying AlleleBalance filter (lowest priority)...")
        
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

        # Identify variants with suspicious AlleleBalance
        suspicious_ab_mask = ~(
            (df['GT_inferred'] == 'hom_ref') |  # AB <= 0.1
            (df['GT_inferred'] == 'het') |      # AB 0.25-0.75  
            (df['GT_inferred'] == 'hom_alt')    # AB >= 0.9
        )
        
        # Apply AlleleBalance filter only to PASS variants (lowest priority)
        ab_apply_mask = suspicious_ab_mask & (df['GATK_FILTER'] == 'PASS')
        df.loc[ab_apply_mask, 'GATK_FILTER'] = 'AlleleBalance'

        ab_filtered = len(df[df['GATK_FILTER'] == 'AlleleBalance'])
        print(f"   AlleleBalance filter: {ab_filtered} variants filtered")
    else:
        print("3. AlleleBalance filter: column not found, skipping")
    
    return df


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
        print("No ranking file provided - skipping ranking annotation")
    
    # Clean any problematic characters in string columns
    str_columns = df.select_dtypes(include=['object']).columns
    for col in str_columns:
        df[col] = df[col].apply(clean_text)
    
    print("\n" + "="*60)
    print("INITIAL DATA")
    print("="*60)
    print(f"Total variants loaded: {len(df)}")
    
    # Step 1: Remove ALL RefCall variants
    print("\n" + "="*60)
    print("STEP 1: REMOVING ALL REFERENCE CALL VARIANTS")
    print("="*60)
    df = remove_refcall_variants(df)
    
    # Step 2: Remove extreme low quality variants (DP <= 5 or GQ <= 10)
    print("\n" + "="*60)
    print(f"STEP 2: REMOVING EXTREME LOW QUALITY VARIANTS")
    print(f"(DP <= {EXTREME_DP_THRESHOLD} or GQ <= {EXTREME_GQ_THRESHOLD})")
    print("="*60)
    df = remove_extreme_low_quality_variants(df)
    
    # Step 3: Apply consequence filter if enabled
    if consequence_filter_enabled:
        print("\n" + "="*60)
        print("STEP 3: APPLYING CONSEQUENCE FILTER")
        print("="*60)
        before_consequence = len(df)
        
        # Consequences to keep (from missense_variant to transcript_ablation)
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
        
        # Keep only variants with consequences in our list
        df = df[df['Consequence'].isin(consequences_to_keep)]
        
        consequence_removed = before_consequence - len(df)
        print(f"Removed {consequence_removed} non-coding variants")
        print(f"Variants remaining: {len(df)}")

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
            print(f"\nMerging ranking data using column: {gene_symbol_col}")
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

    # Step 4: Apply quality filters (variants remain with flags)
    print("\n" + "="*60)
    print("STEP 4: APPLYING QUALITY FILTERS")
    print("(variants remain in results with flags)")
    print("="*60)
    df = apply_quality_filters(df)

    # Step 5: Apply gnomAD filter (remove variants without adding flag)
    if gnomad_filter_enabled:
        print("\n" + "="*60)
        print(f"STEP 5: APPLYING GNOMAD FILTER (AF >= {gnomad_af_threshold})")
        print("="*60)
        before_gnomad = len(df)
        
        # Check for gnomAD columns
        gnomad_cols = [col for col in ['gnomAD_exome_NFE', 'gnomAD_genome_NFE'] if col in df.columns]
        
        if gnomad_cols:
            print(f"Found gnomAD columns: {', '.join(gnomad_cols)}")
            
            # Create mask for variants with gnomAD frequency >= threshold
            gnomad_mask = pd.Series([False] * len(df), index=df.index)
            
            for col in gnomad_cols:
                # Filter variants with frequency >= threshold
                col_mask = df[col] >= gnomad_af_threshold
                gnomad_mask = gnomad_mask | col_mask
                
                filtered_by_col = col_mask.sum()
                if filtered_by_col > 0:
                    print(f"{col}: {filtered_by_col} variants with frequency >= {gnomad_af_threshold}")
            
            # Remove gnomAD filtered variants
            df = df[~gnomad_mask]
            gnomad_removed = before_gnomad - len(df)
            print(f"\nRemoved {gnomad_removed} gnomAD variants (AF >= {gnomad_af_threshold})")
            print(f"Variants remaining: {len(df)}")
        else:
            print("gnomAD filter: no gnomAD columns found in data")

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
        df['Consequence_severity'] = df['Consequence'].map(consequence_order)
        df = df.sort_values('Consequence_severity', ascending=False).drop_duplicates(subset='Chr', keep='first')

    # Print final statistics
    print("\n" + "="*60)
    print("FINAL STATISTICS")
    print("="*60)
    
    # Count variants by filter type
    filter_counts = df['GATK_FILTER'].value_counts()
    
    print("\nFilter distribution (variants in final results):")
    for filter_type in ['DP', 'GQ', 'AlleleBalance', 'PASS']:
        if filter_type in filter_counts:
            print(f"  {filter_type}: {filter_counts[filter_type]} variants")
    
    # Print any additional filters not in the standard list
    other_filters = set(filter_counts.keys()) - set(['DP', 'GQ', 'AlleleBalance', 'PASS'])
    for filt in sorted(other_filters):
        print(f"  {filt}: {filter_counts[filt]} variants")
    
    pass_count = filter_counts.get('PASS', 0)
    filtered_count = len(df) - pass_count
    
    print(f"\nTotal variants in output: {len(df)}")
    print(f"PASS variants: {pass_count}")
    print(f"Filtered variants (with flags): {filtered_count}")
    
    # Summary of removed variants
    print("\n" + "="*60)
    print("REMOVED VARIANTS SUMMARY")
    print("="*60)
    print(f"• Reference calls: Removed completely")
    print(f"• Extreme low quality (DP <= {EXTREME_DP_THRESHOLD} or GQ <= {EXTREME_GQ_THRESHOLD}): Removed completely")
    if consequence_filter_enabled:
        print(f"• Non-coding variants: Removed by consequence filter")
    if gnomad_filter_enabled:
        print(f"• gnomAD variants (AF >= {gnomad_af_threshold}): Removed completely")

    # Final check: Ensure NO RefCall variants remain
    print("\n" + "="*60)
    print("FINAL CHECK: VERIFYING NO REFERENCE CALLS REMAIN")
    print("="*60)
    
    refcall_found = 0
    for col in df.columns:
        col_str = str(col).upper()
        if any(keyword in col_str for keyword in ['ZYG', 'GT', 'GENOTYPE', 'CALL']):
            col_values = df[col].astype(str)
            ref_in_col = (
                col_values.str.contains('^0/0$').sum() +
                col_values.str.contains('^0\|0$').sum() +
                col_values.str.contains('^REF$', case=False).sum()
            )
            if ref_in_col > 0:
                print(f"WARNING: Found {ref_in_col} reference calls in column '{col}'!")
                refcall_found += ref_in_col
    
    if refcall_found == 0:
        print("SUCCESS: No reference call variants found in final results!")
    else:
        print(f"ERROR: Found {refcall_found} reference calls in final results!")
        print("This should not happen. Please check the data.")

    # Print ranking statistics
    if 'rank_gse71613' in df.columns:
        ranked_count = len(df[df['rank_gse71613'].notna()])
        if ranked_count > 0:
            print(f"\nRanking annotation: {ranked_count} variants have ranking information")
            print(f"Ranking range: {df['rank_gse71613'].min():.1f} - {df['rank_gse71613'].max():.1f}")

    # Clean up temporary columns
    cols_to_drop = ['AD1', 'AD2', 'AlleleBalance', 'GT_inferred', 'Consequence_severity']
    if 'HGNC symbol' in df.columns:
        cols_to_drop.extend(['HGNC symbol', 'NCBI gene (formerly Entrezgene) ID', 
                           'Gene name', 'BIOTYPE', 'CANONICAL', 'Gene stable ID'])
    
    # Only drop columns that exist
    cols_to_drop = [col for col in cols_to_drop if col in df.columns]
    df.drop(cols_to_drop, axis=1, inplace=True)

    # Reorder columns: move CoverageDepth, GenotypeQual, AlleleDepth after rsID
    if 'rsID' in df.columns:
        # Get current column order
        columns = list(df.columns)
        
        # Columns to move
        cols_to_move = ['CoverageDepth', 'GenotypeQual', 'AlleleDepth']
        
        # Find rsID position
        rsid_index = columns.index('rsID') if 'rsID' in columns else -1
        
        if rsid_index != -1:
            # Remove columns from their current positions
            for col in cols_to_move:
                if col in columns:
                    columns.remove(col)
            
            # Insert columns after rsID
            insert_position = rsid_index + 1
            for col in reversed(cols_to_move):  # Reverse to maintain order
                if col in df.columns:
                    columns.insert(insert_position, col)
            
            # Reorder dataframe
            df = df[columns]
    
    # Save filtered results
    df.to_csv(output_file, index=False, sep='\t', na_rep='.')
    print(f"\nResults saved to: {output_file}")
    print("="*60)

if __name__ == "__main__":
    with open(log_file, "w") as f:
        sys.stderr = sys.stdout = f
        main()