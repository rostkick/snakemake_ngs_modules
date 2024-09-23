import pandas as pd


input_file = snakemake.input.tsv
mart_file = snakemake.params.mart
output_file = snakemake.output.tsv


df = pd.read_csv(input_file, sep='\t')
df_mart = pd.read_csv(mart_file, sep='\t')

# Get only protein coding
df = df[df['BIOTYPE'] == 'protein_coding']

gnomad_cols = [col for col in df.columns if 'gnomAD' in col]
df[gnomad_cols] = df[gnomad_cols].replace('.', 0)
df[gnomad_cols] = df[gnomad_cols].astype(float)

df = pd.merge(df, df_mart, left_on='RefGene', right_on='HGNC symbol', how='left')

# DP/GQ Filter
df['CoverageDepth'] = pd.to_numeric(df['CoverageDepth'])
df['GenotypeQual'] = pd.to_numeric(df['GenotypeQual'])
df.loc[df['CoverageDepth'] < 10, 'GATK_FILTER'] = 'DP'
df.loc[df['GenotypeQual'] < 20, 'GATK_FILTER'] = 'GQ'

# AlleleBalance Filter
df[['AD1', 'AD2']] = df['AlleleDepth'].str.split(',', expand=True).iloc[:,:2]
df['AD1'] = pd.to_numeric(df['AD1'])
df['AD2'] = pd.to_numeric(df['AD2'])
df['AlleleBalance'] = df['AD1'] / (df['AD1'] + df['AD2'])
df.loc[(df['AlleleBalance'] >= 0.25) & (df['AlleleBalance'] <= 0.75), 'Zyg'] = 'HET'
df.loc[(df['AlleleBalance'] >= 0.1) & (df['AlleleBalance'] < 0.25) | (df['AlleleBalance'] > 0.75) & (df['AlleleBalance'] <= 0.9), 'Zyg'] = 'HOM'
df.loc[(df['AlleleBalance'] < 0.1) | (df['AlleleBalance'] > 0.9), ['Zyg', 'GATK_FILTER']] = ['.', 'AlleleBalance']

# Apply Filters
df = df[df['GATK_FILTER'] == 'PASS']

# Consequence priority in case of multiple canonical transcripts
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

df = df.drop_duplicates(subset='Chr', keep='last')
df['Consequence_severity'] = df['Consequence'].map(consequence_order)
df = df.sort_values('Consequence_severity', ascending=False).drop_duplicates(subset='Chr', keep='first')

# Drop unwanted columns from resulting table
df.drop(['AD1', 'AD2', 'AlleleBalance', 'Consequence_severity',
		 'HGNC symbol', 'NCBI gene (formerly Entrezgene) ID',
		 'Gene name', 'BIOTYPE', 'CANONICAL', 'Gene stable ID',
		 'GATK_FILTER'], axis=1, inplace=True)

df.to_csv(output_file, index=False, sep='\t')
