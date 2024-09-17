import pandas as pd


input_file = snakemake.input.tsv
output_file = snakemake.output.tsv

# Get only protein coding
df = pd.read_csv(input_file, sep='\t')
df = df[df['BIOTYPE'] == 'protein_coding']

# DP/GQ Filter
df['CoverageDepth'] = pd.to_numeric(df['CoverageDepth'])
df['GenotypeQual'] = pd.to_numeric(df['GenotypeQual'])
df.loc[df['CoverageDepth'] < 10, 'GATK_FILTER'] = 'DP'
df.loc[df['GenotypeQual'] < 20, 'GATK_FILTER'] = 'GQ'

# AlleleBalance Filter
df[['AD1', 'AD2']] = df['AlleleDepth'].str.split(',', expand=True)
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
	'intergenic': 0,
	'feature_truncation': 1,
	'feature_elongation': 2,
	'regulatory': 3,
	'TF_binding_site': 4,
	'TFBS': 5,
	'downstream': 6,
	'upstream': 7,
	'non_coding_transcript': 8,
	'non_coding': 9,
	'intron': 10,
	'NMD_transcript': 11,
	'non_coding_transcript_exon': 12,
	'5_prime_utr': 13,
	'3_prime_utr': 14,
	'coding_sequence': 15,
	'mature_miRNA': 16,
	'stop_retained': 17,
	'start_retained': 18,
	'synonymous': 19,
	'incomplete_terminal_codon': 20,
	'splice_region': 21,
	'missense': 22,
	'inframe': 23,
	'protein_altering': 24,
	'transcript_amplification': 25,
	'exon_loss': 26,
	'disruptive': 27,
	'start_lost': 28,
	'stop_lost': 29,
	'stop_gained': 30,
	'frameshift': 31,
	'splice_acceptor': 32,
	'splice_donor': 33,
	'transcript_ablation': 34
}

df = df.drop_duplicates(subset='Chr', keep='last')
df['Consequence_severity'] = df['Consequence'].map(consequence_order)
df = df.sort_values('Consequence_severity', ascending=False).drop_duplicates(subset='Chr', keep='first')

# Drop unwanted columns from resulting table
df.drop(['AD1', 'AD2', 'AlleleBalance', 'Consequence_severity'], axis=1, inplace=True)

df.to_csv(output_file, index=False, sep='\t')
