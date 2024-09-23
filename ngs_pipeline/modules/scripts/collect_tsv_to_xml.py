import pandas as pd
import os


input_files = snakemake.input.tsv
output_file = snakemake.output.xlsx
params_file = snakemake.params.tsv

def merge_tsv(input_files, output_file):
	with pd.ExcelWriter(output_file) as writer:
		for file in input_files:
			df = pd.read_csv(file, sep='\t')
			sample_short = os.path.basename(file).split('.')[0]
			sheet_name= params_file.loc[params_file['sample'] == sample_short, 'full_name'].values[0]
			df.to_excel(writer, sheet_name=sheet_name, index=False, float_format='%.8f')

params_file = pd.read_csv(params_file, sep='\t')
merge_tsv(input_files, output_file)
