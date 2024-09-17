import pandas as pd
import os


input_files = snakemake.input.tsv
output_file = snakemake.output.xlsx

def merge_tsv(input_files, output_file):
	with pd.ExcelWriter(output_file) as writer:
		for file in input_files:
			df = pd.read_csv(file, sep='\t')
			sheet_name = os.path.basename(file).split('.')[0]
			df.to_excel(writer, sheet_name=sheet_name, index=False)

merge_tsv(input_files, output_file)
