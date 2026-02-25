"""Convert a single TSV file to an xlsx workbook with one data sheet.
Used for per-sample full WES results (no Legend sheet needed).

Snakemake script interface:
    snakemake.input.tsv   - input TSV file
    snakemake.output.xlsx - output xlsx file
    snakemake.log[0]      - log file
"""
import sys
import os
import pandas as pd

log_file   = snakemake.log[0]
input_file = snakemake.input.tsv
output_file = snakemake.output.xlsx


def tsv_to_xlsx(input_file, output_file):
    sample_name = os.path.basename(input_file).split('.')[0]
    df = pd.read_csv(input_file, sep='\t')
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name=sample_name, index=False, float_format='%.8f')
        # Auto-width columns
        ws = writer.sheets[sample_name]
        for col in ws.columns:
            max_len = max(len(str(cell.value)) if cell.value is not None else 0 for cell in col)
            ws.column_dimensions[col[0].column_letter].width = min(max_len + 2, 60)


with open(log_file, 'w') as f:
    sys.stderr = sys.stdout = f
    tsv_to_xlsx(input_file, output_file)
