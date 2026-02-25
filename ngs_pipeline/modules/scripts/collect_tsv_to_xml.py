import pandas as pd
import os
import sys
from openpyxl.worksheet.hyperlink import Hyperlink
from openpyxl.styles import Font

log_file = snakemake.log[0]
input_files = snakemake.input.tsv
output_file = snakemake.output.xlsx
params_file = snakemake.params.tsv

def merge_tsv_with_legend(input_files, output_file, params_df):
    # Collect sheet names for legend
    sheet_data = []
    
    # First pass: collect all sheet information
    for file in input_files:
        sample_short = os.path.basename(file).split('.')[0]
        full_name = params_df.loc[params_df['sample'] == sample_short, 'full_name'].values[0]
        sheet_data.append({
            'file': file,
            'sample_short': sample_short,
            'full_name': full_name
        })
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Create legend DataFrame
        legend_data = []
        for data in sheet_data:
            legend_data.append([data['full_name'], data['sample_short']])
        
        legend_df = pd.DataFrame(legend_data, columns=['Sheet Name', 'Sample'])
        legend_df.to_excel(writer, sheet_name='Legend', index=False)
        
        # Create all data sheets with startrow=1 to leave row 1 for back-link
        for data in sheet_data:
            df = pd.read_csv(data['file'], sep='\t')
            df.to_excel(writer, sheet_name=data['full_name'], index=False,
                        float_format='%.8f', startrow=1)
        
        workbook = writer.book

        # --- Legend sheet: format headers + add hyperlinks to sample sheets ---
        legend_ws = workbook['Legend']
        
        for cell in legend_ws[1]:
            cell.font = cell.font.copy(bold=True)
        
        for i, data in enumerate(sheet_data, start=2):
            sheet_name = data['full_name']
            cell = legend_ws.cell(row=i, column=1)
            hyperlink_obj = Hyperlink(ref="", location=f"'{sheet_name}'!A1")
            hyperlink_obj.display = sheet_name
            cell.hyperlink = hyperlink_obj
            cell.value = sheet_name
            cell.font = cell.font.copy(color="0000FF", underline="single")

        # --- Sample sheets: add back-link to Legend in A1 ---
        for data in sheet_data:
            sheet_name = data['full_name']
            ws = workbook[sheet_name]
            back_cell = ws.cell(row=1, column=1)
            back_cell.value = "← Legend"
            hyperlink_obj = Hyperlink(ref="", location="'Legend'!A1")
            hyperlink_obj.display = "← Legend"
            back_cell.hyperlink = hyperlink_obj
            back_cell.font = Font(color="0000FF", underline="single", bold=True)

with open(log_file, "w") as f:
    sys.stderr = sys.stdout = f
    params_df = pd.read_csv(params_file, sep='\t')
    merge_tsv_with_legend(input_files, output_file, params_df)