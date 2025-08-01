import pandas as pd
import os
import sys
from openpyxl.worksheet.hyperlink import Hyperlink

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
        
        # First create all data sheets
        for data in sheet_data:
            df = pd.read_csv(data['file'], sep='\t')
            df.to_excel(writer, sheet_name=data['full_name'], index=False, float_format='%.8f')
        
        # Now add hyperlinks to Legend sheet
        workbook = writer.book
        legend_ws = workbook['Legend']
        
        # Format headers
        for cell in legend_ws[1]:
            cell.font = cell.font.copy(bold=True)
        
        # Add hyperlinks to sheet names using Hyperlink object with text override
        for i, data in enumerate(sheet_data, start=2):
            sheet_name = data['full_name']
            cell = legend_ws.cell(row=i, column=1)
            
            # Create hyperlink object that works for navigation
            hyperlink_obj = Hyperlink(ref="", location=f"'{sheet_name}'!A1")
            # Override the display text to show the sheet name instead of gid reference
            hyperlink_obj.display = sheet_name
            
            # Apply the hyperlink to the cell
            cell.hyperlink = hyperlink_obj
            cell.value = sheet_name  # Ensure cell value is set
            # Style the cell to look like a hyperlink
            cell.font = cell.font.copy(color="0000FF", underline="single")

with open(log_file, "w") as f:
    sys.stderr = sys.stdout = f
    params_df = pd.read_csv(params_file, sep='\t')
    merge_tsv_with_legend(input_files, output_file, params_df)