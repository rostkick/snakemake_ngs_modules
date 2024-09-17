import os
from modules.scripts.params_builder import *
from modules.scripts.get_input import get_final_inputs


configfile: 'data/configure.run_settings.yml'
configfile: 'configure.solid_deps.yml' # will be merged with current settings configuration from 'configure.run_settings.yml'

if config['assembly'] == 'GRCh38':
	configfile: 'configure.references38.yml'
elif config['assembly'] == 'GRCh37':
	configfile: 'configure.references37.yml'

ngs = NGSSetup()
data = ngs.data

if not os.path.exists(f"results/{config['run']}"):
    os.makedirs(f"results/{config['run']}")

data.to_csv(f"results/{config['run']}/run_table.tsv", sep='\t', index=False)

wildcard_constraints:
	sample="|".join(ngs.SAMPLES),
	patient = "|".join(ngs.TMR_PATIENTS)

# === for debugging ===
# print(ngs.__dict__)
# print(f'{ngs.GRM_SAMPLES=}')
# print(f'{ngs.SAMPLES=}')
# print(f'{ngs.GRM_VS_TMR_PATIENTS=}')
# print(f'{ngs.TMR_SAMPLES=}')
# print(f'{ngs.TMR_PATIENTS=}')
# print(f'{ngs.ONLY_TMR_PATIENTS=}')
# === for debugging ===

rule all:
	input: get_final_inputs(ngs)

include: config["snakemake_modules"] + "rules_1.aligning.smk"
include: config["snakemake_modules"] + "rules_2.preprocessing.smk"
include: config["snakemake_modules"] + "rules_3.collect_metrics.smk"
include: config["snakemake_modules"] + "rules_4.germline_calling.smk"
include: config["snakemake_modules"] + "rules_5.somatic_calling.initial.smk"
include: config["snakemake_modules"] + "rules_6.somatic_calling.paired.smk"
include: config["snakemake_modules"] + "rules_7.somatic_calling.tmr_only.smk"
include: config["snakemake_modules"] + "rules_8.annotation.smk"
include: config["snakemake_modules"] + "rules_9.parse_results.smk"
