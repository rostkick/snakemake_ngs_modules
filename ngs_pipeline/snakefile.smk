from modules.scripts.params_builder import *
from modules.scripts.get_input import final_inputs


configfile: 'configure.yml'

ngs = NGSSetup()
data = ngs.data

ngs.data.to_csv('run_table.tsv', sep='\t', index=False)

wildcard_constraints:
	sample="|".join(ngs.SAMPLES),
	patient = "|".join(ngs.TMR_PATIENTS)

rule all:
	input: final_inputs(ngs)

include: config["snakemake_modules"] + "rules_1.aligning.smk"
include: config["snakemake_modules"] + "rules_2.preprocessing.smk"
include: config["snakemake_modules"] + "rules_3.collect_metrics.smk"
include: config["snakemake_modules"] + "rules_4.germline_calling.haplotypecaller.smk"
include: config["snakemake_modules"] + "rules_5.somatic_calling.smk"
include: config["snakemake_modules"] + "rules_6.annotation.smk"
