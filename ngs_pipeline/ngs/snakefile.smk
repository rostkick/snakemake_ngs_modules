from modules.scripts.params_builder import *
from modules.scripts.get_input import final_inputs


# set configuration
configfile: 'pipeline/configure_run.yml'
configfile: 'pipeline/configure_tools.yml'

if config['assembly'] == 'GRCh37':
	configfile: 'pipeline/configure_reference_37.yml'
elif config['assembly'] == 'GRCh38':
	configfile: 'pipeline/configure_reference_38.yml'

# init run env & paths
ngs = NGSSetup()

ngs.data.to_csv('results/bla.tsv', sep='\t', index=False)

wildcard_constraints:
	sample="|".join(ngs.SAMPLES)

rule all:
	input: [f'results/{run}/germline/vcf/cohort.annotated.vcf.gz' for run in [config['run']]]

include: config["snakemake_modules"] + "rules_1.aligning.smk"
include: config["snakemake_modules"] + "rules_2.preprocessing.smk"
include: config["snakemake_modules"] + "rules_3.germline_calling.haplotypecaller.smk"
include: config["snakemake_modules"] + "rules_4.annotation.smk"
