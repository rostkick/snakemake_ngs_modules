from itertools import product
from modules.scripts.params_builder import *
from pprint import pprint

configfile: 'configure.yml'


ngs = NGSSetup()
data = ngs.data
mapping = ngs.mapping

def final_inputs():
	germline_inputs, somatic_inputs, metrics = [], [], []

	# germline
	germline_inputs = [f'results/{run}/germline/vcf/{sample}.annotated.vcf.gz' for run, sample in product([config['run']], ngs.GRM_SAMPLES)]
	if len(ngs.GRM_SAMPLES) > 1:
		germline_inputs = germline_inputs + [f'results/{run}/germline/vcf/cohort.annotated.vcf.gz' for run in [config['run']]]
	# somatic
	somatic_inputs = [f'results/{run}/somatic/{patient}/annotated.vcf.gz' for run, patient in product([config['run']], ngs.ALL_PATIENTS)]
	# metrics
	metrics = [f"results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv" for run, sample in product([config['run']], ngs.SAMPLES)]

	if ngs.GRM and ngs.TMR is False:
		return germline_inputs + metrics
	elif ngs.GRM is False and ngs.TMR:
		return somatic_inputs + metrics
	elif ngs.GRM and ngs.TMR:
		return germline_inputs + somatic_inputs + metrics

wildcard_constraints:
	sample="|".join(ngs.SAMPLES),
	patient = "|".join(ngs.ALL_PATIENTS)

pprint(ngs.__dict__)
rule all:
	# input: final_inputs()
	# input:  [f'results/{run}/germline/vcf/{sample}.annotated.vcf.gz' for run, sample in product([config['run']], ngs.GRM_SAMPLES)]
	input: expand('results/{run}/somatic/{patient}/final.vcf.gz', run=config['run'], patient=ngs.ONLY_TMR_PATIENTS),
	# 		 [f'results/{run}/germline/vcf/cohort.annotated.vcf.gz' for run in [config['run']]]
	# input: 

include: config["snakemake_modules"] + "rules_1.aligning.smk"
include: config["snakemake_modules"] + "rules_2.preprocessing.smk"
include: config["snakemake_modules"] + "rules_3.collect_metrics.smk"
include: config["snakemake_modules"] + "rules_4.germline_calling.deepvariant.smk"
include: config["snakemake_modules"] + "rules_5.somatic_calling.smk"
include: config["snakemake_modules"] + "rules_6.somatic_calling.grm_vs_tmr.smk"
include: config["snakemake_modules"] + "rules_7.somatic_calling.tmr_only.smk"
# include: config["snakemake_modules"] + "rules_8.sv_calling.smk"
# include: config["snakemake_modules"] + "rules_9.annotation.smk"
