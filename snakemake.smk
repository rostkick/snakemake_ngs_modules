import sys
from itertools import product

from modules.scripts.params_builder import NGSSetup


ngs_setup = NGSSetup()
ngs = ngs_setup.create_params()

ngs.wide_df.to_csv('wide_params.tsv', sep='\t')
ngs.long_df.to_csv('long_params.tsv', sep='\t')

ALL_PATIENTS = ngs.wide_df[~ngs.wide_df['tmr_samples'].isnull()]['patients'].to_list()
ONLY_TMR_PATIENTS = ngs.wide_df[ngs.wide_df['grm_samples'].isnull()]['patients'].to_list()
GRM_VS_TMR_PATIENTS = ngs.wide_df.query("~(tmr_samples.isnull() | grm_samples.isnull())")['patients'].to_list()
SAMPLES = ngs.long_df.loc[:, "samples"].tolist()
GRM_SAMPLES = ngs.wide_df.loc[:, "grm_samples"].dropna().tolist()
TMR_SAMPLES = ngs.wide_df.loc[:, "tmr_samples"].dropna().tolist()


wildcard_constraints:
	sample="|".join(SAMPLES),
	patient="|".join(ALL_PATIENTS)

def final_inputs():
	# germline
	germline_inputs = [f'results/{run}/germline/vcf/cohort.annotated.vcf.gz' for run in [config['run']]] + \
					[f'results/{run}/germline/vcf/{sample}.annotated.vcf.gz' for run, sample in product([config['run']], GRM_SAMPLES)]
	# somatic
	somatic_inputs = [f'results/{run}/somatic/{patient}/annotated.vcf.gz' for run, patient in product([config['run']], ALL_PATIENTS)]
	# metrics
	metrics = [f"results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv" for run, sample in product([config['run']], SAMPLES)]
	if ngs_setup.GRM and ngs_setup.TMR is False:
		return germline_inputs + metrics
	elif ngs_setup.GRM is False and ngs_setup.TMR:
		return somatic_inputs + metrics
	elif ngs_setup.GRM and ngs_setup.TMR:
		return germline_inputs + somatic_inputs + metrics


rule all:
	input: final_inputs()

include: config["snakemake_modules"] + "aligning.smk"
include: config["snakemake_modules"] + "preprocessing.smk"
include: config["snakemake_modules"] + "collect_metrics.smk"
include: config["snakemake_modules"] + "germline_calling.deepvariant.smk"
include: config["snakemake_modules"] + "somatic_calling.smk"
include: config["snakemake_modules"] + "annotation.smk"
