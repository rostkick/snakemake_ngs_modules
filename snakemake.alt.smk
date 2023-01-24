import sys


from modules.scripts.functions import *
from modules.scripts.params_builder import RUN, GERMLINE, SOMATIC, NGSData


ngs = NGSData(config['grm_dir'], config['tmr_dir'])

ONLY_TUMOR = ngs.wide_df[ngs.wide_df['grm_samples'].isnull()]['patients'].to_list()
PATIENTS = [i for i in ngs.wide_df.loc[:, "patients"].tolist() if i not in ONLY_TUMOR]

SAMPLES = ngs.long_df.loc[:, "samples"].tolist()
GERMLINE_SAMPLES = ngs.wide_df.loc[:, "grm_samples"].dropna().tolist()
SOMATIC_SAMPLES = ngs.wide_df.loc[:, "tmr_samples"].dropna().tolist()

wildcard_constraints:
	sample="|".join(SAMPLES)

def final_inputs():
	# germline
	germline_inputs = [f'results/{RUN}/germline/vcf/germline.annotated.vcf.gz'] + \
					[f'results/{RUN}/germline/vcf/{sample}.germline.annotated.vcf.gz' for sample in SAMPLES]
	# somatic
	somatic_inputs = [f'results/{RUN}/somatic/{patient}/annotation/somatic.annotated.vcf.gz' for patient in PATIENTS]
	# metrics
	metrics = [f"results/{RUN}/bam/hs_metrics/{sample}.hs_metrics.tsv" for sample in SAMPLES]
	if GERMLINE and SOMATIC is False:
		return germline_inputs + metrics
	elif GERMLINE is False and SOMATIC:
		return somatic_inputs + metrics
	elif GERMLINE and SOMATIC:
		return germline_inputs + somatic_inputs + metrics

# print(ngs.wide_df)
ngs.wide_df.to_csv('wide_param.tsv', sep='\t', index=None)
ngs.long_df.to_csv('long_param.tsv', sep='\t', index=None)
# print(ngs.long_df)

rule all:
	input: final_inputs()

include: config["snakemake_modules"] + "aligning.smk"
include: config["snakemake_modules"] + "preprocessing.smk"
include: config["snakemake_modules"] + "collect_metrics.smk"
include: config["snakemake_modules"] + "germline_calling.deepvariant.smk"
include: config["snakemake_modules"] + "somatic_calling.smk"
include: config["snakemake_modules"] + "sv_calling.smk"
include: config["snakemake_modules"] + "liftover.smk"
include: config["snakemake_modules"] + "annotation.smk"
