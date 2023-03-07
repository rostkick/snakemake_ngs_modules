from itertools import product
from snakemake.workflow import config


def final_inputs(ngs):
	germline_inputs, somatic_inputs, metrics = [], [], []

	# germline
	germline_inputs = [f'results/{run}/germline/vcf/{sample}.annotated.vcf.gz' for run, sample in product([config['run']], ngs.GRM_SAMPLES)]
	if len(ngs.GRM_SAMPLES) > 1:
		germline_inputs = germline_inputs + [f'results/{run}/germline/vcf/cohort.annotated.vcf.gz' for run in [config['run']]]
	# somatic
	somatic_inputs = [f'results/{run}/somatic/{patient}/annotated.vcf.gz' for run, patient in product([config['run']], ngs.TMR_PATIENTS)]
	# metrics
	metrics = [f"results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv" for run, sample in product([config['run']], ngs.SAMPLES)]

	if ngs.GRM and ngs.TMR is False:
		return germline_inputs + metrics
	elif ngs.GRM is False and ngs.TMR:
		return somatic_inputs + metrics
	elif ngs.GRM and ngs.TMR:
		return germline_inputs + somatic_inputs + metrics
