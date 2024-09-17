from itertools import product
from snakemake.workflow import config


def get_final_inputs(ngs):
	"""Return a list of input files for all samples."""
	germline_inputs = []
	somatic_inputs = []
	metrics = []

	if ngs.GRM:
		germline_inputs_cohort = [
			f"results/{config['run']}/germline/vcf/cohort.annotated.vcf.gz"
		]
		germline_inputs_individual = [
			f"results/{config['run']}/germline/xlsx/individual.germline.results.xlsx"
		]
		germline_inputs += germline_inputs_cohort + germline_inputs_individual
		if ngs.TMR:
			somatic_inputs += [
				f"results/{config['run']}/somatic/{patient}/somatic_annotated.vcf.gz"
				for patient in ngs.GRM_VS_TMR_PATIENTS
			]

	if len(ngs.ONLY_TMR_PATIENTS) > 0:
		somatic_inputs += [
			f"results/{config['run']}/somatic/{patient}/somatic_annotated_tonly.vcf.gz"
			for patient in ngs.ONLY_TMR_PATIENTS
		]

	if config['ngs_type'] == 'WES':
		metrics = [
			f"results/{config['run']}/bam/hs_metrics/{sample}.hs_metrics.tsv"
			for sample in ngs.SAMPLES
		]

	input_files = (
		germline_inputs + somatic_inputs + metrics
	)

	return input_files
