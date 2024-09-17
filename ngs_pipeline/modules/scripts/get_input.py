from itertools import product
from snakemake.workflow import config


def final_inputs(ngs):
	germline_inputs, somatic_inputs_grm_vs_tmr, somatic_inputs_tonly, metrics = [], [], [], []

	if ngs.GRM:
		germline_inputs_cohort = [f'results/{run}/germline/vcf/cohort.annotated.vcf.gz' for run in [config['run']]]
		germline_inputs_individual = [f'results/{run}/germline/xlsx/individual.germline.results.xlsx' for run in [config['run']]]
		germline_inputs = germline_inputs + germline_inputs_cohort + germline_inputs_individual
		# somatic
		if ngs.TMR:
			somatic_inputs_grm_vs_tmr = [f'results/{run}/somatic/{patient}/somatic_annotated.vcf.gz' for run, patient in product([config['run']]\
										, ngs.GRM_VS_TMR_PATIENTS)]
	
	if len(ngs.ONLY_TMR_PATIENTS) > 0:
		somatic_inputs_tonly = [f'results/{run}/somatic/{patient}/somatic_annotated_tonly.vcf.gz' for run, patient in product([config['run']], ngs.ONLY_TMR_PATIENTS)]
	
	# metrics
	metrics = [f"results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv" for run, sample in product([config['run']], ngs.SAMPLES)]

	if ngs.GRM and ngs.TMR is False:
		input_files = germline_inputs
	elif ngs.GRM is False and ngs.TMR:
		input_files = somatic_inputs_grm_vs_tmr + somatic_inputs_tonly
	elif ngs.GRM and ngs.TMR:
		input_files = germline_inputs + somatic_inputs_grm_vs_tmr + somatic_inputs_tonly

	if config['ngs_type'] != 'WGS':
		input_files = input_files + metrics

	return input_files
