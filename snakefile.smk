from itertools import product
from modules.scripts.setup_run import setup_run


configfile: 'configure.yml'

mapping, long_df = setup_run()

SAMPLES = long_df['sample'].unique().tolist()
LANES = long_df['lane'].unique().tolist()

GRM_SAMPLES = mapping['sample_grm'].dropna().tolist()
TMR_SAMPLES = mapping['sample_tmr'].dropna().tolist()
ALL_PATIENTS = {i[:-4] for i in TMR_SAMPLES}
GRM_VS_TMR_PATIENTS = mapping.query("~(sample_tmr.isnull() | sample_grm.isnull())")['patient'].to_list()
ONLY_TMR_PATIENTS =  mapping.query("sample_grm.isnull()")['patient'].to_list()

GRM=config['grm_dir'] != ''
TMR=config['tmr_dir'] != ''
PAIR=config['reads_type'] == 'pair'

def final_inputs():
	germline_inputs, somatic_inputs, metrics = [], [], []

	# germline
	germline_inputs = [f'results/{run}/germline/vcf/{sample}.annotated.vcf.gz' for run, sample in product([config['run']], GRM_SAMPLES)]
	if len(GRM_SAMPLES) > 1:
		germline_inputs = germline_inputs + [f'results/{run}/germline/vcf/cohort.annotated.vcf.gz' for run in [config['run']]]
	# somatic
	somatic_inputs = [f'results/{run}/somatic/{patient}/annotated.vcf.gz' for run, patient in product([config['run']], ALL_PATIENTS)]
	# metrics
	metrics = [f"results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv" for run, sample in product([config['run']], SAMPLES)]

	if GRM and TMR is False:
		return germline_inputs + metrics
	elif GRM is False and TMR:
		return somatic_inputs + metrics
	elif GRM and TMR:
		return germline_inputs + somatic_inputs + metrics

wildcard_constraints:
	sample="|".join(SAMPLES),
	patient="|".join(ALL_PATIENTS)

rule all:
	input: final_inputs()

include: config["snakemake_modules"] + "rules_1.aligning.smk"
include: config["snakemake_modules"] + "rules_2.preprocessing.smk"
include: config["snakemake_modules"] + "rules_3.collect_metrics.smk"
include: config["snakemake_modules"] + "rules_4.germline_calling.deepvariant.smk"
include: config["snakemake_modules"] + "rules_5.somatic_calling.smk"
include: config["snakemake_modules"] + "rules_6.somatic_calling.grm_vs_tmr.smk"
include: config["snakemake_modules"] + "rules_7.somatic_calling.tmr_only.smk"
# include: config["snakemake_modules"] + "rules_8.sv_calling.smk"
include: config["snakemake_modules"] + "rules_9.annotation.smk"
