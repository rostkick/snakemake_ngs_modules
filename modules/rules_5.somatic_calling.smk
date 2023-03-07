from collections import defaultdict


def get_somatic_input(wc, data):
	inp = defaultdict(str)
	sample_tumor = data.loc[:, 'sample_tmr'][data['patient']==wc.patient].values[0]
	inp['tumor'] = f'results/{wc.run}/bam/{sample_tumor}.final.bam'
	if ngs.GRM:
		sample_germline = data.loc[:, 'sample_grm'][data.loc[:, 'patient']==wc.patient].values[0]
		inp['germline'] = f'results/{wc.run}/bam/{sample_germline}.final.bam'
	return inp

# here specifies the priority of germline-tumor runs vs tmr-only when both types of samples exist 
ruleorder: r6_mutect2_grm_vs_tmr > r7_mutect2_tumor_only
ruleorder: r6_calculate_contamination_grm_vs_tmr > r7_calculate_contamination_tmr_only
ruleorder: r6_filter_pass_exclude_normal_grm_vs_tmr > r7_filter_mutect_calls_tmr_only


rule r5_collect_fir2_counts:
	input: 
		bam = lambda wc: get_somatic_input(wc, mapping)['tumor']
	output: 
		f1r2 = 'results/{run}/somatic/{patient}/f1r2.tsv'
	log: 
		'results/{run}/logs/somatic/{patient}/CollectF1R2Counts.log'
	params:
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
	shell:
		"""gatk CollectF1R2Counts \
				-R {params.ref} \
				-I {input.bam} \
				-O {output.f1r2} 2>{log}"""


rule r5_learn_read_orientation_model:
	input: 
		f1r2 = rules.r5_collect_fir2_counts.output.f1r2
	output: 
		rom = temp('results/{run}/somatic/{patient}/read-orientation-model.tar.gz')
	log: 
		'results/{run}/logs/somatic/{patient}/LearnReadOrientationModel.log'
	shell: 
		"""gatk LearnReadOrientationModel \
				-I {input.f1r2} \
				-O {output.rom} 2>{log}"""

rule r5_get_pileup_summaries_tmr:
	input: 
		bam = lambda wc: get_somatic_input(wc, mapping)['tumor']
	output: 
		getpileupsum = 'results/{run}/somatic/{patient}/getpileupsummaries_tmr.table'
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_tumor.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		exac = config['references38']['small_exac_common'] if config['assembly'] == 'GRCh38' else config['references37']['small_exac_common']
	shell: """
			{params.gatk} GetPileupSummaries \
				-I {input.bam} \
				-R {params.ref} \
				-V {params.exac} \
				-L {params.exac} \
				-O {output.getpileupsum} 2>{log}"""