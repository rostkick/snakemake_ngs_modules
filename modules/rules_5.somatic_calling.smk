def get_somatic_input(wc, data):
	sample_tumor = data.loc[:, 'tmr_samples'][data['patients']==wc.patient].values[0]
	sample_germline = data.loc[:, 'grm_samples'][data.loc[:, 'patients']==wc.patient].values[0]
	return {'tumor': f'results/{wc.run}/bam/{sample_tumor}.final.bam',
			'germline': f'results/{wc.run}/bam/{sample_germline}.final.bam'}

# here specifies the priority of germline-tumor runs vs tmr-only when both types of samples exist 
ruleorder: mutect2_grm_vs_tmr > mutect2_tumor_only
ruleorder: calculate_contamination_grm_vs_tmr > calculate_contamination_tmr_only
ruleorder: filter_pass_exclude_normal_grm_vs_tmr > filter_mutect_calls_tmr_only


rule collect_fir2_counts:
	wildcard_constraints:
		patient="|".join(ALL_PATIENTS)
	input: 
		bam_tmr=lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor']
	output: 
		f1r2=temp('results/{run}/somatic/{patient}/f1r2.tsv')
	log: 
		'results/{run}/logs/somatic/{patient}/CollectF1R2Counts.log'
	params:
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
	shell:
		"""gatk CollectF1R2Counts \
				-R {params.ref} \
				-I {input.bam_tmr} \
				-O {output.f1r2} 2>{log}"""


rule learn_read_orientation_model:
	wildcard_constraints:
		patient="|".join(ALL_PATIENTS)
	input: 
		f1r2='results/{run}/somatic/{patient}/f1r2.tsv'
	output: 
		rom=temp('results/{run}/somatic/{patient}/read-orientation-model.tar.gz')
	log: 
		'results/{run}/logs/somatic/{patient}/LearnReadOrientationModel.log'
	shell: 
		"""gatk LearnReadOrientationModel \
				-I {input.f1r2} \
				-O {output.rom} 2>{log}"""

rule get_pileup_summaries_tmr:
	wildcard_constraints:
		patient="|".join(ALL_PATIENTS)
	input: 
		bam=lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor']
	output: 
		getpileupsum=temp('results/{run}/somatic/{patient}/getpileupsummaries_tmr.table')
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_tumor.log'
	params: 
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		exac=config['references38']['small_exac_common'] if config['assembly'] == 'GRCh38' else config['references37']['small_exac_common']
	shell: 
		"""gatk GetPileupSummaries \
				-I {input.bam} \
				-R {params.ref} \
				-V {params.exac} \
				-L {params.exac} \
				-O {output.getpileupsum} 2>{log}"""