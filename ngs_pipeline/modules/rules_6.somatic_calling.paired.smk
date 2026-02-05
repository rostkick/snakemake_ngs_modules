rule r6_1_mutect2_paired:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		bam_tmr = lambda wc: get_somatic_input(wc, ngs.data)['tumor'],
		bam_grm = lambda wc: get_somatic_input(wc, ngs.data)['germline'],
		capture = rules.r2_0_bed_to_intervals.output.intervals if config['panel_capture']['target'][-4:] == 'bed' else config['panel_capture']['target']
	output: 
		vcf = temp('results/{run}/somatic/{patient}/raw.vcf.gz'),
		bam = temp('results/{run}/somatic/{patient}/raw.bam')
	log: 
		'results/{run}/logs/somatic/{patient}/Mutect2.log'
	params:
		gatk = config['tools']['gatk'],
		sample_name = lambda wc: wc.patient + '_grm',
		ref = config['references']['genome_fa'],
		grm_res = config['references']['af_only_gnomad'],
		pon = config['references']['snps']
	shell:
		"""{params.gatk} Mutect2 \
				-R {params.ref} \
				-I {input.bam_tmr} \
				-O {output.vcf} \
				-I {input.bam_grm} \
				-normal {params.sample_name} \
				-bamout {output.bam} \
				--germline-resource {params.grm_res} \
				--panel-of-normals {params.pon} 2>{log}"""

use rule r5_3_get_pileup_summaries_tmr as r6_2_get_pileup_summaries_grm with:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		bam = lambda wc: get_somatic_input(wc, ngs.data)['germline']
	output: 
		getpileupsum = temp('results/{run}/somatic/{patient}/getpileupsummaries_grm.table')
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_germline.log'

rule r6_3_get_pileup_summaries_grm:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		getpileupsum_tmr = rules.r5_3_get_pileup_summaries_tmr.output.getpileupsum,
		getpileupsum_grm = rules.r6_2_get_pileup_summaries_grm.output.getpileupsum
	output: 
		contamination = 'results/{run}/somatic/{patient}/contamination.table'
	params:
		gatk = config['tools']['gatk']
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell:"""
		{params.gatk} CalculateContamination \
				-I {input.getpileupsum_tmr} \
				-O {output.contamination} \
				-matched {input.getpileupsum_grm} 2>{log}"""

rule r6_4_filter_mutect_calls_paired:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		vcf = rules.r6_1_mutect2_paired.output.vcf,
		contamination = rules.r6_3_get_pileup_summaries_grm.output.contamination,
		rom = rules.r5_2_learn_read_orientation_model.output.rom
	output: 
		vcf = 'results/{run}/somatic/{patient}/filtered.vcf.gz'
	log: 
		'results/{run}/logs/somatic/{patient}/FilterMutectCalls.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa']
	shell:"""
			{params.gatk} FilterMutectCalls \
				-R {params.ref} \
				-V {input.vcf} \
				-O {output.vcf} \
				--contamination-table {input.contamination} \
				--orientation-bias-artifact-priors {input.rom} 2>{log}
				"""

rule r6_5_filter_pass_exclude_normal_paired:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		vcf = rules.r6_4_filter_mutect_calls_paired.output.vcf
	output: 
		vcf = 'results/{run}/somatic/{patient}/final.vcf.gz'
	params: 
		bcftools = config['tools']['bcftools']
	# threads: 
	# 	workflow.cores/max(len(ngs.GRM_VS_TMR_PATIENTS), 1)
	shell:"""
			tumor_sample=$(zcat {input.vcf} | \
					grep -oP '##tumor_sample=\K.+') &&\
					{params.bcftools} view -i "%FILTER='PASS'" \
					-s $tumor_sample \
					-Oz -o {output.vcf} {input.vcf}"""
