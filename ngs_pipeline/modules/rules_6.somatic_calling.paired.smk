rule r6_1_mutect2_paired:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		bam_tmr = lambda wc: get_somatic_input(wc, ngs.data)['tumor'],
		bam_grm = lambda wc: get_somatic_input(wc, ngs.data)['germline'],
		capture = rules.r2_1_bed_to_intervals.output.intervals if config['panel_capture']['target'][-4:] == 'bed' else config['panel_capture']['target']
	output: 
		vcf = temp('results/{run}/somatic/{patient}/raw.vcf.gz'),
		bam = temp('results/{run}/somatic/{patient}/raw.bam')
	log: 
		'results/{run}/logs/somatic/{patient}/Mutect2.log'
	priority: 35
	params:
		gatk = config['tools']['gatk'],
		sample_name = lambda wc: wc.patient + '_grm',
		ref = config['references']['genome_fa'],
		grm_res = config['references']['af_only_gnomad'],
		pon = config['references']['pon'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.6)}m"
	resources:
		mem_mb={'panel': 8000, 'WES': 16000, 'WGS': 20000}.get(config['ngs_type'], 16000),
		runtime_min={'panel': 1440, 'WES': 5760, 'WGS': 11520}.get(config['ngs_type'], 5760)
	shell:
		"""{params.gatk} --java-options "{params.java_opts}" Mutect2 \
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
	resources:
		mem_mb={'panel': 4000, 'WES': 8000, 'WGS': 12000}.get(config['ngs_type'], 8000),
		runtime_min={'panel': 240, 'WES': 1440, 'WGS': 2880}.get(config['ngs_type'], 1440)
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_germline.log'
	priority: 35

rule r6_3_get_pileup_summaries_grm:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		getpileupsum_tmr = rules.r5_3_get_pileup_summaries_tmr.output.getpileupsum,
		getpileupsum_grm = rules.r6_2_get_pileup_summaries_grm.output.getpileupsum
	output: 
		contamination = temp('results/{run}/somatic/{patient}/contamination.table')
	params:
		gatk = config['tools']['gatk'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	resources:
		mem_mb={'panel': 2000, 'WES': 4000, 'WGS': 4000}.get(config['ngs_type'], 4000),
		runtime_min=120
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	priority: 35
	shell:"""
		{params.gatk} --java-options "{params.java_opts}" CalculateContamination \
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
		vcf = temp('results/{run}/somatic/{patient}/filtered.vcf.gz')
	log: 
		'results/{run}/logs/somatic/{patient}/FilterMutectCalls.log'
	priority: 35
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	resources:
		mem_mb={'panel': 2000, 'WES': 4000, 'WGS': 6000}.get(config['ngs_type'], 4000),
		runtime_min={'panel': 60, 'WES': 240, 'WGS': 480}.get(config['ngs_type'], 240)
	shell:"""
			{params.gatk} --java-options "{params.java_opts}" FilterMutectCalls \
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
		vcf = temp('results/{run}/somatic/{patient}/final.vcf.gz')
	params: 
		bcftools = config['tools']['bcftools']
	resources:
		mem_mb=2000,
		runtime_min=60
	priority: 35
	# threads: 
	# 	workflow.cores/max(len(ngs.GRM_VS_TMR_PATIENTS), 1)
	shell:"""
			set -euo pipefail
			tumor_sample=$(zcat {input.vcf} | \
					grep -oP '##tumor_sample=\K.+') &&\
					{params.bcftools} view -i "%FILTER='PASS'" \
					-s $tumor_sample \
					-Oz -o {output.vcf} {input.vcf}"""
