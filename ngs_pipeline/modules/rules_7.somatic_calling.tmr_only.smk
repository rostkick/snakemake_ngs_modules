# ===========================================
# ===== Tumor-only mode variant calling =====
# ===========================================
rule r7_1_mutect2_tmr_only:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		bam = lambda wc: get_somatic_input(wc, ngs.data)['tumor'],
		capture = rules.r2_8_bed_to_intervals.output.intervals if config['panel_capture']['target'][-4:] == 'bed' else config['panel_capture']['target']
	output: 
		vcf = 'results/{run}/somatic/{patient}/raw_tonly.vcf.gz',
		bam = 'results/{run}/somatic/{patient}/raw_tonly.bam'
	log: 
		'results/{run}/logs/somatic/{patient}/Mutect2.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		grm_res = config['references38']['af_only_gnomad'] if config['assembly'] == 'GRCh38' else config['references37']['af_only_gnomad'],
		pon = config['references38']['snps'] if config['assembly'] == 'GRCh38' else config['references37']['snps']
	shell: """
			{params.gatk} Mutect2 \
				-R {params.ref} \
				-I {input.bam} \
				-O {output.vcf} \
				-bamout {output.bam} \
				--germline-resource {params.grm_res} \
				--panel-of-normals {params.pon} 2>{log}"""

rule r7_2_calculate_contamination_tmr_only:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		getpileupsum_tmr = rules.r5_3_get_pileup_summaries_tmr.output.getpileupsum
	output: 
		contamination = 'results/{run}/somatic/{patient}/contamination_tonly.table'
	params:
		gatk = config['tools']['gatk']
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell:"""
			{params.gatk} CalculateContamination \
				-I {input.getpileupsum_tmr} \
				-O {output.contamination} 2>{log}"""


use rule r6_4_filter_mutect_calls_paired as r7_3_filter_mutect_calls_tmr_only with:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input:
		vcf = rules.r7_1_mutect2_tmr_only.output.vcf,
		contamination = rules.r7_2_calculate_contamination_tmr_only.output.contamination,
		rom = rules.r5_2_learn_read_orientation_model.output.rom
	output:
		vcf = 'results/{run}/somatic/{patient}/final_tonly.vcf.gz'
