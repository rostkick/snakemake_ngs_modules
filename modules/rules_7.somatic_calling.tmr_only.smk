
rule r7_mutect2_tumor_only:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		bam = lambda wc: get_somatic_input(wc, mapping)['tumor'],
		capture = rules.r2_bed_to_intervals.output.intervals
	output: 
		vcf = 'results/{run}/somatic/{patient}/raw.vcf.gz',
		bam = 'results/{run}/somatic/{patient}/raw.bam'
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

rule r7_calculate_contamination_tmr_only:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		getpileupsum_tmr = rules.r5_get_pileup_summaries_tmr.output.getpileupsum
	output: 
		contamination = 'results/{run}/somatic/{patient}/contamination.table'
	params:
		gatk = config['tools']['gatk']
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell:"""
			{params.gatk} CalculateContamination \
				-I {input.getpileupsum_tmr} \
				-O {output.contamination} 2>{log}"""


use rule r6_filter_mutect_calls_grm_vs_tmr as r7_filter_mutect_calls_tmr_only with:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		vcf = rules.r7_mutect2_tumor_only.output.vcf,
		contamination = rules.r7_calculate_contamination_tmr_only.output.contamination,
		rom = rules.r5_learn_read_orientation_model.output.rom
	output:
		vcf = 'results/{run}/somatic/{patient}/final.vcf.gz'
