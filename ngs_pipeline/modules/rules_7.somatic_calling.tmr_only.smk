# ===========================================
# ===== Tumor-only mode variant calling =====
# ===========================================
rule r7_1_mutect2_tmr_only:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		bam = lambda wc: get_somatic_input(wc, ngs.data)['tumor'],
		capture = rules.r2_1_bed_to_intervals.output.intervals if config['panel_capture']['target'][-4:] == 'bed' else config['panel_capture']['target']
	output: 
		vcf = temp('results/{run}/somatic/{patient}/raw_tonly.vcf.gz'),
		bam = temp('results/{run}/somatic/{patient}/raw_tonly.bam')
	log: 
		'results/{run}/logs/somatic/{patient}/Mutect2.log'
	priority: 35
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		grm_res = config['references']['af_only_gnomad'],
		pon = config['references']['pon'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.6)}m"
	resources:
		mem_mb={'panel': 8000, 'WES': 16000, 'WGS': 20000}.get(config['ngs_type'], 16000),
		runtime_min={'panel': 1440, 'WES': 5760, 'WGS': 11520}.get(config['ngs_type'], 5760)
	shell: """
			{params.gatk} --java-options "{params.java_opts}" Mutect2 \
				-R {params.ref} \
				-I {input.bam} \
				-O {output.vcf} \
				-bamout {output.bam} \
				--germline-resource {params.grm_res} \
				--af-of-alleles-not-in-resource 0.0001 \
				--panel-of-normals {params.pon} 2>{log}"""

rule r7_2_calculate_contamination_tmr_only:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		getpileupsum_tmr = rules.r5_3_get_pileup_summaries_tmr.output.getpileupsum
	output: 
		contamination = temp('results/{run}/somatic/{patient}/contamination_tonly.table')
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
				-O {output.contamination} 2>{log}"""


use rule r6_4_filter_mutect_calls_paired as r7_3_filter_mutect_calls_tmr_only with:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input:
		vcf = rules.r7_1_mutect2_tmr_only.output.vcf,
		contamination = rules.r7_2_calculate_contamination_tmr_only.output.contamination,
		rom = rules.r5_2_learn_read_orientation_model.output.rom
	output:
		vcf = temp('results/{run}/somatic/{patient}/final_tonly.vcf.gz')
	resources:
		mem_mb={'panel': 2000, 'WES': 4000, 'WGS': 6000}.get(config['ngs_type'], 4000),
		runtime_min={'panel': 60, 'WES': 240, 'WGS': 480}.get(config['ngs_type'], 240)
	priority: 35
