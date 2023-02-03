
rule mutect2_tumor_only:
	wildcard_constraints:
		patient="|".join(ONLY_TMR_PATIENTS)
	input: 
		bam_tmr=lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor'],
		capture="results/{run}/capture.intervals"
	output: 
		vcf_raw=temp('results/{run}/somatic/{patient}/raw.vcf.gz'),
		bam=temp('results/{run}/somatic/{patient}/raw.bam')
	log: 
		'results/{run}/logs/somatic/{patient}/Mutect2.log'
	params:
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		grm_res=config['references38']['af_only_gnomad'] if config['assembly'] == 'GRCh38' else config['references37']['af_only_gnomad'],
		pon=config['references38']['snps'] if config['assembly'] == 'GRCh38' else config['references37']['snps']
	shell: 
		"""gatk Mutect2 \
				-R {params.ref} \
				-I {input.bam_tmr} \
				-O {output.vcf_raw} \
				-bamout {output.bam} \
				--germline-resource {params.grm_res} \
				--panel-of-normals {params.pon} 2>{log}"""

rule calculate_contamination_tmr_only:
	wildcard_constraints:
		patient="|".join(ONLY_TMR_PATIENTS)
	input: 
		getpileupsum_tmr='results/{run}/somatic/{patient}/getpileupsummaries_tmr.table'
	output: 
		contamination=temp('results/{run}/somatic/{patient}/contamination.table')
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell: 
		"""gatk CalculateContamination \
				-I {input.getpileupsum_tmr} \
				-O {output.contamination} 2>{log}"""


use rule filter_mutect_calls_grm_vs_tmr as filter_mutect_calls_tmr_only with:
	wildcard_constraints:
		patient="|".join(ONLY_TMR_PATIENTS)
	output:
		vcf_filt=temp('results/{run}/somatic/{patient}/final.vcf.gz')
