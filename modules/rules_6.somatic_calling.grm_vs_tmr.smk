rule mutect2_grm_vs_tmr:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		bam_tmr=lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor'],
		bam_grm=lambda wc: get_somatic_input(wc, ngs.wide_df)['germline'],
		capture="results/{run}/capture.intervals"
	output: 
		vcf_raw=temp('results/{run}/somatic/{patient}/raw.vcf.gz'),
		bam=temp('results/{run}/somatic/{patient}/raw.bam')
	log: 
		'results/{run}/logs/somatic/{patient}/Mutect2.log'
	params:
		sample_name=lambda wc: ngs.wide_df.loc[:, 'grm_samples'][ngs.wide_df.loc[:, 'patients']==wc.patient].values[0],
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		grm_res=config['references38']['af_only_gnomad'] if config['assembly'] == 'GRCh38' else config['references37']['af_only_gnomad'],
		pon=config['references38']['snps'] if config['assembly'] == 'GRCh38' else config['references37']['snps']
	shell:
		"""gatk Mutect2 \
				-R {params.ref} \
				-I {input.bam_tmr} \
				-O {output.vcf_raw} \
				-I {input.bam_grm} \
				-normal {params.sample_name} \
				-bamout {output.bam} \
				--germline-resource {params.grm_res} \
				--panel-of-normals {params.pon} 2>{log}"""

use rule get_pileup_summaries_tmr as get_pileup_summaries_grm with:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		bam=lambda wc: get_somatic_input(wc, ngs.wide_df)['germline']
	output: 
		getpileupsum=temp('results/{run}/somatic/{patient}/getpileupsummaries_grm.table')
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_germline.log'


rule calculate_contamination_grm_vs_tmr:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		getpileupsum_tmr='results/{run}/somatic/{patient}/getpileupsummaries_tmr.table',
		getpileupsum_grm='results/{run}/somatic/{patient}/getpileupsummaries_grm.table'
	output: 
		contamination=temp('results/{run}/somatic/{patient}/contamination.table')
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell: 
		"""gatk CalculateContamination \
				-I {input.getpileupsum_tmr} \
				-O {output.contamination} \
				-matched {input.getpileupsum_grm} 2>{log}"""

rule filter_mutect_calls_grm_vs_tmr:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		vcf_raw='results/{run}/somatic/{patient}/raw.vcf.gz',
		contamination='results/{run}/somatic/{patient}/contamination.table',
		rom='results/{run}/somatic/{patient}/read-orientation-model.tar.gz'
	output: 
		vcf_filt=temp('results/{run}/somatic/{patient}/filtered.vcf.gz')
	log: 
		'results/{run}/logs/somatic/{patient}/FilterMutectCalls.log'
	params:
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa']
	shell: 
		"""gatk FilterMutectCalls \
				-R {params.ref} \
				-V {input.vcf_raw} \
				-O {output.vcf_filt} \
				--contamination-table {input.contamination} \
				--orientation-bias-artifact-priors {input.rom} 2>{log}
				"""

rule filter_pass_exclude_normal_grm_vs_tmr:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		vcf_filt='results/{run}/somatic/{patient}/filtered.vcf.gz'
	output: 
		vcf_final=temp('results/{run}/somatic/{patient}/final.vcf.gz')
	params: 
		bcftools=config['tools']['bcftools']
	threads: 
		workflow.cores/8
	shell: 
		"""tumor_sample=$(zcat {input.vcf_filt} | \
					grep -oP '##tumor_sample=\K.+') &&\
					{params.bcftools} view -i "%FILTER='PASS'" \
					-s $tumor_sample -t {threads} \
					-Oz -o {output.vcf_final} {input.vcf_filt}"""

rule generate_somatic_plots:
	input: 
		'results/{run}/somatic/{patient}/filtered.pass.vcf.gz'
	output: 
		'results/{run}/somatic/{patient}/plots/plot.pdf'
	log: 
		'results/{run}/logs/somatic/{patient}/plotting.log'
	script: 
		"scripts/somatic_plots.R"
