def get_somatic_input(wc, data):
	sample_tumor = data.loc[:, 'tmr_samples'][data['patients']==wc.patient].values[0]
	sample_germline = data.loc[:, 'grm_samples'][data.loc[:, 'patients']==wc.patient].values[0]
	return {'tumor': f'results/{wc.run}/bam/{sample_tumor}.final.bam',
			'germline': f'results/{wc.run}/bam/{sample_germline}.final.bam'}


rule Mutect2_tumor_vs_normal:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		bam_tmr=lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor'],
		bam_grm=lambda wc: get_somatic_input(wc, ngs.wide_df)['germline'],
		capture="results/{run}/capture.intervals"
	output: 
		vcf_raw='results/{run}/somatic/{patient}/raw.vcf.gz',
		bam='results/{run}/somatic/{patient}/raw.bam'
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
				-L {input.capture} \
				--germline-resource {params.grm_res} \
				--panel-of-normals {params.pon} 2>{log}"""

rule Mutect2_tumor_only:
	wildcard_constraints:
		patient="|".join(ONLY_TMR_PATIENTS)
	input: 
		bam_tmr=lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor'],
		capture="results/{run}/capture.intervals"
	output: 
		vcf_raw='results/{run}/somatic/{patient}/raw.vcf.gz',
		bam='results/{run}/somatic/{patient}/raw.bam'
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
				-L {input.capture} \
				--germline-resource {params.grm_res} \
				--panel-of-normals {params.pon} 2>{log}"""

rule CollectF1R2Counts:
	wildcard_constraints:
		patient="|".join(ALL_PATIENTS)
	input: 
		bam_tmr=lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor']
	output: 
		f1r2='results/{run}/somatic/{patient}/f1r2.tsv'
	log: 
		'results/{run}/logs/somatic/{patient}/CollectF1R2Counts.log'
	params:
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
	shell:
		"""gatk CollectF1R2Counts \
				-R {params.ref} \
				-I {input.bam_tmr} \
				-O {output.f1r2} 2>{log}"""


rule LearnReadOrientationModel:
	wildcard_constraints:
		patient="|".join(ALL_PATIENTS)
	input: 
		f1r2='results/{run}/somatic/{patient}/f1r2.tsv'
	output: 
		rom='results/{run}/somatic/{patient}/read-orientation-model.tar.gz'
	log: 
		'results/{run}/logs/somatic/{patient}/LearnReadOrientationModel.log'
	shell: 
		"""gatk LearnReadOrientationModel \
				-I {input.f1r2} \
				-O {output.rom} 2>{log}"""

rule GetPileupSummaries_tmr:
	wildcard_constraints:
		patient="|".join(ALL_PATIENTS)
	input: 
		bam_tmr=lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor']
	output: 
		getpileupsum_tmr='results/{run}/somatic/{patient}/getpileupsummaries_tmr.table'
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_tumor.log'
	params: 
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		exac=config['references38']['small_exac_common'] if config['assembly'] == 'GRCh38' else config['references37']['small_exac_common']
	shell: 
		"""gatk GetPileupSummaries \
				-I {input.bam_tmr} \
				-R {params.ref} \
				-V {params.exac} \
				-L {params.exac} \
				-O {output.getpileupsum} 2>{log}"""

rule GetPileupSummaries_grm:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		bam_grm=lambda wc: get_somatic_input(wc, ngs.wide_df)['germline']
	output: 
		getpileupsum_grm='results/{run}/somatic/{patient}/getpileupsummaries_grm.table'
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_germline.log'
	params: 
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		exac=config['references38']['small_exac_common'] if config['assembly'] == 'GRCh38' else config['references37']['small_exac_common']
	shell: 
		"""gatk GetPileupSummaries \
				-I {input.bam_grm} \
				-R {params.ref} \
				-V {params.exac} \
				-L {params.exac} \
				-O {output.getpileupsum} 2>{log}"""

rule CalculateContamination_grm_vs_tmr:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		getpileupsum_tmr='results/{run}/somatic/{patient}/getpileupsummaries_tmr.table',
		getpileupsum_grm='results/{run}/somatic/{patient}/getpileupsummaries_grm.table'
	output: 
		contamination='results/{run}/somatic/{patient}/contamination.table'
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell: 
		"""gatk CalculateContamination \
				-I {input.getpileupsum_tmr} \
				-O {output.contamination} \
				-matched {input.getpileupsum_grm} 2>{log}"""

rule CalculateContamination_tmr_only:
	wildcard_constraints:
		patient="|".join(ONLY_TMR_PATIENTS)
	input: 
		getpileupsum_tmr='results/{run}/somatic/{patient}/getpileupsummaries_tmr.table'
	output: 
		contamination='results/{run}/somatic/{patient}/contamination.table'
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell: 
		"""gatk CalculateContamination \
				-I {input.getpileupsum_tmr} \
				-O {output.contamination} 2>{log}"""

rule FilterMutectCalls:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		vcf_raw='results/{run}/somatic/{patient}/raw.vcf.gz',
		contamination='results/{run}/somatic/{patient}/contamination.table',
		rom='results/{run}/somatic/{patient}/read-orientation-model.tar.gz'
	output: 
		vcf_filt='results/{run}/somatic/{patient}/filtered.vcf.gz'
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

use rule FilterMutectCalls as FilterMutectCalls_tmr_only with:
	wildcard_constraints:
		patient="|".join(ONLY_TMR_PATIENTS)
	output:
		vcf_filt='results/{run}/somatic/{patient}/final.vcf.gz'

rule FilterPASS_exclude_normal:
	wildcard_constraints:
		patient="|".join(GRM_VS_TMR_PATIENTS)
	input: 
		vcf_filt='results/{run}/somatic/{patient}/filtered.vcf.gz'
	output: 
		vcf_final='results/{run}/somatic/{patient}/final.vcf.gz'
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
