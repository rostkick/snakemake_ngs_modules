def get_somatic_input(wc, data):
	sample_tumor = data.loc[:, 'tmr_samples'][data['patients']==wc.patient].values[0]
	sample_germline = data.loc[:, 'grm_samples'][data.loc[:, 'patients']==wc.patient].values[0]
	return {'tumor': f'results/{wc.run}/bam/{sample_tumor}.final.bam',
			'germline': f'results/{wc.run}/bam/{sample_germline}.final.bam'}

rule CollectF1R2Counts:
	input: lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor']
	output: 'results/{run}/somatic/{patient}/tumor-artifact-prior-table.tsv'
	log: 'results/{run}/logs/somatic/{patient}/CollectF1R2Counts.log'
	params: reference=config['references']['genome_fa']
	shell: """gatk CollectF1R2Counts \
				-R {params.reference} \
				-I {input} \
				-O {output} 2>{log}"""

rule LearnReadOrientationModel:
	input: 'results/{run}/somatic/{patient}/tumor-artifact-prior-table.tsv'
	output: 'results/{run}/somatic/{patient}/tumor-sample-artifact-prior.tsv.tar.gz'
	log: 'results/{run}/logs/somatic/{patient}/LearnReadOrientationModel.log'
	shell: """gatk LearnReadOrientationModel \
				-I {input} \
				-O {output} 2>{log}"""

rule Mutect2_tumor_vs_normal:
	input: 
		bam_tumor=lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor'],
		bam_germline=lambda wc: get_somatic_input(wc, ngs.wide_df)['germline'],
		capture="results/{run}/capture.intervals"
	output: 
		mutect2_raw='results/{run}/somatic/{patient}/mutect2.vcf.gz',
		bam='results/{run}/somatic/{patient}/mutect2.bamout.bam'
	log: 'results/{run}/logs/somatic/{patient}/Mutect2.log'
	params:
		reference_fasta=config['references']['genome_fa'],
		sample_name=lambda wc: ngs.wide_df.loc[:, 'grm_samples'][ngs.wide_df.loc[:, 'patients']==wc.patient].values[0],
		germline_res=config['references']['af_only_gnomad'],
		panel_of_normals=config['references']['snps']
	shell: """gatk Mutect2 \
				-R {params.reference_fasta} \
				-I {input.bam_tumor} \
				-O {output.mutect2_raw} \
				-I {input.bam_germline} \
				-normal {params.sample_name} \
				-bamout {output.bam} \
				-L {input.capture} \
				--germline-resource {params.germline_res} \
				--panel-of-normals {params.panel_of_normals} 2>{log}"""

rule GetPileupSummaries_tumor:
	input: lambda wc: get_somatic_input(wc, ngs.wide_df)['tumor']
	output: 'results/{run}/somatic/{patient}/tumor.pileups.table'
	log: 'results/{run}/logs/somatic/{patient}/GetPileupSummaries_tumor.log'
	params: 
		reference_fasta=config['references']['genome_fa'],
		small_exac_common=config['references']['small_exac_common']
	shell: """gatk GetPileupSummaries \
				-I {input} \
				-R {params.reference_fasta} \
				-V {params.small_exac_common} \
				-L {params.small_exac_common} \
				-O {output} 2>{log}"""

rule GetPileupSummaries_germline:
	input: lambda wc: get_somatic_input(wc, ngs.wide_df)['germline']
	output: 'results/{run}/somatic/{patient}/germline.pileups.table'
	log: 'results/{run}/logs/somatic/{patient}/GetPileupSummaries_germline.log'
	params: 
		reference_fasta=config['references']['genome_fa'],
		small_exac_common=config['references']['small_exac_common']
	shell: """gatk GetPileupSummaries \
				-I {input} \
				-R {params.reference_fasta} \
				-V {params.small_exac_common} \
				-L {params.small_exac_common} \
				-O {output} 2>{log}"""

rule CalculateContamination:
	input: 
		tumor='results/{run}/somatic/{patient}/tumor.pileups.table',
		germline='results/{run}/somatic/{patient}/germline.pileups.table'
	output: 'results/{run}/somatic/{patient}/contamination.table'
	log: 'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell: """gatk CalculateContamination \
				-I {input.tumor} \
				-O {output} \
				-matched {input.germline} 2>{log}"""

rule FilterMutectCalls:
	input: 
		mutect2_raw='results/{run}/somatic/{patient}/mutect2.vcf.gz',
		contamination_table='results/{run}/somatic/{patient}/contamination.table',
		orientation_artifacts='results/{run}/somatic/{patient}/tumor-sample-artifact-prior.tsv.tar.gz'
	output: 'results/{run}/somatic/{patient}/mutect2.filtered.vcf.gz'
	log: 'results/{run}/logs/somatic/{patient}/FilterMutectCalls.log'
	params: 
		reference_fasta=config['references']['genome_fa']
	shell: """gatk FilterMutectCalls \
				-R {params.reference_fasta} \
				-V {input.mutect2_raw} \
				-O {output} \
				--contamination-table {input.contamination_table} \
				--orientation-bias-artifact-priors {input.orientation_artifacts} 2>{log}
				"""

rule FilterPASS_exclude_normal:
	input: 'results/{run}/somatic/{patient}/mutect2.filtered.vcf.gz'
	output: 'results/{run}/somatic/{patient}/mutect2.final.vcf.gz'
	params: bcftools=config['tools']['bcftools']
	threads: workflow.cores/8
	shell: """tumor_sample=$(zcat {input} | \
					grep -oP '##tumor_sample=\K.+') &&\
					{params.bcftools} view -i "%FILTER='PASS'" \
					-s $tumor_sample -t {threads} \
					-Oz -o {output} {input}"""

rule generate_somatic_plots:
	input: 'results/{run}/somatic/{patient}/mutect2.filtered.pass.vcf.gz'
	output: 'results/{run}/somatic/{patient}/plots/plot.pdf'
	log: 'results/{run}/logs/somatic/{patient}/plotting.log'
	script: "scripts/somatic_plots.R"
