rule CollectF1R2Counts:
	input: lambda wc: get_somatic_input(wc, data)['tumor']
	output: 'results/{run}/somatic/{patient}/tumor-artifact-prior-table.tsv'
	log: 'results/{run}/logs/somatic/{patient}/CollectF1R2Counts.log'
	params: reference=config['GRCh38']['GATK_b38']['reference_fasta']
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

rule Mutect2:
	input: 
		bam_tumor=lambda wc: get_somatic_input(wc, data)['tumor'],
		bam_germline=lambda wc: get_somatic_input(wc, data)['germline']
	output: 
		mutect2_raw='results/{run}/somatic/{patient}/mutect2.vcf.gz',
		bam='results/{run}/somatic/{patient}/mutect2.bamout.bam'
	log: 'results/{run}/logs/somatic/{patient}/Mutect2.log'
	params:
		reference_fasta=config['GRCh38']['GATK_b38']['reference_fasta'],
		sample_name=lambda wc: data.loc[:, 'germline_samples'][data.loc[:, 'patient']==wc.patient].values[0],
		germline_res=config['GRCh38']['somatic_hg38']['af_only_gnomad'],
		panel_of_normals=config['GRCh38']['GATK_b38']['snps']
	shell: """gatk Mutect2 \
				-R {params.reference_fasta} \
				-I {input.bam_tumor} \
				-O {output.mutect2_raw} \
				-I {input.bam_germline} \
				-normal {params.sample_name} \
				-bamout {output.bam} \
				--germline-resource {params.germline_res} \
				--panel-of-normals {params.panel_of_normals} 2>{log}"""

rule GetPileupSummaries_tumor:
	input: lambda wc: get_somatic_input(wc, data)['tumor']
	output: 'results/{run}/somatic/{patient}/tumor.pileups.table'
	log: 'results/{run}/logs/somatic/{patient}/GetPileupSummaries_tumor.log'
	params: 
		reference_fasta=config['GRCh38']['GATK_b38']['reference_fasta'],
		small_exac_common=config['GRCh38']['somatic_hg38']['small_exac_common']
	shell: """gatk GetPileupSummaries \
				-I {input} \
				-R {params.reference_fasta} \
				-V {params.small_exac_common} \
				-L {params.small_exac_common} \
				-O {output} 2>{log}"""

rule GetPileupSummaries_germline:
	input: lambda wc: get_somatic_input(wc, data)['germline']
	output: 'results/{run}/somatic/{patient}/germline.pileups.table'
	log: 'results/{run}/logs/somatic/{patient}/GetPileupSummaries_germline.log'
	params: 
		reference_fasta=config['GRCh38']['GATK_b38']['reference_fasta'],
		small_exac_common=config['GRCh38']['somatic_hg38']['small_exac_common']
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
		reference_fasta=config['GRCh38']['GATK_b38']['reference_fasta']
	shell: """gatk FilterMutectCalls \
				-R {params.reference_fasta} \
				-V {input.mutect2_raw} \
				-O {output} \
				--contamination-table {input.contamination_table} \
				--orientation-bias-artifact-priors {input.orientation_artifacts} 2>{log}
				"""

rule FilterPASS_exclude_normal:
	input: 'results/{run}/somatic/{patient}/mutect2.filtered.vcf.gz'
	output: 'results/{run}/somatic/{patient}/mutect2.filtered.pass.vcf.gz'
	threads: workflow.cores/8
	shell: """tumor_sample=$(zcat {input} | \
					grep -oP '##tumor_sample=\K.+') &&\
				bcftools view -i "%FILTER='PASS'" \
					-s $tumor_sample -t {threads} \
					-Oz -o {output} {input}"""

rule generate_somatic_plots:
	input: 'results/{run}/somatic/{patient}/mutect2.filtered.pass.vcf.gz'
	output: 'results/{run}/somatic/{patient}/plots/plot.pdf'
	log: 'results/{run}/logs/somatic/{patient}/plotting.log'
	script: "scripts/somatic_plots.R"
