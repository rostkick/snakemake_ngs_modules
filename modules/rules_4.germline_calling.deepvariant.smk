rule r4_deepvariant:
	input: 
		bam = rules.r2_apply_bqsr.output.bam
	output: 
		vcf = "results/{run}/germline/vcf/{sample}.vcf.gz",
		gvcf = "results/{run}/germline/vcf/{sample}.gvcf.gz"
	log: 'results/{run}/logs/germline_calling/{sample}.deepvariant.log'
	params:
		singularity = config['tools']['singularity'],
		deepvariant = config['tools']['deepvariant'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		panel_capture = config['panel_capture']['target']
	threads: workflow.cores/2
	shell: """
			{params.singularity} run \
				-B /mnt:/mnt \
				{params.deepvariant} /opt/deepvariant/bin/run_deepvariant \
				--model_type=WES \
				--output_vcf={output.vcf} \
				--output_gvcf={output.gvcf} \
				--reads={input.bam} \
				--ref={params.ref} \
				--regions {params.panel_capture} \
				--num_shards={threads} &>{log} || true"""

rule r4_make_gvcf_list:
	input: 
		gvcfs = expand("results/{{run}}/germline/vcf/{sample}.gvcf.gz", sample=GRM_SAMPLES)
	output: 
		gvcfs_list = "results/{run}/germline/vcf/gvcfs.list"
	shell: "echo {input.gvcfs} | sed 's/\\s/\\n/g' > {output.gvcfs_list}"

rule r4_merge_glnexus:
	input: 
		gvcfs_list = rules.r4_make_gvcf_list.output.gvcfs_list
	output: 
		vcf = "results/{run}/germline/vcf/cohort.vcf.gz"
	log: 'results/{run}/logs/germline_calling/glnexus.log'
	params:
		glnexus_cli = config['tools']['glnexus_cli'],
		run_dir = "results/{run}/germline/vcf/GLnexus.DB",
		bcftools = config['tools']['bcftools']
	threads: workflow.cores/2
	shell: '''rm -rf {params.run_dir} &&\
				{params.glnexus_cli} \
				--config gatk_unfiltered \
				--list --trim-uncalled-alleles\
				--dir {params.run_dir} \
				-t {threads} \
				{input} 2>{log} | \
				{params.bcftools} view --min-ac 1 -i "%FILTER=='.'" -Oz -o {output} - && \
				rm -rf {params.run_dir}'''
