rule manta_germline_joint_calling:
	input: expand('results/{{run}}/bam/{sample}.final.bam', sample=GERMLINE_SAMPLES)
	output: 'results/{run}/germline/sv/results/variants/diploidSV.vcf.gz'
	params: 
		manta=config['tools']['manta']['configManta'],
		input_bams=lambda wc, input: [f'--bam {i}' for i in input],
		reference_fasta=config['GRCh38']['GATK_b38']['reference_fasta'],
		rundir='results/{run}/germline/sv'
	log: 'results/{run}/logs/germline/manta.log'
	conda: "env/manta.yml"
	threads: workflow.cores/8
	shell: """rm -rf {params.rundir} && \
			{params.manta} \
				{params.input_bams} \
				--referenceFasta {params.reference_fasta} \
				--runDir {params.rundir} 2>{log} && \
			python {params.rundir}/runWorkflow.py -j {threads} 2>>{log}"""

rule convertInversion:
	input: 'results/{run}/germline/sv/results/variants/diploidSV.vcf.gz'
	output: 'results/{run}/germline/sv/results/variants/diploidSV.inv_converted.vcf.gz'
	params:
		convertInversion=config['tools']['manta']['convertInversion'],
		samtools=config['tools']['samtools'],
		reference_fasta=config['GRCh38']['GATK_b38']['reference_fasta']
	conda: "env/manta.yml"
	shell: "python {params.convertInversion} {params.samtools} {params.reference_fasta} {input} > {output}"

# rule sv_somatic_single_calling:
#     input: lambda wc: get_somatic_input(wc, data)['tumor']
#     output: 'results/{run}/somatic/{patient}/sv/sv.vcf.gz'
#     params: 
#         manta=config['tools']['manta'],
#         input_bams=lambda wc, input: [f'--bam {i}' for i in input],
#         reference_fasta=config['GRCh38']['GATK_b38']['reference_fasta'],
#         rundir='results/{run}/somatic/{patient}/sv'
#     conda: "modules/envs/manta.yml"
#     shell: """{params.manta} \
#             {params.input_bams} \
#             --referenceFasta {params.reference_fasta} \
#             --runDir {params.rundir}"""
