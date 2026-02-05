rule r1_1_read_alignment:
	wildcard_constraints:
		sample='[\w\d]+',
		lane='[lL0-9M]+'
	input: 
		fr = lambda wc: ngs.data.loc[(ngs.data["sample"]==wc.sample) & (ngs.data["lane"]==wc.lane) & (ngs.data["reads_orientation"]=='R1'), 'fastq'].tolist(),
		rr = lambda wc: ngs.data.loc[(ngs.data["sample"]==wc.sample) & (ngs.data["lane"]==wc.lane) & (ngs.data["reads_orientation"]=='R2'), 'fastq'].tolist()
	output: 
		sam = pipe('results/{run}/bam/{sample}.{lane}.sam')
	log:
		'results/{run}/logs/aligning/{sample}.{lane}.bwa_mem2.log'
	params: 
		bwa_mem2 = config['tools']['bwa_mem2'],
		reference = config['references']['genome_fa'],
		samtools = config['tools']['samtools']
	benchmark:
		"results/{run}/benchmarks/bam/{sample}.{lane}.bwa_mem2.bm"
	# threads: max(1, workflow.cores // max(len(ngs.SAMPLES), 1)) if config['ngs_type'] == 'WES' else max(1, workflow.cores // 2)
	threads:
		4
	shell: """{params.bwa_mem2} mem \
				-o {output.sam} \
				-M -t {threads} -R '@RG\\tID:{wildcards.lane}\\tSM:{wildcards.sample}\\tLB:1\\tPL:ILLUMINA' \
				{params.reference} {input.fr} {input.rr} 2>{log}"""
