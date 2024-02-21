rule r1_read_alignment:
	wildcard_constraints:
		sample='[\w\d]+',
		lane='[lL0-9M]+'
	input: 
		fr = lambda wc: ngs.data.loc[(ngs.data["sample"]==wc.sample) & (ngs.data["lane"]==wc.lane) & (ngs.data["reads_orientation"]=='R1'), 'fastq'].tolist(),
		rr = lambda wc: ngs.data.loc[(ngs.data["sample"]==wc.sample) & (ngs.data["lane"]==wc.lane) & (ngs.data["reads_orientation"]=='R2'), 'fastq'].tolist()
	output: 
		sam = 'results/{run}/bam/{sample}.{lane}.sam'
	log: 
		'results/{run}/logs/aligning/{sample}.{lane}.bwa_mem2.log'
	params: 
		bwa_mem2 = config['tools']['bwa_mem2'],
		reference = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		samtools = config['tools']['samtools']
	threads: workflow.cores/len(ngs.SAMPLES) if config['ngs_type'] == 'WES' else workflow.cores/2
	shell: """{params.bwa_mem2} mem \
				-o {output.sam} \
				-M -t {threads} -R '@RG\\tID:{wildcards.lane}\\tSM:{wildcards.sample}\\tLB:1\\tPL:ILLUMINA' \
				{params.reference} {input.fr} {input.rr} 2>{log}"""
