rule r1_read_alignment:
	input: 
		fr = lambda wc: ngs.data.loc[(ngs.data["sample"]==wc.sample[:-4]) & (ngs.data["lane"]==wc.lane) & (ngs.data["reads_orientation"]=='R1'), 'fastq'].tolist()[0],
		rr = lambda wc: ngs.data.loc[(ngs.data["sample"]==wc.sample[:-4]) & (ngs.data["lane"]==wc.lane) & (ngs.data["reads_orientation"]=='R2'), 'fastq'].tolist()[0]
	output: 
		bam = 'results/{run}/bam/{sample}.{lane}.for_sort1.bam'
	log: 
		'results/{run}/logs/aligning/{sample}.{lane}.bwa_mem2.log'
	params: 
		bwa_mem2 = config['tools']['bwa_mem2'],
		reference = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		samtools = config['tools']['samtools']
	threads: workflow.cores/len(ngs.SAMPLES)
	shell: """{params.bwa_mem2} mem \
				-M -t {threads} -R '@RG\\tID:{wildcards.lane}\\tSM:{wildcards.sample}\\tLB:1\\tPL:ILLUMINA' \
				{params.reference} {input.fr} {input.rr} 2>{log} |\
			{params.samtools} view -bS -o {output.bam} -
			"""
