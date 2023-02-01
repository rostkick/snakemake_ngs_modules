rule read_alignment:
	input: 
		fr = lambda wc: ngs.long_df.loc[ngs.long_df.loc[:, 'samples']==wc.sample, 'fastq_forward'].values[0],
		rr = lambda wc: ngs.long_df.loc[ngs.long_df.loc[:, 'samples']==wc.sample, 'fastq_reverse'].values[0]
	output: touch('results/{run}/bam/{sample}.raw.bam')
	log: 'results/{run}/logs/aligning/{sample}.bwa_mem2.log'
	params: 
		bwa_mem2=config['tools']['bwa_mem2'],
		reference=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		samtools=config['tools']['samtools']
	threads: workflow.cores/len(SAMPLES)
	shell: """{params.bwa_mem2} mem \
				-M -t {threads} -R '@RG\\tID:lane\\tSM:{wildcards.sample}\\tLB:1\\tPL:ILLUMINA' \
				{params.reference} {input.fr} {input.rr} 2>{log} |\
			{params.samtools} view -bS -o {output} -
			"""
