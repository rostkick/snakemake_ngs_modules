def input_align(wc) -> str:
	fastq = ngs.data.loc[(ngs.data["sample"]==wc.sample) & (ngs.data["lane"]==wc.lane), 'fastq'].sort_values().tolist()
	fastq = ' '.join(fastq)
	return fastq

rule r1_read_alignment:
	wildcard_constraints:
		sample='[\w\d]+',
		lane='[lL0-9M]+'
	input: 
		fastq = lambda wc: ngs.data.loc[(ngs.data["sample"]==wc.sample) & (ngs.data["lane"]==wc.lane) & (ngs.data["reads_orientation"]=='R1'), 'fastq'].tolist()
	output: 
		sam = 'results/{run}/bam/{sample}.{lane}.sam'
	log: 
		'results/{run}/logs/aligning/{sample}.{lane}.bwa_mem2.log'
	params:
		fastq = lambda wc: input_align(wc),
		bwa_mem2 = config['tools']['bwa_mem2'],
		reference = config['references']['genome_fa'],
		samtools = config['tools']['samtools']
	threads: workflow.cores/len(ngs.SAMPLES) if config['ngs_type'] == 'WES' else workflow.cores/2
	shell: """{params.bwa_mem2} mem \
				-o {output.sam} \
				-M -t {threads} -R '@RG\\tID:{wildcards.lane}\\tSM:{wildcards.sample}\\tLB:1\\tPL:ILLUMINA' \
				{params.reference} {params.fastq} 2>{log}"""
