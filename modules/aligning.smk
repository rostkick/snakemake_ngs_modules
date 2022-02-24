rule align:
	input: 
		fr = lambda wc: long_data.loc[long_data.loc[:, 'samples']==wc.sample, 'fastq_forward'].values[0],
		rr = lambda wc: long_data.loc[long_data.loc[:, 'samples']==wc.sample, 'fastq_reverse'].values[0]
	output: touch('results/{run}/bam/{sample}.raw.bam')
	log: 'results/{run}/logs/aligning/{sample}.bwa_mem2.log'
	params: reference=config['GRCh38']['GATK_b38']['reference_fasta']
	threads: workflow.cores/len(SAMPLES)
	shell: """bwa-mem2 mem \
				-M -t {threads} -R '@RG\\tID:lane\\tSM:{wildcards.sample}\\tLB:1\\tPL:ILLUMINA' \
				{params.reference} {input.fr} {input.rr} 2>{log} |\
			samtools view -bS -o {output} -
			"""
