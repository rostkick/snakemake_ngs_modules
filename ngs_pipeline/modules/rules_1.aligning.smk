rule r1_0_stage_fastq:
	"""Copy FASTQ from slow NFS to local disk before alignment.
	This decouples NFS read pressure from the CPU-intensive alignment step.
	"""
	wildcard_constraints:
		sample='[\w\d]+',
		lane='[lL0-9M]+'
	input:
		fr = lambda wc: ngs.data.loc[(ngs.data["sample"]==wc.sample) & (ngs.data["lane"]==wc.lane) & (ngs.data["reads_orientation"]=='R1'), 'fastq'].tolist(),
		rr = lambda wc: ngs.data.loc[(ngs.data["sample"]==wc.sample) & (ngs.data["lane"]==wc.lane) & (ngs.data["reads_orientation"]=='R2'), 'fastq'].tolist()
	output:
		fr = temp('results/{run}/staged_fastq/{sample}.{lane}.R1.fastq.gz'),
		rr = temp('results/{run}/staged_fastq/{sample}.{lane}.R2.fastq.gz')
	log:
		'results/{run}/logs/aligning/{sample}.{lane}.stage_fastq.log'
	benchmark:
		'results/{run}/benchmarks/bam/{sample}.{lane}.stage_fastq.bm'
	priority: 10
	threads: 1
	resources:
		mem_mb=100,
		runtime_min=720,
		staging_slots=1
	params:
		lock_dir=".staging_locks",
		max_staged=6,
		wait_sec=30,
		part_ttl_min=120
	script:
		"scripts/stage_fastq_with_lock.py"

rule r1_1_read_alignment:
	wildcard_constraints:
		sample='[\w\d]+',
		lane='[lL0-9M]+'
	input: 
		fr = rules.r1_0_stage_fastq.output.fr,
		rr = rules.r1_0_stage_fastq.output.rr
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
	priority: 10
	threads: 20
	resources:
		mem_mb={'panel': 32000, 'WES': 40000, 'WGS': 48000}.get(config['ngs_type'], 40000),
		runtime_min={'panel': 2880, 'WES': 8640, 'WGS': 17280}.get(config['ngs_type'], 8640),
		nfs_io=0
	shell: """{params.bwa_mem2} mem \
				-o {output.sam} \
				-M -t {threads} -R '@RG\\tID:{wildcards.lane}\\tSM:{wildcards.sample}\\tLB:1\\tPL:ILLUMINA' \
				{params.reference} {input.fr} {input.rr} 2>{log}"""
