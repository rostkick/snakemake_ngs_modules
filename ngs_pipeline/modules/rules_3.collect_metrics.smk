rule r3_1_calculate_hs_metrics:
	input: 
		bam = rules.r2_9_apply_bqsr.output.bam,
		intervals = rules.r2_1_bed_to_intervals.output.intervals if config['panel_capture']['target'][-4:] == '.bed' else config['panel_capture']['target']
	output: 
		tsv = "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	benchmark:
		'results/{run}/benchmarks/bam/hs_metrics/{sample}.hs_metrics.bm'
	log: 'results/{run}/logs/hs_metrics/{sample}.hs_metrics.log'
	priority: 20
	resources:
		mem_mb=4000,
		runtime_min=240
	params:
		picard = config['tools']['picard_old'],
		fasta_reference = config['references']['genome_fa'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	shell: """java -jar {params.java_opts} {params.picard} CalculateHsMetrics \
				INPUT={input.bam} \
				OUTPUT={output.tsv} \
				REFERENCE_SEQUENCE={params.fasta_reference} \
				BAIT_INTERVALS={input.intervals} \
				TARGET_INTERVALS={input.intervals} \
				NEAR_DISTANCE=0 2>{log}"""