rule r3_1_calculate_hs_metrics:
	input: 
		bam = rules.r2_9_apply_bqsr.output.bam,
		intervals = rules.r2_1_bed_to_intervals.output.intervals if config['panel_capture']['target'][-4:] == '.bed' else config['panel_capture']['target']
	output: 
		tsv = "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	benchmark:
		'results/{run}/benchmarks/bam/hs_metrics/{sample}.hs_metrics.bm'
	log: 'results/{run}/logs/hs_metrics/{sample}.hs_metrics.log'
	params:
		picard = config['tools']['picard_old'],
		fasta_reference = config['references']['genome_fa']
	shell: """java -jar {params.picard} CalculateHsMetrics \
				INPUT={input.bam} \
				OUTPUT={output.tsv} \
				BAIT_INTERVALS={input.intervals} \
				TARGET_INTERVALS={input.intervals} \
				NEAR_DISTANCE=0 2>{log}"""