rule r3_calculateHsMetrics:
	input: 
		bam = rules.r2_apply_bqsr.output.bam,
		intervals = rules.r2_bed_to_intervals.output.intervals if config['panel_capture']['target'][-4:] == 'bed' else config['panel_capture']['target']
	output: 
		tsv = "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	log: 'results/{run}/logs/hs_metrics/{sample}.hs_metrics.log'
	params:
		picard = config['tools']['picard_old'],
		fasta_reference = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa']
	shell: """java -jar {params.picard} CalculateHsMetrics \
				INPUT={input.bam} \
				OUTPUT={output.tsv} \
				BAIT_INTERVALS={input.intervals} \
				TARGET_INTERVALS={input.intervals} \
				NEAR_DISTANCE=0 2>{log}"""
