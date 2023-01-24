rule CalculateHsMetrics:
	input: "results/{run}/bam/{sample}.final.bam"
	output: "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	log: 'results/{run}/logs/hs_metrics/{sample}.hs_metrics.log'
	params:
		picard = config['tools']['picard_old'],
		fasta_reference = config['references']['genome_fa'],
		intervals = config['panel_capture']['intervals']
	shell: """java -jar {params.picard} CalculateHsMetrics \
				INPUT={input} \
				OUTPUT={output} \
				BAIT_INTERVALS={params.intervals} \
				TARGET_INTERVALS={params.intervals} \
				NEAR_DISTANCE=0 2>{log}"""

rule hs_metrics_status:
	input: "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	output: touch("results/{run}/bam/hs_metrics/{sample}.metrics_status")