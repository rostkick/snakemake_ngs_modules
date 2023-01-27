rule BedToIntervalList:
	input: config["panel_capture"]["target"]
	output: "results/{run}/capture.intervals"
	params:
		picard_old=config['tools']['picard_old'],
		ref_dict=config['references']['dict']
	shell: "java -jar {params.picard_old} BedToIntervalList I={input} SD={params.ref_dict} O={output}"

rule CalculateHsMetrics:
	input: 
		bam="results/{run}/bam/{sample}.final.bam",
		intervals="results/{run}/capture.intervals"
	output: "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	log: 'results/{run}/logs/hs_metrics/{sample}.hs_metrics.log'
	params:
		picard = config['tools']['picard_old'],
		fasta_reference = config['references']['genome_fa']
	shell: """java -jar {params.picard} CalculateHsMetrics \
				INPUT={input.bam} \
				OUTPUT={output} \
				BAIT_INTERVALS={input.intervals} \
				TARGET_INTERVALS={input.intervals} \
				NEAR_DISTANCE=0 2>{log}"""

rule hs_metrics_status:
	input: "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	output: touch("results/{run}/bam/hs_metrics/{sample}.metrics_status")