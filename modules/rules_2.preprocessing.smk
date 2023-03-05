rule r2_sort_premerged_bams:
	input: 
		bam = rules.r1_read_alignment.output.bam
	output: 
		bam = touch('results/{run}/bam/{sample}.{lane}.for_merge.bam')
	params:
		samtools = config['tools']['samtools']
	threads:
		workflow.cores/(len(SAMPLES)*len(LANES))
	shell: "{params.samtools} sort -@ {threads} -o {output.bam} {input.bam}"

rule r2_merge_bams:
	input: 
		bams = expand("results/{{run}}/bam/{{sample}}.{lane}.for_merge.bam", lane=LANES, allow_missing=True)
	output: 
		bam = touch('results/{run}/bam/{sample}.for_sort2.bam')
	params:
		samtools = config['tools']['samtools']
	shell: "{params.samtools} merge {output.bam} {input.bams}"

rule r2_sorting_bam:
	input: 
		bam = rules.r2_merge_bams.output.bam
	output: 
		bam = temp("results/{run}/bam/{sample}.for_dedup.bam")
	threads: 
		workflow.cores/len(SAMPLES)
	params: 
		samtools = config['tools']['samtools']
	shell: "{params.samtools} sort -@ {threads} -o {output.bam} {input.bam}"

rule r2_mark_duplicates:
	input: 
		bam = rules.r2_sorting_bam.output.bam
	output: 
		bam = temp("results/{run}/bam/{sample}.for_bqsr.bam")
	params:
		gatk = config['tools']['gatk'],
		samtools = config['tools']['samtools']
	log: 
		log1 = "results/{run}/logs/prep/{sample}.dedup.log1",
		log2 = "results/{run}/logs/prep/{sample}.dedup.log2"
	shell: """
		{params.gatk} MarkDuplicates \
					-I {input.bam} -O {output.bam} -M {log.log1} 2>{log.log2} \
					--ASSUME_SORTED true"""
		"{params.samtools} index {output.bam}"

rule r2_prepare_bqsr:
	input: 
		bam = rules.r2_mark_duplicates.output.bam
	output: 
		bqsr = temp("results/{run}/bam/{sample}.for_bqsr.recal.table")
	log: 
		'results/{run}/logs/prep/{sample}.bqsr_recal.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		snps = config['references38']['snps'] if config['assembly'] == 'GRCh38' else config['references37']['snps'],
		indels = config['references38']['indels'] if config['assembly'] == 'GRCh38' else config['references37']['indels'],
		wgs_calling_regions = config['references38']['wgs_calling_regions'] if config['assembly'] == 'GRCh38' else config['references37']['wgs_calling_regions']
	shell: """{params.gatk} BaseRecalibrator \
				-R {params.ref} \
				-I {input.bam} \
				--known-sites {params.snps} \
				--known-sites {params.indels} \
				-L {params.wgs_calling_regions} \
				-O {output} &>{log}"""

rule r2_apply_bqsr:
	input: 
		bam = rules.r2_mark_duplicates.output.bam,
		bqsr = rules.r2_prepare_bqsr.output.bqsr
	output: 
		bam = "results/{run}/bam/{sample}.final.bam"
	log: 
		'results/{run}/logs/prep/{sample}.bqsr_apply.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		wgs_calling_regions = config['references38']['wgs_calling_regions'] if config['assembly'] == 'GRCh38' else config['references37']['wgs_calling_regions']
	threads: workflow.cores/len(SAMPLES)
	shell: """{params.gatk} ApplyBQSR \
				-R {params.ref} \
				-I {input.bam} \
				-L {params.wgs_calling_regions} \
				-bqsr {input.bqsr} \
				-O {output} &>{log}"""

rule r2_bed_to_intervals:
	input: 
		bed = config["panel_capture"]["target"]
	output: 
		intervals = "results/{run}/capture.intervals"
	params:
		picard_old = config['tools']['picard_old'],
		ref_dict = config['references38']['dict'] if config['assembly'] == 'GRCh38' else config['references37']['dict'],
	shell: "java -jar {params.picard_old} BedToIntervalList I={input.bed} SD={params.ref_dict} O={output.intervals}"
