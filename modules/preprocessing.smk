rule sort_bam:
	input: "results/{run}/bam/{sample}.raw.bam"
	output: "results/{run}/bam/{sample}.sorted.bam"
	threads: workflow.cores/len(SAMPLES)
	params: samtools=config['tools']['samtools']
	shell: '{params.samtools} sort -@ {threads} -o {output} {input}'

rule dedup:
	input: "results/{run}/bam/{sample}.sorted.bam"
	output: "results/{run}/bam/{sample}.dedup.bam"
	log: 
		log1='results/{run}/logs/prep/{sample}.dedup.log1',
		log2='results/{run}/logs/prep/{sample}.dedup.log2'
	shell: """gatk MarkDuplicates \
				-I {input} -O {output} -M {log.log1} 2>{log.log2}\
				--ASSUME_SORTED true"""

rule index_dedup:
	input: "results/{run}/bam/{sample}.dedup.bam"
	output: "results/{run}/bam/{sample}.dedup.bam.bai"
	params: samtools=config['tools']['samtools']
	shell: '{params.samtools} index {input}'

rule prep_bqsr:
	input: 
		bam="results/{run}/bam/{sample}.dedup.bam",
		bai="results/{run}/bam/{sample}.dedup.bam.bai"
	output: 'results/{run}/bam/{sample}.bqsr.recal.table'
	log: 'results/{run}/logs/prep/{sample}.bqsr_recal.log'
	params: 
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		snps=config['references38']['snps'] if config['assembly'] == 'GRCh38' else config['references37']['snps'],
		indels=config['references38']['indels'] if config['assembly'] == 'GRCh38' else config['references37']['indels'],
		wgs_calling_regions=config['references38']['wgs_calling_regions'] if config['assembly'] == 'GRCh38' else config['references37']['wgs_calling_regions']
	shell: """gatk BaseRecalibrator \
				-R {params.ref} \
				-I {input.bam} \
				--known-sites {params.snps} \
				--known-sites {params.indels} \
				-L {params.wgs_calling_regions} \
				-O {output} &>{log}"""

rule apply_bqsr:
	input: 
		bam="results/{run}/bam/{sample}.dedup.bam",
		bqsr='results/{run}/bam/{sample}.bqsr.recal.table'
	output: "results/{run}/bam/{sample}.final.bam"
	log: 'results/{run}/logs/prep/{sample}.bqsr_apply.log'
	params: 
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		wgs_calling_regions=config['references38']['wgs_calling_regions'] if config['assembly'] == 'GRCh38' else config['references37']['wgs_calling_regions']
	threads: workflow.cores/len(SAMPLES)
	shell: """gatk ApplyBQSR \
				-R {params.ref} \
				-I {input.bam} \
				-L {params.wgs_calling_regions} \
				-bqsr {input.bqsr} \
				-O {output} &>{log}"""

rule BedToIntervalList:
	input: config["panel_capture"]["target"]
	output: "results/{run}/capture.intervals"
	params:
		picard_old=config['tools']['picard_old'],
		ref_dict=config['references38']['dict'] if config['assembly'] == 'GRCh38' else config['references37']['dict'],
	shell: "java -jar {params.picard_old} BedToIntervalList I={input} SD={params.ref_dict} O={output}"
