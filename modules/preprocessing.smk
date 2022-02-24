rule sort_bam:
	input: "results/{run}/bam/{sample}.raw.bam"
	output: "results/{run}/bam/{sample}.sorted.bam"
	threads: workflow.cores/len(SAMPLES)
	shell: 'samtools sort -@ {threads} -o {output} {input}'

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
	shell: 'samtools index {input}'

rule prep_bqsr:
	input: 
		bam="results/{run}/bam/{sample}.dedup.bam",
		bai="results/{run}/bam/{sample}.dedup.bam.bai"
	output: 'results/{run}/bam/{sample}.bqsr.recal.table'
	log: 'results/{run}/logs/prep/{sample}.bqsr_recal.log'
	params: 
		ref=config['GRCh38']['GATK_b38']['reference_fasta'],
		snps=config['GRCh38']['GATK_b38']['snps'],
		indels=config['GRCh38']['GATK_b38']['indels'],
		wgs_calling_regions=config['GRCh38']['GATK_b38']['wgs_calling_regions'] # <-- WGS intervals for exome?
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
		ref=config['GRCh38']['GATK_b38']['reference_fasta'],
		wgs_calling_regions=config['GRCh38']['GATK_b38']['wgs_calling_regions']
	threads: workflow.cores/len(SAMPLES)
	shell: """gatk ApplyBQSR \
				-R {params.ref} \
				-I {input.bam} \
				-L {params.wgs_calling_regions} \
				-bqsr {input.bqsr} \
				-O {output} &>{log}"""

rule CalculateHsMetrics:
	input: "results/{run}/bam/{sample}.final.bam"
	output: "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	log: 'results/{run}/logs/hs_metrics/{sample}.hs_metrics.log'
	params:
		picard = config['tools']['picard_old'],
		fasta_reference = config['GRCh38']['GATK_b38']['reference_fasta'],
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
