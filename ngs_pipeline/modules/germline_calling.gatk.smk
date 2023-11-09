rule HaplotypeCaller:
	input: 
		bam='results/{run}/bam/{sample}.final.bam',
		metrics="results/{run}/bam/hs_metrics/{sample}.metrics_status"
	output: "results/{run}/germline/vcf/{sample}.g.vcf.gz"
	log: 'results/{run}/logs/germline_calling/{sample}.haplotypecaller.log'
	params:
		fasta_reference=config['GRCh38']['GATK_b38']['reference_fasta'],
		intervals=config['panel_capture']['target']
	threads: workflow.cores/8
	shell: """gatk --java-options "-Xmx{threads}g -XX:ParallelGCThreads=1" \
		HaplotypeCaller \
		-R {params.fasta_reference} \
		-I {input.bam} \
		-O {output} \
		-ERC GVCF \
		-L {params.intervals} \
		--native-pair-hmm-threads {threads} 2>{log}"""

rule CombineGVCFs:
	input: expand("results/{{run}}/germline/vcf/{sample}.g.vcf.gz", sample=GERMLINE_SAMPLES)
	output: "results/{run}/germline/vcf/cohort.g.vcf.gz"
	log: 'results/{run}/logs/germline_calling/combinegvcfs.log'
	params:
		fasta_reference=config['GRCh38']['GATK_b38']['reference_fasta'],
		gvcf=lambda wc, input: [f'--variant {i}' for i in input]
	shell:"""gatk CombineGVCFs \
		-R {params.fasta_reference} \
		{params.gvcf} \
		-O {output} 2>{log}"""

rule GenotypeGVCFs:
	input: "results/{run}/germline/vcf/cohort.g.vcf.gz"
	output: "results/{run}/germline/vcf/cohort.vcf.gz"
	log: 'results/{run}/logs/germline_calling/genotypegvcfs.log'
	params:
		fasta_reference=config['GRCh38']['GATK_b38']['reference_fasta'],
		intervals=config['panel_capture']['target']
	shell:"""gatk GenotypeGVCFs \
		-R {params.fasta_reference} \
		-V {input} \
		-L {params.intervals} \
		-O {output} 2>{log}"""

rule SelectVariants_SNP:
	input: "results/{run}/germline/vcf/cohort.vcf.gz"
	output: "results/{run}/germline/vcf/snps.vcf.gz"
	log: 'results/{run}/logs/germline_calling/select_snp.log'
	shell: """gatk SelectVariants \
		-V {input} \
		-select-type SNP \
		-O {output} 2>{log}"""

rule SelectVariants_INDEL:
	input: "results/{run}/germline/vcf/cohort.vcf.gz"
	output: "results/{run}/germline/vcf/indel.vcf.gz"
	log: 'results/{run}/logs/germline_calling/select_indel.log'
	shell: """gatk SelectVariants \
		-V {input} \
		-select-type INDEL \
		-O {output} 2>{log}"""

rule VariantFiltration_SNP:
	input: "results/{run}/germline/vcf/snps.vcf.gz"
	output: "results/{run}/germline/vcf/snps.filtered.vcf.gz"
	log: 'results/{run}/logs/germline_calling/variantfilttation_snp.log'
	shell: """gatk VariantFiltration \
		-V {input} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "SOR > 3.0" --filter-name "SOR3" \
		-filter "FS > 60.0" --filter-name "FS60" \
		-filter "MQ < 40.0" --filter-name "MQ40" \
		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
		-O {output} 2>{log}"""

rule VariantFiltration_INDEL:
	input: "results/{run}/germline/vcf/indel.vcf.gz"
	output: "results/{run}/germline/vcf/indel.filtered.vcf.gz"
	log: 'results/{run}/logs/germline_calling/variantfiltration_indel.log'
	shell:"""gatk VariantFiltration \
		-V {input} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "FS > 200.0" --filter-name "FS200" \
		-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
		-O {output} 2>{log}"""

rule MergeVcfs:
	input: 
		snps="results/{run}/germline/vcf/snps.filtered.vcf.gz",
		indels="results/{run}/germline/vcf/indel.filtered.vcf.gz"
	output: "results/{run}/germline/vcf/cohort.filtered.vcf.gz"
	log: 'results/{run}/logs/germline_calling/mergevcfs.log'
	shell:"""gatk MergeVcfs \
		-I {input.snps} \
		-I {input.indels} \
		-O {output} 2>{log}"""
