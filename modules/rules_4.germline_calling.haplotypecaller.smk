rule r4_haplotypecaller:
	input: 
		bam = rules.r2_apply_bqsr.output.bam
	output:
		gvcf = "results/{run}/germline/vcf/{sample}.gvcf.gz"
	log: 'results/{run}/logs/germline/{sample}.haplotypecaller.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		panel_capture = config['panel_capture']['target']
	threads: 8
	shell: """{params.gatk} --java-options "-Xmx{threads}g -XX:ParallelGCThreads=1" \
		HaplotypeCaller \
		-R {params.ref} \
		-I {input.bam} \
		-O {output.gvcf} \
		-L {params.panel_capture} \
		-ERC GVCF \
		--native-pair-hmm-threads {threads} 2>{log} || true"""

rule r4_haplotypecaller_wgs:
	input: 
		bam = rules.r2_apply_bqsr.output.bam
	output:
		gvcf = "results/{run}/germline/vcf/{sample}.wgs.gvcf.gz"
	log: 'results/{run}/logs/germline/{sample}.haplotypecaller.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
	threads: 8
	shell: """{params.gatk} --java-options "-Xmx{threads}g -XX:ParallelGCThreads=1" \
		HaplotypeCaller \
		-R {params.ref} \
		-I {input.bam} \
		-O {output.gvcf} \
		-ERC GVCF \
		--native-pair-hmm-threads {threads} 2>{log} || true"""

rule r4_combinegvcfs:
	input: 
		gvcfs = expand("results/{{run}}/germline/vcf/{sample}.gvcf.gz", sample=ngs.GRM_SAMPLES) \
					if config['ngs_type'] != 'WGS' else \
						expand("results/{{run}}/germline/vcf/{sample}.wgs.gvcf.gz", sample=ngs.GRM_SAMPLES)
	output: 
		gvcf = "results/{run}/germline/vcf/cohort.g.vcf.gz"
	log: 'results/{run}/logs/germline/combinegvcfs.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		gvcf=lambda wc, input: [f'--variant {i}' for i in input]
	shell:"""{params.gatk} CombineGVCFs \
		-R {params.ref} \
		{params.gvcf} \
		-O {output.gvcf} 2>{log}"""

rule r4_genotypegvcfs:
	input: 
		gvcf = rules.r4_combinegvcfs.output.gvcf
	output: 
		vcf = "results/{run}/germline/vcf/cohort.to_filters.vcf.gz"
	log: 'results/{run}/logs/germline/genotypegvcfs.log'
	params:
		gatk = config['tools']['gatk'],
		fasta_reference=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		panel_capture = config['panel_capture']['target']
	shell:"""{params.gatk} GenotypeGVCFs \
		-R {params.fasta_reference} \
		-L {params.panel_capture} \
		-V {input.gvcf} \
		-O {output.vcf} 2>{log}"""

rule r4_genotypegvcfs_wgs:
	input: 
		gvcf = rules.r4_combinegvcfs.output.gvcf
	output: 
		vcf = "results/{run}/germline/vcf/cohort.to_filters.wgs.vcf.gz"
	log: 'results/{run}/logs/germline/genotypegvcfs.log'
	params:
		gatk = config['tools']['gatk'],
		fasta_reference=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa']
	shell:"""{params.gatk} GenotypeGVCFs \
		-R {params.fasta_reference} \
		-V {input.gvcf} \
		-O {output.vcf} 2>{log}"""


rule r4_selectvariants_snp:
	input:
		vcf = rules.r4_genotypegvcfs.output.vcf if config['ngs_type'] != 'WGS' else rules.r4_genotypegvcfs_wgs.output.vcf
	output: 
		vcf = "results/{run}/germline/vcf/snps.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	log: 'results/{run}/logs/germline/select_snp.log'
	shell: """{params.gatk} SelectVariants \
		-V {input.vcf} \
		-select-type SNP \
		-O {output.vcf} 2>{log}"""

rule r4_selectvariants_indel:
	input:
		vcf = rules.r4_genotypegvcfs.output.vcf if config['ngs_type'] != 'WGS' else rules.r4_genotypegvcfs_wgs.output.vcf
	output:
		vcf = "results/{run}/germline/vcf/indel.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	log: 'results/{run}/logs/germline/select_indel.log'
	shell: """{params.gatk} SelectVariants \
		-V {input.vcf} \
		-select-type INDEL \
		-O {output.vcf} 2>{log}"""

rule r4_variantfiltration_snp:
	input: 
		vcf = rules.r4_selectvariants_snp.output.vcf
	output: 
		vcf = "results/{run}/germline/vcf/snps.filtered.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	log: 'results/{run}/logs/germline/variantfilttation_snp.log'
	shell: """{params.gatk} VariantFiltration \
		-V {input.vcf} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "SOR > 3.0" --filter-name "SOR3" \
		-filter "FS > 60.0" --filter-name "FS60" \
		-filter "MQ < 40.0" --filter-name "MQ40" \
		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
		-O {output.vcf} 2>{log}"""

rule r4_variantfiltration_indel:
	input: 
		vcf = rules.r4_selectvariants_indel.output.vcf
	output: 
		vcf = "results/{run}/germline/vcf/indel.filtered.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	log: "results/{run}/logs/germline/variantfiltration_indel.log"
	shell:"""{params.gatk} VariantFiltration \
		-V {input.vcf} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "FS > 200.0" --filter-name "FS200" \
		-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
		-O {output.vcf} 2>{log}"""

rule r4_mergevcfs:
	input: 
		snps = rules.r4_variantfiltration_snp.output.vcf,
		indels = rules.r4_variantfiltration_indel.output.vcf
	output: 
		vcf = "results/{run}/germline/vcf/cohort.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	log: "results/{run}/logs/germline/mergevcfs.log"
	shell:"""{params.gatk} MergeVcfs \
		-I {input.snps} \
		-I {input.indels} \
		-O {output.vcf} 2>{log}"""
