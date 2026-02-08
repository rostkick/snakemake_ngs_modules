rule r4_1_haplotypecaller:
	input: 
		bam = rules.r2_7_apply_bqsr.output.bam
	output:
		gvcf = "results/{run}/germline/vcf/{sample}.gatk.gvcf.gz"
	benchmark:
		'results/{run}/benchmarks/germline/vcf/{sample}.haplotypecaller.bm'
	log: 'results/{run}/logs/germline/{sample}.haplotypecaller.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		panel_capture = config['panel_capture']['target']
	threads: 8
	resources:
		mem_mb=16000
	shell: """
		{params.gatk} --java-options "-Xmx{threads}g -XX:ParallelGCThreads=1" \
			HaplotypeCaller \
			-R {params.ref} \
			-I {input.bam} \
			-O {output.gvcf} \
			-L {params.panel_capture} \
			-ERC GVCF \
			--native-pair-hmm-threads {threads} 2>{log}
	"""

rule r4_1_haplotypecaller_wgs:
	input: 
		bam = rules.r2_7_apply_bqsr.output.bam
	output:
		gvcf = "results/{run}/germline/vcf/{sample}.wgs.gatk.gvcf.gz"
	log: 'results/{run}/logs/germline/{sample}.haplotypecaller_wgs.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa']
	threads: 8
	resources:
		mem_mb=16000
	shell: """
		{params.gatk} --java-options "-Xmx{threads}g -XX:ParallelGCThreads=1" \
			HaplotypeCaller \
			-R {params.ref} \
			-I {input.bam} \
			-O {output.gvcf} \
			-ERC GVCF \
			--native-pair-hmm-threads {threads} 2>{log} || true
	"""

rule r4_2_combinegvcfs:
	input: 
		gvcfs = expand("results/{{run}}/germline/vcf/{sample}.gatk.gvcf.gz", sample=ngs.GRM_SAMPLES) \
					if config['ngs_type'] != 'WGS' else \
						expand("results/{{run}}/germline/vcf/{sample}.wgs.gatk.gvcf.gz", sample=ngs.GRM_SAMPLES)
	output: 
		gvcf = "results/{run}/germline/vcf/cohort.g.vcf.gz"
	benchmark:
		'results/{run}/benchmarks/germline/vcf/combine_gvcf.bm'
	log: 'results/{run}/logs/germline/combinegvcfs.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		gvcf=lambda wc, input: [f'--variant {i}' for i in input]
	threads: 1
	resources:
		mem_mb=8000
	shell: """
		{params.gatk} CombineGVCFs \
			-R {params.ref} \
			{params.gvcf} \
			-O {output.gvcf} 2>{log}
	"""

rule r4_3_genotypegvcfs:
	input: 
		gvcf = rules.r4_2_combinegvcfs.output.gvcf
	output: 
		vcf = "results/{run}/germline/vcf/cohort.to_filters.vcf.gz"
	benchmark:
		'results/{run}/benchmarks/germline/vcf/genotype_gvcf.bm'
	log: 'results/{run}/logs/germline/genotypegvcfs.log'
	params:
		gatk = config['tools']['gatk'],
		fasta_reference=config['references']['genome_fa'],
		panel_capture = config['panel_capture']['target']
	threads: 1
	resources:
		mem_mb=8000
	shell: """
		{params.gatk} GenotypeGVCFs \
			-R {params.fasta_reference} \
			-L {params.panel_capture} \
			-V {input.gvcf} \
			-O {output.vcf} 2>{log}
	"""

rule r4_3_genotypegvcfs_wgs:
	input: 
		gvcf = rules.r4_2_combinegvcfs.output.gvcf
	output: 
		vcf = "results/{run}/germline/vcf/cohort.to_filters.wgs.vcf.gz"
	log: 'results/{run}/logs/germline/genotypegvcfs_wgs.log'
	params:
		gatk = config['tools']['gatk'],
		fasta_reference=config['references']['genome_fa']
	threads: 1
	resources:
		mem_mb=8000
	shell: """
		{params.gatk} GenotypeGVCFs \
			-R {params.fasta_reference} \
			-V {input.gvcf} \
			-O {output.vcf} 2>{log}
	"""

rule r4_4_selectvariants_snp:
	input:
		vcf = rules.r4_3_genotypegvcfs.output.vcf if config['ngs_type'] != 'WGS' else rules.r4_3_genotypegvcfs_wgs.output.vcf
	output: 
		vcf = "results/{run}/germline/vcf/snps.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	benchmark:
		'results/{run}/benchmarks/germline/vcf/select_snp.bm'
	log: 'results/{run}/logs/germline/select_snp.log'
	threads: 1
	resources:
		mem_mb=4000
	shell: """
		{params.gatk} SelectVariants \
			-V {input.vcf} \
			-select-type SNP \
			-O {output.vcf} 2>{log}
	"""

rule r4_5_selectvariants_indel:
	input:
		vcf = rules.r4_3_genotypegvcfs.output.vcf if config['ngs_type'] != 'WGS' else rules.r4_3_genotypegvcfs_wgs.output.vcf
	output:
		vcf = "results/{run}/germline/vcf/indel.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	log: 'results/{run}/logs/germline/select_indel.log'
	benchmark:
		'results/{run}/benchmarks/germline/vcf/select_indels.bm'
	threads: 1
	resources:
		mem_mb=4000
	shell: """
		{params.gatk} SelectVariants \
			-V {input.vcf} \
			-select-type INDEL \
			-O {output.vcf} 2>{log}
	"""

rule r4_6_variantfiltration_snp:
	input: 
		vcf = rules.r4_4_selectvariants_snp.output.vcf
	output: 
		vcf = "results/{run}/germline/vcf/snps.filtered.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	benchmark:
		'results/{run}/benchmarks/germline/vcf/apply_hard_filters_snp.bm'
	log: 'results/{run}/logs/germline/variantfilttation_snp.log'
	threads: 1
	resources:
		mem_mb=4000
	shell: """
		{params.gatk} VariantFiltration \
			-V {input.vcf} \
			-filter "QD < 2.0" --filter-name "QD2" \
			-filter "QUAL < 30.0" --filter-name "QUAL30" \
			-filter "SOR > 3.0" --filter-name "SOR3" \
			-filter "FS > 60.0" --filter-name "FS60" \
			-filter "MQ < 40.0" --filter-name "MQ40" \
			-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
			-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
			-O {output.vcf} 2>{log}
	"""

rule r4_7_variantfiltration_indel:
	input:
		vcf = rules.r4_5_selectvariants_indel.output.vcf
	output: 
		vcf = "results/{run}/germline/vcf/indel.filtered.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	benchmark:
		'results/{run}/benchmarks/germline/vcf/apply_hard_filters_indels.bm'
	log: "results/{run}/logs/germline/variantfiltration_indel.log"
	threads: 1
	resources:
		mem_mb=4000
	shell: """
		{params.gatk} VariantFiltration \
			-V {input.vcf} \
			-filter "QD < 2.0" --filter-name "QD2" \
			-filter "QUAL < 30.0" --filter-name "QUAL30" \
			-filter "FS > 200.0" --filter-name "FS200" \
			-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
			-O {output.vcf} 2>{log}
	"""

rule r4_8_mergevcfs:
	input: 
		snps = rules.r4_6_variantfiltration_snp.output.vcf,
		indels = rules.r4_7_variantfiltration_indel.output.vcf
	output: 
		vcf = "results/{run}/germline/vcf/cohort.vcf.gz"
	params:
		gatk = config['tools']['gatk']
	benchmark:
		'results/{run}/benchmarks/germline/vcf/merge_vcfs.bm'
	log: "results/{run}/logs/germline/mergevcfs.log"
	threads: 1
	resources:
		mem_mb=4000
	shell: """
		{params.gatk} MergeVcfs \
			-I {input.snps} \
			-I {input.indels} \
			-O {output.vcf} 2>{log}
	"""

rule r4_9_deepvariant:
	input: 
		bam = rules.r2_7_apply_bqsr.output.bam,
		bed = rules.r2_0_prepare_bed_to_bed.output.bed
	output: 
		vcf = "results/{run}/germline/vcf/{sample}.vcf.gz",
		gvcf = "results/{run}/germline/vcf/{sample}.gvcf.gz"
	benchmark:
		'results/{run}/benchmarks/germline/vcf/{sample}.dv_calling.bm'
	log:
		'results/{run}/logs/germline/{sample}.deepvariant.log'
	params:
		singularity = config['tools']['singularity'],
		deepvariant = config['tools']['deepvariant'],
		ngs_type = config['ngs_type'],
		ref = config['references']['genome_fa']
	threads: 16
	resources:
		mem_mb=20000
	shell: """
		{params.singularity} run \
			-B /ngs_pipeline:/ngs_pipeline \
			{params.deepvariant} /opt/deepvariant/bin/run_deepvariant \
			--model_type={params.ngs_type} \
			--output_vcf={output.vcf} \
			--output_gvcf={output.gvcf} \
			--reads={input.bam} \
			--ref={params.ref} \
			--regions {input.bed} \
			--num_shards={threads} &>{log} || true
	"""

rule r4_10_archive_vcfs:
	input:
		vcf = rules.r4_9_deepvariant.output.vcf,
		gvcf = rules.r4_9_deepvariant.output.gvcf
	output:
		vcf = "archive/{run}/germline/vcf/{sample}.vcf.gz",
		vcf_tbi = "archive/{run}/germline/vcf/{sample}.vcf.gz.tbi",
		gvcf = "archive/{run}/germline/vcf/{sample}.gvcf.gz",
		gvcf_tbi = "archive/{run}/germline/vcf/{sample}.gvcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=1000
	benchmark:
		'results/{run}/benchmarks/germline/vcf/{sample}.archive_vcf.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_vcf.log"
	shell: """
		mkdir -p $(dirname {output.vcf})
		rsync -av {input.vcf}* $(dirname {output.vcf})/ &>>{log}
		rsync -av {input.gvcf}* $(dirname {output.gvcf})/ &>>{log}
		echo "Archived {wildcards.sample} (VCF + GVCF with indices) to network storage" >> {log}
	"""