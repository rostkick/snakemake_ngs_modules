rule r2_1_sam_to_bam:
	input: 
		sam = rules.r1_read_alignment.output.sam
	output:
		bam = temp('results/{run}/bam/{sample}.{lane}.for_sort1.bam')
	params:
		samtools = config['tools']['samtools']
	threads:
		workflow.cores/(len(ngs.SAMPLES) * max(len(ngs.LANES), 1)) if config['ngs_type'] == 'WES' else workflow.cores/2
	shell: "{params.samtools} view -@ {threads} -bS -o {output.bam} {input.sam}"

rule r2_2_sort_premerged_bams:
	input: 
		bam = rules.r2_1_sam_to_bam.output.bam
	output: 
		bam = temp('results/{run}/bam/{sample}.{lane}.for_merge.bam')
	params:
		samtools = config['tools']['samtools']
	threads:
		workflow.cores/(len(ngs.SAMPLES)*max(len(ngs.LANES), 1))
	shell: "{params.samtools} sort -@ {threads} -o {output.bam} {input.bam}"

# Function to create filtered combinations of wildcards, based on the presence of input files.
ngs_comb = set(ngs.data.loc[:, ['sample', 'lane']].itertuples(index=False, name=None))
allow_combs = []
for s, l in ngs_comb:
	tup = (("sample", s), ("lane", l))
	allow_combs.append(tup)

def filter_combinator(whitelist):
	def filtered_combinator(*args, **kwargs):
		for wc_comb in product(*args, **kwargs):
			for ac in allow_combs:
				if(wc_comb[0:3] == ac):
					print ("SUCCESS")
					yield(wc_comb)
					break
	return filtered_combinator
filtered_product = filter_combinator(allow_combs)

def get_merged_bams(wc):
	result_output = []
	for run, sample, lane in product([config['run']], [wc.sample], ngs.LANES):
		line = ngs.data.loc[(ngs.data["sample"]==sample) & (ngs.data["lane"]==lane), 'fastq'].tolist()
		if line:
			result_output.append(f'results/{run}/bam/{sample}.{lane}.for_merge.bam')
	return result_output

rule r2_3_merge_bams:
	input:
		bams = get_merged_bams
	output: 
		bam = temp('results/{run}/bam/{sample}.for_sort2.bam')
	params:
		samtools = config['tools']['samtools']
	run:
		if len(input.bams)>1:
			print(f"{params.samtools} merge {output.bam} {input.bams}")
			shell("{params.samtools} merge {output.bam} {input.bams}")
		else:
			print(f"ln -s {input.bams} {output.bam}")
			shell("cp {input.bams} {output.bam}")

rule r2_4_sorting_bam:
	input: 
		bam = rules.r2_3_merge_bams.output.bam
	output: 
		bam = temp("results/{run}/bam/{sample}.for_dedup.bam")
	threads: 
		workflow.cores/len(ngs.SAMPLES)
	params: 
		samtools = config['tools']['samtools']
	shell: "{params.samtools} sort -@ {threads} -o {output.bam} {input.bam}"

rule r2_5_mark_duplicates:
	input: 
		bam = rules.r2_4_sorting_bam.output.bam
	output: 
		bam = temp("results/{run}/bam/{sample}.for_bqsr.bam")
	params:
		gatk = config['tools']['gatk'],
		samtools = config['tools']['samtools']
	log: 
		log1 = "results/{run}/logs/prep/{sample}.dedup.log1",
		log2 = "results/{run}/logs/prep/{sample}.dedup.log2"
	run: 
		if config['hs']:
			shell("""
			{params.gatk} MarkDuplicates \
						-I {input.bam} -O {output.bam} -M {log.log1} 2>{log.log2} \
						--ASSUME_SORTED true
						{params.samtools} index {output.bam}""")
		else:
			shell('ln -s {input.bam} {output.bam} && {params.samtools} index {output.bam}')

rule r2_6_prepare_bqsr:
	input: 
		bam = rules.r2_5_mark_duplicates.output.bam
	output: 
		bqsr = temp("results/{run}/bam/{sample}.for_bqsr.recal.table")
	log: 
		'results/{run}/logs/prep/{sample}.bqsr_recal.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		snps = config['references']['snps'],
		indels = config['references']['indels'],
		wgs_calling_regions = config['references']['wgs_calling_regions']
	shell: """{params.gatk} BaseRecalibrator \
				-R {params.ref} \
				-I {input.bam} \
				--known-sites {params.snps} \
				--known-sites {params.indels} \
				-L {params.wgs_calling_regions} \
				-O {output} &>{log}"""

rule r2_7_apply_bqsr:
	input: 
		bam = rules.r2_5_mark_duplicates.output.bam,
		bqsr = rules.r2_6_prepare_bqsr.output.bqsr
	output: 
		bam = "results/{run}/bam/{sample}.final.bam"
	log: 
		'results/{run}/logs/prep/{sample}.bqsr_apply.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		wgs_calling_regions = config['references']['wgs_calling_regions']
	threads: workflow.cores/len(ngs.SAMPLES)
	shell: """{params.gatk} ApplyBQSR \
				-R {params.ref} \
				-I {input.bam} \
				-L {params.wgs_calling_regions} \
				-bqsr {input.bqsr} \
				-O {output} &>{log}"""

rule r2_8_bed_to_intervals:
	input: 
		bed = config["panel_capture"]["target"]
	output: 
		intervals = "results/{run}/capture.intervals"
	params:
		picard_old = config['tools']['picard_old'],
		ref_dict = config['references']['dict']
	log: 'results/{run}/logs/prep/intervals.log'
	shell: "java -jar {params.picard_old} BedToIntervalList I={input.bed} SD={params.ref_dict} O={output.intervals} 2> {log}"
