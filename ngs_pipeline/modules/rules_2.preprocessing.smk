rule r2_1_bed_to_intervals:
	input: 
		bed = config["panel_capture"]["target"]
	output: 
		intervals = "results/{run}/capture.intervals"
	params:
		picard = config['tools']['picard'],
		ref_dict = config['references']['dict']
	benchmark:
		'results/{run}/benchmarks/bam/bed2intervals.bm'
	log: 
		'results/{run}/logs/prep/intervals.log'
	resources:
		mem_mb=2000,
		runtime_min=60
	shell: 
		"_JAVA_OPTIONS='-Xmx1600m' {params.picard} BedToIntervalList I={input.bed} SD={params.ref_dict} O={output.intervals} 2> {log}"

rule r2_2_prepare_bed_to_bed:
	input: 
		bed = config["panel_capture"]["target"]
	output: 
		bed = "results/{run}/capture.clean.bed"
	log: 
		'results/{run}/logs/prep/prepare_bed.log'
	shell: """
		grep -v "^browser" {input.bed} | grep -vP "^track|^browser" | grep -v "^#" | cut -f1-3 > {output.bed} 2>{log}
		"""

rule r2_3_sam_to_bam:
	input: 
		sam = rules.r1_1_read_alignment.output.sam
	output:
		bam = temp('results/{run}/bam/{sample}.{lane}.for_sort1.bam')
	params:
		samtools = config['tools']['samtools'],
		release_lock_script = "modules/scripts/release_staged_lock.py"
	threads: 4
	resources:
		mem_mb=2000,
		runtime_min={'panel': 2880, 'WES': 8640, 'WGS': 17280}.get(config['ngs_type'], 8640)
	benchmark:
		'results/{run}/benchmarks/bam/{sample}.{lane}.sam2bam.bm'
	log:
		'results/{run}/logs/prep/{sample}.{lane}.sam2bam.log'
	priority: 30
	shell: """
		set -euo pipefail
		{params.samtools} view -@ {threads} -bS -o {output.bam} {input.sam} 2>{log}
		python {params.release_lock_script} --lock-dir .staging_locks --key {wildcards.sample}.{wildcards.lane} --log {log}
	"""

rule r2_4_sort_premerged_bams:
	input: 
		bam = rules.r2_3_sam_to_bam.output.bam
	output: 
		bam = temp('results/{run}/bam/{sample}.{lane}.for_merge.bam')
	params:
		samtools = config['tools']['samtools']
	threads: {'panel': 2, 'WES': 4, 'WGS': 8}.get(config['ngs_type'], 4)
	resources:
		mem_mb={'panel': 6000, 'WES': 12000, 'WGS': 24000}.get(config['ngs_type'], 12000),
		runtime_min={'panel': 240, 'WES': 720, 'WGS': 1440}.get(config['ngs_type'], 720)
	benchmark:
		'results/{run}/benchmarks/bam/{sample}.{lane}.sort1.bm'
	log:
		'results/{run}/logs/prep/{sample}.{lane}.sort1.log'
	priority: 30
	shell: "{params.samtools} sort -@ {threads} -m 2G -o {output.bam} {input.bam} 2>{log}"

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

rule r2_5_merge_bams:
	input:
		bams = get_merged_bams
	output: 
		bam = temp('results/{run}/bam/{sample}.for_sort2.bam')
	params:
		samtools = config['tools']['samtools']
	threads: {'panel': 2, 'WES': 4, 'WGS': 4}.get(config['ngs_type'], 4)
	resources:
		mem_mb={'panel': 2000, 'WES': 4000, 'WGS': 8000}.get(config['ngs_type'], 4000),
		runtime_min={'panel': 120, 'WES': 480, 'WGS': 960}.get(config['ngs_type'], 480)
	benchmark:
		'results/{run}/benchmarks/bam/{sample}.merge_bam.bm'
	log:
		'results/{run}/logs/prep/{sample}.merge_bam.log'
	priority: 30
	shell: """
		set -euo pipefail
		if [ $(echo {input.bams} | wc -w) -gt 1 ]; then
			{params.samtools} merge -@ {threads} {output.bam} {input.bams} 2>{log}
		else
			cp {input.bams} {output.bam} 2>{log}
		fi
	"""

rule r2_6_sorting_bam:
	input: 
		bam = rules.r2_5_merge_bams.output.bam
	output: 
		bam = temp("results/{run}/bam/{sample}.for_dedup.bam")
	threads: {'panel': 2, 'WES': 4, 'WGS': 8}.get(config['ngs_type'], 4)
	resources:
		mem_mb={'panel': 6000, 'WES': 12000, 'WGS': 24000}.get(config['ngs_type'], 12000),
		runtime_min={'panel': 240, 'WES': 720, 'WGS': 1440}.get(config['ngs_type'], 720)
	params: 
		samtools = config['tools']['samtools']
	benchmark:
		'results/{run}/benchmarks/bam/{sample}.sort2.bm'
	log:
		'results/{run}/logs/prep/{sample}.sort2.log'
	priority: 30
	shell: "{params.samtools} sort -@ {threads} -m 2G -o {output.bam} {input.bam} 2>{log}"

rule r2_7_mark_duplicates:
	input: 
		bam = rules.r2_6_sorting_bam.output.bam
	output: 
		bam = temp("results/{run}/bam/{sample}.for_bqsr.bam")
	params:
		gatk = config['tools']['gatk'],
		samtools = config['tools']['samtools'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	resources:
		mem_mb={'panel': 4000, 'WES': 8000, 'WGS': 20000}.get(config['ngs_type'], 8000),
		runtime_min={'panel': 480, 'WES': 1440, 'WGS': 2880}.get(config['ngs_type'], 1440)
	benchmark:
		'results/{run}/benchmarks/bam/{sample}.mark_dedup.bm'
	log: 
		log1 = "results/{run}/logs/prep/{sample}.dedup.log1",
		log2 = "results/{run}/logs/prep/{sample}.dedup.log2"
	priority: 30
	shell: """
		set -euo pipefail
		if [ "{config[hs]}" = "True" ] || [ "{config[hs]}" = "true" ]; then
			{params.gatk} --java-options "{params.java_opts}" MarkDuplicates \
				-I {input.bam} -O {output.bam} -M {log.log1} 2>{log.log2} \
				--ASSUME_SORTED true
			{params.samtools} index {output.bam}
		else
			cp {input.bam} {output.bam}
			{params.samtools} index {output.bam}
		fi
	"""

rule r2_8_prepare_bqsr:
	input: 
		bam = rules.r2_7_mark_duplicates.output.bam
	output: 
		bqsr = temp("results/{run}/bam/{sample}.for_bqsr.recal.table")
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		snps = config['references']['snps'],
		indels = config['references']['indels'],
		wgs_calling_regions = config['references']['wgs_calling_regions'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	resources:
		mem_mb={'panel': 4000, 'WES': 8000, 'WGS': 10000}.get(config['ngs_type'], 8000),
		runtime_min={'panel': 480, 'WES': 1440, 'WGS': 2880}.get(config['ngs_type'], 1440)
	benchmark:
		'results/{run}/benchmarks/bam/{sample}.bqsr_recal.bm'
	log: 
		'results/{run}/logs/prep/{sample}.bqsr_recal.log'
	priority: 30
	shell: """{params.gatk} --java-options "{params.java_opts}" BaseRecalibrator \
				-R {params.ref} \
				-I {input.bam} \
				--known-sites {params.snps} \
				--known-sites {params.indels} \
				-L {params.wgs_calling_regions} \
				-O {output} &>{log}"""

rule r2_9_apply_bqsr:
	input: 
		bam = rules.r2_7_mark_duplicates.output.bam,
		bqsr = rules.r2_8_prepare_bqsr.output.bqsr
	output: 
		bam = temp("results/{run}/bam/{sample}.final.bam"),
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		panel_capture = config['panel_capture']['target'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	resources:
		mem_mb={'panel': 4000, 'WES': 8000, 'WGS': 10000}.get(config['ngs_type'], 8000),
		runtime_min={'panel': 240, 'WES': 720, 'WGS': 1440}.get(config['ngs_type'], 720)
	benchmark:
		'results/{run}/benchmarks/bam/{sample}.apply_bqsr.bm'
	log: 
		'results/{run}/logs/prep/{sample}.bqsr_apply.log'
	priority: 30
	shell: """{params.gatk} --java-options "{params.java_opts}" ApplyBQSR \
				-R {params.ref} \
				-I {input.bam} \
				-L {params.panel_capture} \
				-bqsr {input.bqsr} \
				-O {output.bam} &>{log}"""

rule r2_10_index_bam:
	input: rules.r2_9_apply_bqsr.output.bam
	output: "results/{run}/bam/{sample}.final.bam.bai"
	params:
		samtools = config['tools']['samtools']
	resources:
		mem_mb=1000,
		runtime_min=60
	priority: 30
	shell: "{params.samtools} index {input}"