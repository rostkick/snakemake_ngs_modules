_glnexus_config = 'DeepVariantWGS' if config['ngs_type'] == 'WGS' else 'DeepVariantWES'


rule r4_1_deepvariant:
	"""DeepVariant germline variant calling on dedup BAM without BQSR.
	DeepVariant is trained on raw base quality scores — BQSR degrades accuracy
	by distorting the quality score patterns the model was trained on.
	"""
	input:
		bam = rules.r2_7_mark_duplicates.output.bam,
		bai = rules.r2_10_index_dedup_bam.output,
		bed = rules.r2_2_prepare_bed_to_bed.output.bed
	output:
		vcf      = temp("results/{run}/germline/vcf/{sample}.vcf.gz"),
		vcf_tbi  = temp("results/{run}/germline/vcf/{sample}.vcf.gz.tbi"),
		gvcf     = temp("results/{run}/germline/vcf/{sample}.gvcf.gz"),
		gvcf_tbi = temp("results/{run}/germline/vcf/{sample}.gvcf.gz.tbi")
	benchmark:
		'results/{run}/benchmarks/germline/vcf/{sample}.deepvariant.bm'
	log:
		'results/{run}/logs/germline/{sample}.deepvariant.log'
	priority: 35
	params:
		singularity = config['tools']['singularity'],
		deepvariant = config['tools']['deepvariant'],
		model_type  = lambda wc: 'WES' if config['ngs_type'] in ['WES', 'panel'] else config['ngs_type'],
		ref         = config['references']['genome_fa']
	threads: {'panel': 4, 'WES': 8, 'WGS': 16}.get(config['ngs_type'], 8)
	resources:
		mem_mb      = {'panel': 8000,  'WES': 16000, 'WGS': 32000}.get(config['ngs_type'], 16000),
		runtime_min = {'panel': 720,   'WES': 8640,  'WGS': 17280}.get(config['ngs_type'], 8640)
	shell: """
		{params.singularity} run \
			-B /ngs_pipeline:/ngs_pipeline \
			{params.deepvariant} /opt/deepvariant/bin/run_deepvariant \
			--model_type={params.model_type} \
			--output_vcf={output.vcf} \
			--output_gvcf={output.gvcf} \
			--reads={input.bam} \
			--ref={params.ref} \
			--regions {input.bed} \
			--num_shards={threads} &>{log}
	"""


rule r4_2_glnexus_joint:
	input:
		gvcf     = expand(rules.r4_1_deepvariant.output.gvcf,     run=config['run'], sample=ngs.GRM_SAMPLES),
		gvcf_tbi = expand(rules.r4_1_deepvariant.output.gvcf_tbi, run=config['run'], sample=ngs.GRM_SAMPLES)
	output:
		bcf = temp("results/{run}/germline/vcf/cohort.glnexus.bcf")
	params:
		singularity    = config['tools']['singularity'],
		glnexus_sif    = config['tools']['glnexus'],
		glnexus_config = _glnexus_config,
		db_dir         = lambda wc: f"/tmp/glnexus_db_{wc.run}",
		extra_gvcfs    = " ".join(
			sorted(__import__('glob').glob(
				config['references']['wes_joint_gvcfs'].rstrip('/') + '/*.gvcf.gz'
			))
		) if config['references'].get('wes_joint_gvcfs', '') else ''
	threads: 8
	resources:
		mem_mb      = {'panel': 8000,  'WES': 16000, 'WGS': 32000}.get(config['ngs_type'], 16000),
		runtime_min = {'panel': 240,   'WES': 720,   'WGS': 2880}.get(config['ngs_type'], 720)
	benchmark:
		'results/{run}/benchmarks/germline/vcf/glnexus_joint.bm'
	log:
		'results/{run}/logs/germline/glnexus_joint.log'
	shell: """
		set -euo pipefail
		rm -rf {params.db_dir}
		{params.singularity} exec \
			-B /ngs_pipeline:/ngs_pipeline \
			{params.glnexus_sif} \
			glnexus_cli \
			--config  {params.glnexus_config} \
			--dir     {params.db_dir} \
			--threads {threads} \
			{params.extra_gvcfs} \
			{input.gvcf} \
			> {output.bcf} \
			2> {log}
		rm -rf {params.db_dir}
	"""


rule r4_3_bcf_to_vcf:
	input:
		bcf = rules.r4_2_glnexus_joint.output.bcf
	output:
		vcf = temp("results/{run}/germline/vcf/cohort.vcf.gz"),
		tbi = temp("results/{run}/germline/vcf/cohort.vcf.gz.tbi")
	params:
		bcftools = config['tools']['bcftools'],
		ref      = config['references']['genome_fa']
	threads: 4
	resources:
		mem_mb      = 4000,
		runtime_min = 120
	benchmark:
		'results/{run}/benchmarks/germline/vcf/bcf_to_vcf.bm'
	log:
		'results/{run}/logs/germline/bcf_to_vcf.log'
	shell: """
		set -euo pipefail
		{params.bcftools} view {input.bcf} 2>>{log} | \
		{params.bcftools} norm \
			-m -any \
			-f {params.ref} \
			--output-type z \
			--threads {threads} \
			-o {output.vcf} 2>>{log}
		{params.bcftools} index -t {output.vcf} 2>>{log}
	"""