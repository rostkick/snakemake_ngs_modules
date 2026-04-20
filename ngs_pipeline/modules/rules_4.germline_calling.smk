_calling_mode = config.get('calling_mode', 'joint')


rule r4_1_haplotypecaller:
	"""HaplotypeCaller germline variant calling in GVCF mode on BQSR BAM.
	Produces per-sample gVCF for joint genotyping (joint mode) or
	single-sample genotyping (individual mode).
	"""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		bam = rules.r2_9_apply_bqsr.output.bam,
		bai = rules.r2_11_index_bqsr_bam.output,
		bed = rules.r2_2_prepare_bed_to_bed.output.bed
	output:
		gvcf     = temp("results/{run}/germline/vcf/{sample}.gvcf.gz"),
		gvcf_tbi = temp("results/{run}/germline/vcf/{sample}.gvcf.gz.tbi")
	benchmark:
		'results/{run}/benchmarks/germline/vcf/{sample}.haplotypecaller.bm'
	log:
		'results/{run}/logs/germline/{sample}.haplotypecaller.log'
	priority: 35
	params:
		gatk      = config['tools']['gatk'],
		ref       = config['references']['genome_fa'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	threads: {'panel': 4, 'WES': 8, 'WGS': 16}.get(config['ngs_type'], 8)
	resources:
		mem_mb      = {'panel': 8000,  'WES': 16000, 'WGS': 32000}.get(config['ngs_type'], 16000),
		runtime_min = {'panel': 720,   'WES': 8640,  'WGS': 17280}.get(config['ngs_type'], 8640)
	shell: """
		set -euo pipefail
		{params.gatk} --java-options "{params.java_opts}" HaplotypeCaller \
			-R {params.ref} \
			-I {input.bam} \
			-O {output.gvcf} \
			-L {input.bed} \
			-ERC GVCF \
			-G StandardAnnotation \
			-G AS_StandardAnnotation \
			--native-pair-hmm-threads {threads} \
			&>{log}
		{params.gatk} --java-options "{params.java_opts}" IndexFeatureFile \
			-I {output.gvcf} &>>{log}
	"""


if _calling_mode == 'joint':

	rule r4_2_genomicsdbimport:
		"""Import all per-sample gVCFs into a single GenomicsDB workspace.
		extra_gvcfs allows adding archived gVCFs from previous runs (WES only).
		"""
		input:
			gvcf     = expand(rules.r4_1_haplotypecaller.output.gvcf,     run=config['run'], sample=ngs.GRM_SAMPLES),
			gvcf_tbi = expand(rules.r4_1_haplotypecaller.output.gvcf_tbi, run=config['run'], sample=ngs.GRM_SAMPLES)
		output:
			db = directory("results/{run}/germline/genomicsdb/cohort")
		params:
			gatk       = config['tools']['gatk'],
			db_dir     = lambda wc: f"results/{wc.run}/germline/genomicsdb/cohort",
			intervals  = config['references']['wgs_calling_regions'],
			sample_map = lambda wc, input: "\n".join(
				f"{s}\t{g}" for s, g in zip(ngs.GRM_SAMPLES, input.gvcf)
			),
			extra_map  = "\n".join(
				f"{__import__('os').path.basename(g).replace('.gvcf.gz','')}\t{g}"
				for g in sorted(__import__('glob').glob(
					config['references']['wes_joint_gvcfs'].rstrip('/') + '/*.gvcf.gz'
				))
			) if config['references'].get('wes_joint_gvcfs', '') else '',
			java_opts  = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m -Xms4g"
		threads: 8
		resources:
			mem_mb      = {'panel': 8000, 'WES': 16000, 'WGS': 32000}.get(config['ngs_type'], 16000),
			runtime_min = {'panel': 240,  'WES': 720,   'WGS': 2880}.get(config['ngs_type'], 720)
		benchmark:
			'results/{run}/benchmarks/germline/vcf/genomicsdbimport.bm'
		log:
			'results/{run}/logs/germline/genomicsdbimport.log'
		shell: """
			set -euo pipefail
			SAMPLE_MAP=$(mktemp)
			trap "rm -f $SAMPLE_MAP" EXIT
			printf "{params.sample_map}\n{params.extra_map}" | grep -v '^$' > $SAMPLE_MAP
			rm -rf {params.db_dir}
			{params.gatk} --java-options "{params.java_opts}" GenomicsDBImport \
				--sample-name-map  $SAMPLE_MAP \
				--genomicsdb-workspace-path {params.db_dir} \
				-L {params.intervals} \
				--reader-threads {threads} \
				--batch-size 50 \
				&>{log}
		"""


	rule r4_3_genotypegvcfs:
		"""Joint genotyping from GenomicsDB. Produces cohort VCF, left-normalised."""
		input:
			db = rules.r4_2_genomicsdbimport.output.db
		output:
			vcf = "results/{run}/germline/vcf/cohort.vcf.gz",
			tbi = "results/{run}/germline/vcf/cohort.vcf.gz.tbi"
		params:
			gatk      = config['tools']['gatk'],
			bcftools  = config['tools']['bcftools'],
			ref       = config['references']['genome_fa'],
			java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
		threads: 4
		resources:
			mem_mb      = {'panel': 4000, 'WES': 8000, 'WGS': 16000}.get(config['ngs_type'], 8000),
			runtime_min = {'panel': 120,  'WES': 480,  'WGS': 1440}.get(config['ngs_type'], 480)
		benchmark:
			'results/{run}/benchmarks/germline/vcf/genotypegvcfs.bm'
		log:
			'results/{run}/logs/germline/genotypegvcfs.log'
		shell: """
			set -euo pipefail
			{params.gatk} --java-options "{params.java_opts}" GenotypeGVCFs \
				-R {params.ref} \
				-V gendb://{input.db} \
				-O /tmp/cohort.raw.vcf.gz \
				-G StandardAnnotation \
				-G AS_StandardAnnotation \
				&>{log}
			{params.bcftools} norm \
				-m -any \
				-f {params.ref} \
				--output-type z \
				--threads {threads} \
				-o {output.vcf} \
				/tmp/cohort.raw.vcf.gz &>>{log}
			{params.bcftools} index -t {output.vcf} &>>{log}
			rm -f /tmp/cohort.raw.vcf.gz /tmp/cohort.raw.vcf.gz.tbi
		"""


	rule r4_5_vqsr_snp:
		"""VQSR model training for SNPs.
		Requires sufficient variant counts — use joint cohorts of ≥10 WGS or ≥30 WES samples.
		"""
		input:
			vcf = rules.r4_3_genotypegvcfs.output.vcf,
			tbi = rules.r4_3_genotypegvcfs.output.tbi
		output:
			recal    = temp("results/{run}/germline/vcf/cohort.snp.recal"),
			tranches = temp("results/{run}/germline/vcf/cohort.snp.tranches")
		params:
			gatk         = config['tools']['gatk'],
			ref          = config['references']['genome_fa'],
			hapmap       = config['references']['vqsr']['hapmap'],
			omni         = config['references']['vqsr']['omni'],
			thousandg    = config['references']['vqsr']['thousandg_snps'],
			dbsnp        = config['references']['snps'],
			java_opts    = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
		threads: 4
		resources:
			mem_mb      = 16000,
			runtime_min = 480
		benchmark:
			'results/{run}/benchmarks/germline/vcf/vqsr_snp.bm'
		log:
			'results/{run}/logs/germline/vqsr_snp.log'
		shell: """
			set -euo pipefail
			{params.gatk} --java-options "{params.java_opts}" VariantRecalibrator \
				-R {params.ref} \
				-V {input.vcf} \
				--resource:hapmap,known=false,training=true,truth=true,prior=15.0  {params.hapmap} \
				--resource:omni,known=false,training=true,truth=true,prior=12.0    {params.omni} \
				--resource:1000G,known=false,training=true,truth=false,prior=10.0  {params.thousandg} \
				--resource:dbsnp,known=true,training=false,truth=false,prior=7.0   {params.dbsnp} \
				-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
				-mode SNP \
				--dont-run-rscript true \
				--tranches-file {output.tranches} \
				-O {output.recal} \
				&>{log}
		"""


	rule r4_6_apply_vqsr_snp:
		"""Apply SNP VQSR at 99.5% sensitivity tranche."""
		input:
			vcf     = rules.r4_3_genotypegvcfs.output.vcf,
			tbi     = rules.r4_3_genotypegvcfs.output.tbi,
			recal   = rules.r4_5_vqsr_snp.output.recal,
			tranches= rules.r4_5_vqsr_snp.output.tranches
		output:
			vcf = temp("results/{run}/germline/vcf/cohort.snp_vqsr.vcf.gz"),
			tbi = temp("results/{run}/germline/vcf/cohort.snp_vqsr.vcf.gz.tbi")
		params:
			gatk      = config['tools']['gatk'],
			ref       = config['references']['genome_fa'],
			java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
		threads: 2
		resources:
			mem_mb      = 8000,
			runtime_min = 120
		benchmark:
			'results/{run}/benchmarks/germline/vcf/apply_vqsr_snp.bm'
		log:
			'results/{run}/logs/germline/apply_vqsr_snp.log'
		shell: """
			set -euo pipefail
			{params.gatk} --java-options "{params.java_opts}" ApplyVQSR \
				-R {params.ref} \
				-V {input.vcf} \
				--recal-file    {input.recal} \
				--tranches-file {input.tranches} \
				--truth-sensitivity-filter-level 99.5 \
				--create-output-variant-index true \
				-mode SNP \
				-O {output.vcf} \
				&>{log}
		"""


	rule r4_7_vqsr_indel:
		"""VQSR model training for indels."""
		input:
			vcf = rules.r4_6_apply_vqsr_snp.output.vcf,
			tbi = rules.r4_6_apply_vqsr_snp.output.tbi
		output:
			recal    = temp("results/{run}/germline/vcf/cohort.indel.recal"),
			tranches = temp("results/{run}/germline/vcf/cohort.indel.tranches")
		params:
			gatk      = config['tools']['gatk'],
			ref       = config['references']['genome_fa'],
			mills     = config['references']['vqsr']['mills'],
			axiom     = config['references']['vqsr']['axiom_poly'],
			dbsnp     = config['references']['snps'],
			java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
		threads: 4
		resources:
			mem_mb      = 16000,
			runtime_min = 480
		benchmark:
			'results/{run}/benchmarks/germline/vcf/vqsr_indel.bm'
		log:
			'results/{run}/logs/germline/vqsr_indel.log'
		shell: """
			set -euo pipefail
			{params.gatk} --java-options "{params.java_opts}" VariantRecalibrator \
				-R {params.ref} \
				-V {input.vcf} \
				--resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mills} \
				--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 {params.axiom} \
				--resource:dbsnp,known=true,training=false,truth=false,prior=2.0  {params.dbsnp} \
				-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
				-mode INDEL \
				--max-gaussians 4 \
				--dont-run-rscript true \
				--tranches-file {output.tranches} \
				-O {output.recal} \
				&>{log}
		"""


	rule r4_8_apply_vqsr_indel:
		"""Apply indel VQSR at 99.0% sensitivity tranche.
		Output is the final filtered cohort VCF, ready for annotation.
		"""
		input:
			vcf      = rules.r4_6_apply_vqsr_snp.output.vcf,
			tbi      = rules.r4_6_apply_vqsr_snp.output.tbi,
			recal    = rules.r4_7_vqsr_indel.output.recal,
			tranches = rules.r4_7_vqsr_indel.output.tranches
		output:
			vcf = "results/{run}/germline/vcf/cohort.filtered.vcf.gz",
			tbi = "results/{run}/germline/vcf/cohort.filtered.vcf.gz.tbi"
		params:
			gatk      = config['tools']['gatk'],
			ref       = config['references']['genome_fa'],
			java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
		threads: 2
		resources:
			mem_mb      = 8000,
			runtime_min = 120
		benchmark:
			'results/{run}/benchmarks/germline/vcf/apply_vqsr_indel.bm'
		log:
			'results/{run}/logs/germline/apply_vqsr_indel.log'
		shell: """
			set -euo pipefail
			{params.gatk} --java-options "{params.java_opts}" ApplyVQSR \
				-R {params.ref} \
				-V {input.vcf} \
				--recal-file    {input.recal} \
				--tranches-file {input.tranches} \
				--truth-sensitivity-filter-level 99.0 \
				--create-output-variant-index true \
				-mode INDEL \
				-O {output.vcf} \
				&>{log}
		"""


else:  # individual mode

	rule r4_2_genotypegvcfs_individual:
		"""Single-sample genotyping from gVCF. Used in individual calling mode."""
		wildcard_constraints:
			sample = "|".join(ngs.GRM_SAMPLES)
		input:
			gvcf     = rules.r4_1_haplotypecaller.output.gvcf,
			gvcf_tbi = rules.r4_1_haplotypecaller.output.gvcf_tbi
		output:
			vcf = temp("results/{run}/germline/vcf/{sample}.raw.vcf.gz"),
			tbi = temp("results/{run}/germline/vcf/{sample}.raw.vcf.gz.tbi")
		params:
			gatk      = config['tools']['gatk'],
			ref       = config['references']['genome_fa'],
			java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
		threads: 2
		resources:
			mem_mb      = {'panel': 4000, 'WES': 8000, 'WGS': 16000}.get(config['ngs_type'], 8000),
			runtime_min = {'panel': 60,   'WES': 240,  'WGS': 720}.get(config['ngs_type'], 240)
		benchmark:
			'results/{run}/benchmarks/germline/vcf/{sample}.genotypegvcfs.bm'
		log:
			'results/{run}/logs/germline/{sample}.genotypegvcfs.log'
		shell: """
			set -euo pipefail
			{params.gatk} --java-options "{params.java_opts}" GenotypeGVCFs \
				-R {params.ref} \
				-V {input.gvcf} \
				-O {output.vcf} \
				-G StandardAnnotation \
				-G AS_StandardAnnotation \
				&>{log}
		"""


	rule r4_3_hard_filter_individual:
		"""Hard filtering for single-sample calls (no VQSR — insufficient variants).
		SNP and indel filters follow GATK best practices.
		Output VCF name matches what downstream rules expect as {sample}.vcf.gz.
		"""
		wildcard_constraints:
			sample = "|".join(ngs.GRM_SAMPLES)
		input:
			vcf = rules.r4_2_genotypegvcfs_individual.output.vcf,
			tbi = rules.r4_2_genotypegvcfs_individual.output.tbi
		output:
			vcf = "results/{run}/germline/vcf/{sample}.vcf.gz",
			tbi = "results/{run}/germline/vcf/{sample}.vcf.gz.tbi"
		params:
			gatk      = config['tools']['gatk'],
			bcftools  = config['tools']['bcftools'],
			ref       = config['references']['genome_fa'],
			java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
		threads: 2
		resources:
			mem_mb      = 4000,
			runtime_min = 60
		benchmark:
			'results/{run}/benchmarks/germline/vcf/{sample}.hard_filter.bm'
		log:
			'results/{run}/logs/germline/{sample}.hard_filter.log'
		shell: """
			set -euo pipefail

			# Filter SNPs
			{params.gatk} --java-options "{params.java_opts}" SelectVariants \
				-R {params.ref} -V {input.vcf} --select-type-to-include SNP \
				-O /tmp/{wildcards.sample}.snp.vcf.gz &>{log}

			{params.gatk} --java-options "{params.java_opts}" VariantFiltration \
				-R {params.ref} -V /tmp/{wildcards.sample}.snp.vcf.gz \
				--filter-expression "QD < 2.0"          --filter-name "QD2" \
				--filter-expression "FS > 60.0"          --filter-name "FS60" \
				--filter-expression "MQ < 40.0"          --filter-name "MQ40" \
				--filter-expression "MQRankSum < -12.5"  --filter-name "MQRankSum-12.5" \
				--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
				-O /tmp/{wildcards.sample}.snp.filtered.vcf.gz &>>{log}

			# Filter indels
			{params.gatk} --java-options "{params.java_opts}" SelectVariants \
				-R {params.ref} -V {input.vcf} --select-type-to-include INDEL \
				-O /tmp/{wildcards.sample}.indel.vcf.gz &>>{log}

			{params.gatk} --java-options "{params.java_opts}" VariantFiltration \
				-R {params.ref} -V /tmp/{wildcards.sample}.indel.vcf.gz \
				--filter-expression "QD < 2.0"               --filter-name "QD2" \
				--filter-expression "FS > 200.0"              --filter-name "FS200" \
				--filter-expression "ReadPosRankSum < -20.0"  --filter-name "ReadPosRankSum-20" \
				-O /tmp/{wildcards.sample}.indel.filtered.vcf.gz &>>{log}

			# Merge, left-normalise
			{params.bcftools} concat --allow-overlaps \
				/tmp/{wildcards.sample}.snp.filtered.vcf.gz \
				/tmp/{wildcards.sample}.indel.filtered.vcf.gz | \
			{params.bcftools} sort | \
			{params.bcftools} norm \
				-m -any -f {params.ref} \
				--output-type z -o {output.vcf} &>>{log}
			{params.bcftools} index -t {output.vcf} &>>{log}

			rm -f /tmp/{wildcards.sample}.snp.vcf.gz* \
			      /tmp/{wildcards.sample}.snp.filtered.vcf.gz* \
			      /tmp/{wildcards.sample}.indel.vcf.gz* \
			      /tmp/{wildcards.sample}.indel.filtered.vcf.gz*
		"""