rule r0_1_archive_dedup_bams:
	"""Archive dedup BAM used as DeepVariant input (no BQSR)."""
	input:
		bam = "results/{run}/bam/{sample}.dedup.bam",
		bai = "results/{run}/bam/{sample}.dedup.bam.bai"
	output:
		bam = "archive/{run}/bam/{sample}.dedup.bam",
		bai = "archive/{run}/bam/{sample}.dedup.bam.bai"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 1000,
		runtime_min = 240,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_dedup_bam.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_dedup_bam.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.bam})
		rsync -av {input.bam} {output.bam} &>>{log}
		rsync -av {input.bai} {output.bai} &>>{log}
	"""


rule r0_1b_archive_bqsr_bams:
	"""Archive BQSR BAM used as Mutect2 input."""
	input:
		bam = "results/{run}/bam/{sample}.final.bqsr.bam",
		bai = "results/{run}/bam/{sample}.final.bqsr.bam.bai"
	output:
		bam = "archive/{run}/bam/{sample}.final.bqsr.bam",
		bai = "archive/{run}/bam/{sample}.final.bqsr.bam.bai"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 1000,
		runtime_min = 240,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_bqsr_bam.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_bqsr_bam.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.bam})
		rsync -av {input.bam} {output.bam} &>>{log}
		rsync -av {input.bai} {output.bai} &>>{log}
	"""


rule r0_2_archive_vcfs:
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		vcf      = "results/{run}/germline/vcf/{sample}.vcf.gz",
		vcf_tbi  = "results/{run}/germline/vcf/{sample}.vcf.gz.tbi",
		gvcf     = "results/{run}/germline/vcf/{sample}.gvcf.gz",
		gvcf_tbi = "results/{run}/germline/vcf/{sample}.gvcf.gz.tbi"
	output:
		vcf      = "archive/{run}/germline/vcf/{sample}.vcf.gz",
		vcf_tbi  = "archive/{run}/germline/vcf/{sample}.vcf.gz.tbi",
		gvcf     = "archive/{run}/germline/vcf/{sample}.gvcf.gz",
		gvcf_tbi = "archive/{run}/germline/vcf/{sample}.gvcf.gz.tbi"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 1000,
		runtime_min = 120,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_vcf.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_vcf.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.vcf})
		rsync -av {input.vcf}      {output.vcf}      &>>{log}
		rsync -av {input.vcf_tbi}  {output.vcf_tbi}  &>>{log}
		rsync -av {input.gvcf}     {output.gvcf}     &>>{log}
		rsync -av {input.gvcf_tbi} {output.gvcf_tbi} &>>{log}
	"""


rule r0_3_archive_tsv:
	"""Archive per-sample TSV for panel and WGS modes."""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		tsv = "results/{run}/germline/tsv/{sample}.tsv"
	output:
		tsv = "archive/{run}/germline/tsv/{sample}.tsv"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 60,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_tsv.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_tsv.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.tsv})
		rsync -av {input.tsv} {output.tsv} &>>{log}
	"""


rule r0_4_archive_tsv_wes_panel:
	"""Archive per-sample per-panel TSV for WES mode."""
	wildcard_constraints:
		sample     = "|".join(ngs.GRM_SAMPLES),
		panel_name = "|".join([p['name'] for p in config['references'].get('wes_gene_panels', [])]) or "NOPANEL"
	input:
		tsv = rules.r9_6_filter_tsv.output.tsv
	output:
		tsv = "archive/{run}/germline/tsv/{sample}.{panel_name}.tsv"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 60,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.{panel_name}.archive_tsv.bm'
	log:
		"results/{run}/logs/archive/{sample}.{panel_name}.archive_tsv.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.tsv})
		rsync -av {input.tsv} {output.tsv} &>>{log}
	"""


rule r0_5_archive_tsv_wes_clinical:
	"""Archive per-sample full WES TSV."""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		tsv = rules.r9_7_filter_tsv_wes_clinical.output.tsv
	output:
		tsv = "archive/{run}/germline/tsv/{sample}.wes_clinical.tsv"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 60,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.wes_clinical.archive_tsv.bm'
	log:
		"results/{run}/logs/archive/{sample}.wes_clinical.archive_tsv.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.tsv})
		rsync -av {input.tsv} {output.tsv} &>>{log}
	"""


rule r0_6_archive_xlsx:
	input:
		xlsx = "results/{run}/germline/xlsx/individual.{run}.germline.results.xlsx"
	output:
		xlsx = "archive/{run}/germline/xlsx/individual.{run}.germline.results.xlsx"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 60,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/archive_xlsx.bm'
	log:
		"results/{run}/logs/archive/xlsx.{run}.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.xlsx})
		rsync -av {input.xlsx} {output.xlsx} &>>{log}
	"""


rule r0_7_fastq_checksums:
	output:
		md5 = "results/{run}/provenance/fastq_checksums.md5"
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 480,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/fastq_checksums.bm'
	log:
		"results/{run}/logs/archive/fastq_checksums.log"
	script: "scripts/fastq_checksums.py"


rule r0_8_archive_hs_metrics:
	input:
		tsv = "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	output:
		tsv = "archive/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 60,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_hs_metrics.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_hs_metrics.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.tsv})
		rsync -av {input.tsv} {output.tsv} &>>{log}
	"""


rule r0_9_archive_tsv_raw:
	"""Archive unfiltered raw TSV files to network storage."""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		tsv = "results/{run}/germline/tsv/{sample}.raw.tsv"
	output:
		tsv = "archive/{run}/germline/tsv/{sample}.raw.tsv"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 60,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_tsv_raw.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_tsv_raw.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.tsv})
		rsync -av {input.tsv} {output.tsv} &>>{log}
	"""


rule r0_10_archive_xlsx_wes_panel:
	"""Archive per-panel WES xlsx to network storage."""
	wildcard_constraints:
		panel_name = "|".join([p['name'] for p in config['references'].get('wes_gene_panels', [])]) or "NOPANEL"
	input:
		xlsx = "results/{run}/germline/xlsx/panel.{panel_name}.{run}.germline.results.xlsx"
	output:
		xlsx = "archive/{run}/germline/xlsx/panel.{panel_name}.{run}.germline.results.xlsx"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 60,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/panel.{panel_name}.archive_xlsx.bm'
	log:
		"results/{run}/logs/archive/panel.{panel_name}.archive_xlsx.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.xlsx})
		rsync -av {input.xlsx} {output.xlsx} &>>{log}
	"""


rule r0_11_archive_xlsx_wes_clinical:
	"""Archive wes_clinical xlsx to network storage."""
	input:
		xlsx = "results/{run}/germline/xlsx/wes_clinical.{run}.germline.results.xlsx"
	output:
		xlsx = "archive/{run}/germline/xlsx/wes_clinical.{run}.germline.results.xlsx"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 60,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/wes_clinical.archive_xlsx.bm'
	log:
		"results/{run}/logs/archive/wes_clinical.{run}.archive_xlsx.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.xlsx})
		rsync -av {input.xlsx} {output.xlsx} &>>{log}
	"""


rule r0_12_archive_gvcfs_joint:
	"""Archive per-sample gVCF files for joint calling mode.
	Essential for future joint calling with additional samples
	without re-running DeepVariant.
	"""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		gvcf     = rules.r4_1_deepvariant.output.gvcf,
		gvcf_tbi = rules.r4_1_deepvariant.output.gvcf_tbi
	output:
		gvcf     = "archive/{run}/germline/gvcf/{sample}.gvcf.gz",
		gvcf_tbi = "archive/{run}/germline/gvcf/{sample}.gvcf.gz.tbi"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 500,
		runtime_min = 120,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_gvcf_joint.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_gvcf_joint.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.gvcf})
		rsync -av {input.gvcf}     {output.gvcf}     &>>{log}
		rsync -av {input.gvcf_tbi} {output.gvcf_tbi} &>>{log}
	"""


rule r0_13_archive_cohort_vcf:
	"""Archive annotated cohort VCF to network storage."""
	input:
		vcf = rules.r8_3_restore_genotypes.output.vcf,
		tbi = rules.r8_3_restore_genotypes.output.tbi
	output:
		vcf = "archive/{run}/germline/vcf/cohort.annotated.vcf.gz",
		tbi = "archive/{run}/germline/vcf/cohort.annotated.vcf.gz.tbi"
	priority: 50
	threads: 1
	resources:
		mem_mb      = 1000,
		runtime_min = 120,
		nfs_io      = 1
	benchmark:
		'results/{run}/benchmarks/archive/cohort.annotated.archive_vcf.bm'
	log:
		"results/{run}/logs/archive/cohort.annotated.archive_vcf.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.vcf})
		rsync -av {input.vcf} {output.vcf} &>>{log}
		rsync -av {input.tbi} {output.tbi} &>>{log}
	"""