rule r0_1_archive_bams:
	"""Archive final BAM files with indices to network storage"""
	input:
		bam = "results/{run}/bam/{sample}.final.bam",
		bai = "results/{run}/bam/{sample}.final.bam.bai"
	output:
		bam = "archive/{run}/bam/{sample}.final.bam",
		bai = "archive/{run}/bam/{sample}.final.bam.bai"
	priority: 50
	threads: 1
	resources:
		mem_mb=1000,
		runtime_min=240,
		nfs_io=1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_bam.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_bam.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.bam})
		rsync -av {input.bam} {output.bam} &>>{log}
		rsync -av {input.bai} {output.bai} &>>{log}
		echo "Archived {wildcards.sample} (BAM + BAI) to network storage" >> {log}
	"""

rule r0_2_archive_vcfs:
	"""Archive germline VCF files with indices to network storage"""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		vcf = "results/{run}/germline/vcf/{sample}.vcf.gz",
		vcf_tbi = "results/{run}/germline/vcf/{sample}.vcf.gz.tbi",
		gvcf = "results/{run}/germline/vcf/{sample}.gvcf.gz",
		gvcf_tbi = "results/{run}/germline/vcf/{sample}.gvcf.gz.tbi"
	output:
		vcf = "archive/{run}/germline/vcf/{sample}.vcf.gz",
		vcf_tbi = "archive/{run}/germline/vcf/{sample}.vcf.gz.tbi",
		gvcf = "archive/{run}/germline/vcf/{sample}.gvcf.gz",
		gvcf_tbi = "archive/{run}/germline/vcf/{sample}.gvcf.gz.tbi"
	priority: 50
	threads: 1
	resources:
		mem_mb=1000,
		runtime_min=120,
		nfs_io=1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_vcf.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_vcf.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.vcf})
		rsync -av {input.vcf} {output.vcf} &>>{log}
		rsync -av {input.vcf_tbi} {output.vcf_tbi} &>>{log}
		rsync -av {input.gvcf} {output.gvcf} &>>{log}
		rsync -av {input.gvcf_tbi} {output.gvcf_tbi} &>>{log}
		echo "Archived {wildcards.sample} (VCF + GVCF with indices) to network storage" >> {log}
	"""

rule r0_3_archive_tsv_individual:
	"""Archive individual TSV results to network storage"""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		tsv = "results/{run}/germline/tsv/{sample}.tsv"
	output:
		tsv = "archive/{run}/germline/tsv/{sample}.tsv"
	priority: 50
	threads: 1
	resources:
		mem_mb=500,
		runtime_min=60,
		nfs_io=1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_tsv.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_tsv.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.tsv})
		rsync -av {input.tsv} {output.tsv} &>>{log}
		echo "Archived {wildcards.sample} TSV to network storage" >> {log}
	"""

rule r0_4_archive_xlsx:
	"""Archive consolidated XLSX results to network storage"""
	input:
		xlsx = "results/{run}/germline/xlsx/individual.{run}.germline.results.xlsx"
	output:
		xlsx = "archive/{run}/germline/xlsx/individual.{run}.germline.results.xlsx"
	priority: 50
	threads: 1
	resources:
		mem_mb=500,
		runtime_min=60,
		nfs_io=1
	benchmark:
		'results/{run}/benchmarks/archive/archive_xlsx.bm'
	log:
		"results/{run}/logs/archive/xlsx.{run}.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.xlsx})
		rsync -av {input.xlsx} {output.xlsx} &>>{log}
		echo "Archived individual.{wildcards.run}.germline.results.xlsx to network storage" >> {log}
	"""

_gatk_gvcf_suffix = '.wgs.gatk.gvcf.gz' if config['ngs_type'] == 'WGS' else '.gatk.gvcf.gz'

rule r0_5_archive_gatk_gvcfs:
	"""Archive per-sample GATK HaplotypeCaller GVCF files to network storage.
	These are valuable for future joint calling with additional samples
	without re-running HaplotypeCaller.
	"""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		gvcf = "results/{run}/germline/vcf/{sample}" + _gatk_gvcf_suffix,
		gvcf_tbi = "results/{run}/germline/vcf/{sample}" + _gatk_gvcf_suffix + ".tbi"
	output:
		gvcf = "archive/{run}/germline/gvcf/{sample}" + _gatk_gvcf_suffix,
		gvcf_tbi = "archive/{run}/germline/gvcf/{sample}" + _gatk_gvcf_suffix + ".tbi"
	priority: 50
	threads: 1
	resources:
		mem_mb=500,
		runtime_min=120,
		nfs_io=1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_gatk_gvcf.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_gatk_gvcf.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.gvcf})
		rsync -av {input.gvcf} {output.gvcf} &>>{log}
		rsync -av {input.gvcf_tbi} {output.gvcf_tbi} &>>{log}
		echo "Archived {wildcards.sample} GATK GVCF to network storage" >> {log}
	"""

rule r0_7_archive_hs_metrics:
	"""Archive HS metrics to network storage"""
	input:
		tsv = "results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	output:
		tsv = "archive/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
	priority: 50
	threads: 1
	resources:
		mem_mb=500,
		runtime_min=60,
		nfs_io=1
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_hs_metrics.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_hs_metrics.log"
	shell: """
		set -euo pipefail
		mkdir -p $(dirname {output.tsv})
		rsync -av {input.tsv} {output.tsv} &>>{log}
	"""

rule r0_6_fastq_checksums:
	"""Compute MD5 checksums for all input FASTQ files in this run.
	Ensures full traceability of input data.
	"""
	output:
		md5 = "results/{run}/provenance/fastq_checksums.md5"
	threads: 1
	resources:
		mem_mb=500,
		runtime_min=480,
		nfs_io=1
	benchmark:
		'results/{run}/benchmarks/archive/fastq_checksums.bm'
	log:
		"results/{run}/logs/archive/fastq_checksums.log"
	script: "scripts/fastq_checksums.py"