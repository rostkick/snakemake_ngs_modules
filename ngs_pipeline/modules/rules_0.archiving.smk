rule r0_1_archive_bams:
	"""Archive final BAM files with indices to network storage"""
	input:
		bam = "results/{run}/bam/{sample}.final.bam",
		bai = "results/{run}/bam/{sample}.final.bam.bai"
	output:
		bam = "archive/{run}/bam/{sample}.final.bam",
		bai = "archive/{run}/bam/{sample}.final.bam.bai"
	threads: 1
	resources:
		mem_mb=1000
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_bam.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_bam.log"
	shell: """
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
	threads: 1
	resources:
		mem_mb=1000
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_vcf.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_vcf.log"
	shell: """
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
	threads: 1
	resources:
		mem_mb=500
	benchmark:
		'results/{run}/benchmarks/archive/{sample}.archive_tsv.bm'
	log:
		"results/{run}/logs/archive/{sample}.archive_tsv.log"
	shell: """
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
	threads: 1
	resources:
		mem_mb=500
	benchmark:
		'results/{run}/benchmarks/archive/archive_xlsx.bm'
	log:
		"results/{run}/logs/archive/xlsx.{run}.log"
	shell: """
		mkdir -p $(dirname {output.xlsx})
		rsync -av {input.xlsx} {output.xlsx} &>>{log}
		echo "Archived individual.{wildcards.run}.germline.results.xlsx to network storage" >> {log}
	"""