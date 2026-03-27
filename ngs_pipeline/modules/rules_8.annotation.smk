rule r8_1_strip_genotypes:
	"""Remove sample genotypes before VEP — annotate sites only, then restore genotypes.
	Significantly reduces VEP runtime on large cohorts.
	"""
	input:
		vcf = rules.r4_3_bcf_to_vcf.output.vcf
	output:
		vcf = temp("results/{run}/germline/vcf/cohort.sites.vcf.gz")
	params:
		bcftools = config['tools']['bcftools']
	resources:
		mem_mb      = 2000,
		runtime_min = 30
	benchmark:
		'results/{run}/benchmarks/germline/vcf/strip_genotypes.bm'
	log:
		'results/{run}/logs/germline/strip_genotypes.log'
	shell: """
		set -euo pipefail
		{params.bcftools} view -G {input.vcf} -Oz -o {output.vcf} 2>{log}
		{params.bcftools} index -t {output.vcf} 2>>{log}
	"""


rule r8_2_vep_germline_joint:
	input:
		vcf = rules.r8_1_strip_genotypes.output.vcf
	output:
		vcf = temp("results/{run}/germline/vcf/cohort.sites.annotated.vcf.gz")
	params:
		singularity    = config['tools']['singularity'],
		assembly       = config['assembly'],
		ref            = config['references']['genome_fa'],
		vep            = config['tools']['vep']['path'],
		cache          = config['tools']['vep']['cache'],
		plugins        = config['tools']['vep']['plugins'],
		cadd_data      = config['references']['vep_plugins_data']['CADD'],
		alpha_missense = config['references']['vep_plugins_data']['AlphaMissense'],
		exacpli        = config['references']['vep_plugins_data']['ExACpLI'],
		clinvar        = config['references']['vep_plugins_data']['custom']['ClinVar'],
		snpred         = config['references']['vep_plugins_data']['custom']['SNPred'],
		phenotypes     = config['references']['vep_plugins_data']['Phenotypes'],
		revel          = config['references']['vep_plugins_data'].get('REVEL', ''),
		mpc            = config['references']['vep_plugins_data'].get('MPC', ''),
		spliceai_snv   = config['references']['vep_plugins_data'].get('SpliceAI_SNV', ''),
		spliceai_indel = config['references']['vep_plugins_data'].get('SpliceAI_INDEL', '')
	benchmark:
		'results/{run}/benchmarks/germline/vcf/anno_joint.bm'
	log:
		'results/{run}/logs/germline/annotation.log'
	threads: 4
	resources:
		mem_mb      = 20000,
		runtime_min = 2880
	shell: """{params.singularity} run -B /ngs_pipeline:/ngs_pipeline {params.vep} \
				/opt/vep/src/ensembl-vep/vep \
				--cache \
				--offline \
				--refseq \
				--format vcf \
				--vcf \
				--force_overwrite \
				--force \
				--assembly {params.assembly} \
				--af_gnomade \
				--af_gnomadg \
				--clin_sig_allele 1 \
				--canonical \
				--hgvs \
				--hgvsg \
				--no_escape \
				--protein \
				--sift b \
				--polyphen b \
				--humdiv \
				--pubmed \
				--domains \
				--plugin CADD,{params.cadd_data} \
				--plugin AlphaMissense,file={params.alpha_missense} \
				--plugin pLI,{params.exacpli} \
				--plugin Phenotypes,file={params.phenotypes} \
				--plugin NMD \
				$([ -n "{params.revel}" ] && echo "--plugin REVEL,file={params.revel}") \
				$([ -n "{params.mpc}" ] && echo "--plugin MPC,{params.mpc}") \
				$([ -n "{params.spliceai_snv}" ] && [ -n "{params.spliceai_indel}" ] && echo "--plugin SpliceAI,snv={params.spliceai_snv},indel={params.spliceai_indel}") \
				--custom {params.clinvar},ClinVar,vcf,exact,0,CLINSIG,CLNDN,CLNREVSTAT \
				--custom {params.snpred},SNPred,vcf,exact,0,SNPred_score \
				--compress_output bgzip \
				--use_given_ref \
				--fasta {params.ref} \
				--dir_cache {params.cache} \
				--dir_plugins {params.plugins} \
				--input_file {input.vcf} \
				--output_file {output.vcf} \
				--fork {threads} 2>{log}
			"""


rule r8_3_restore_genotypes:
	"""Restore sample genotypes by merging annotated sites back into the cohort VCF."""
	input:
		annotated = rules.r8_2_vep_germline_joint.output.vcf,
		cohort    = rules.r4_3_bcf_to_vcf.output.vcf
	output:
		vcf = "results/{run}/germline/vcf/cohort.annotated.vcf.gz",
		tbi = "results/{run}/germline/vcf/cohort.annotated.vcf.gz.tbi"
	params:
		bcftools = config['tools']['bcftools']
	resources:
		mem_mb      = 4000,
		runtime_min = 60
	benchmark:
		'results/{run}/benchmarks/germline/vcf/restore_genotypes.bm'
	log:
		'results/{run}/logs/germline/restore_genotypes.log'
	shell: """
		set -euo pipefail
		{params.bcftools} index -ft {input.annotated} 2>{log}
		{params.bcftools} index -ft {input.cohort} 2>>{log}
		{params.bcftools} annotate \
			-a {input.annotated} \
			-c INFO \
			{input.cohort} \
			-Oz -o {output.vcf} 2>>{log}
		{params.bcftools} index -t {output.vcf} 2>>{log}
	"""


use rule r8_2_vep_germline_joint as r8_4_vep_germline_individual with:
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		vcf = rules.r4_1_deepvariant.output.vcf
	output:
		vcf = temp("results/{run}/germline/vcf/{sample}.annotated.vcf.gz")
	resources:
		mem_mb      = 16000,
		runtime_min = 2880
	priority: 40
	benchmark:
		'results/{run}/benchmarks/germline/vcf/{sample}.anno_individual.bm'
	log:
		'results/{run}/logs/germline/{sample}.annotation.log'


use rule r8_2_vep_germline_joint as r8_5_vep_somatic_paired with:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input:
		vcf = rules.r6_5_filter_pass_exclude_normal_paired.output.vcf
	output:
		vcf = 'results/{run}/somatic/{patient}/somatic_annotated.vcf.gz'
	resources:
		mem_mb      = 16000,
		runtime_min = 2880
	priority: 40
	benchmark:
		'results/{run}/benchmarks/somatic/{patient}/anno_tmr_paired.bm'
	log:
		'results/{run}/logs/somatic/{patient}/annotation.log'


use rule r8_2_vep_germline_joint as r8_6_vep_somatic_tmr_only with:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input:
		vcf = rules.r7_3_filter_mutect_calls_tmr_only.output.vcf
	output:
		vcf = 'results/{run}/somatic/{patient}/somatic_annotated_tonly.vcf.gz'
	resources:
		mem_mb      = 16000,
		runtime_min = 2880
	priority: 40
	benchmark:
		'results/{run}/benchmarks/somatic/{patient}/anno_only_tmr.bm'
	log:
		'results/{run}/logs/somatic/{patient}/annotation.log'