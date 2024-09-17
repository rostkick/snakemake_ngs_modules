rule r8_1_vep_germline_joint:
	input:
		vcf = rules.r4_8_mergevcfs.output.vcf
	output:
		vcf = 'results/{run}/germline/vcf/cohort.annotated.vcf.gz'
	params: 
		singularity = config['tools']['singularity'],
		assembly = config['assembly'],
		ref = config['references']['genome_fa'],
		vep = config['tools']['vep']['path'],
		cache = config['tools']['vep']['cache'],
		plugins = config['tools']['vep']['plugins'],
		cadd_data = config['references']['vep_plugins_data']['CADD'],
		alpha_missense = config['references']['vep_plugins_data']['AlphaMissense']
	log:
		'results/{run}/logs/germline/annotation.log'
	threads:
		workflow.cores
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
				--canonical \
				--hgvs \
				--hgvsg \
				--no_escape \
				--protein \
				--sift b \
				--polyphen b \
				--humdiv \
				--domains \
				--plugin CADD,{params.cadd_data} \
				--plugin AlphaMissense,file={params.alpha_missense} \
				--compress_output bgzip \
				--use_given_ref \
				--fasta {params.ref} \
				--dir_cache {params.cache} \
				--dir_plugins {params.plugins} \
				--input_file {input.vcf} \
				--output_file {output.vcf} \
				--fork {threads} 2>{log}
			"""

use rule r8_1_vep_germline_joint as r8_2_vep_germline_individual with:
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input: 
		vcf = rules.r4_9_deepvariant.output.vcf
	output:
		vcf = "results/{run}/germline/vcf/{sample}.annotated.vcf.gz"
	log:
		'results/{run}/logs/germline/{sample}.annotation.log'

use rule r8_1_vep_germline_joint as r8_3_vep_somatic_paired with:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		vcf = rules.r6_5_filter_pass_exclude_normal_paired.output.vcf
	output:
		vcf = 'results/{run}/somatic/{patient}/somatic_annotated.vcf.gz'
	log:
		'results/{run}/logs/somatic/{patient}/annotation.log'

use rule r8_1_vep_germline_joint as r8_4_vep_somatic_tmr_only with:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		vcf = rules.r7_3_filter_mutect_calls_tmr_only.output.vcf
	output:
		vcf = 'results/{run}/somatic/{patient}/somatic_annotated_tonly.vcf.gz'
	log:
		'results/{run}/logs/somatic/{patient}/annotation.log'
