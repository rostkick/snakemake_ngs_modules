rule r8_1_vep_germline_joint:
	input:
		vcf = rules.r4_8_mergevcfs.output.vcf
	output:
		vcf = 'results/{run}/germline/vcf/cohort.annotated.vcf.gz'
	params: 
		singularity = config['tools']['singularity'],
		assembly = config['assembly'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		vep = config['tools']['vep']['path'],
		cache = config['tools']['vep']['cache'],
		plugins = config['tools']['vep']['plugins'],
		mpc_data = config['tools']['vep']['plugins_data']['MPC']
	log:
		'results/{run}/logs/germline/annotation.log'
	threads: 
		workflow.cores
	shell: """{params.singularity} run -B /mnt:/mnt {params.vep} \
				/opt/vep/src/ensembl-vep/vep \
				--cache \
				--offline \
				--format vcf \
				--vcf \
				--force_overwrite \
				--force \
				--assembly {params.assembly} \
				--af \
				--af_gnomad \
				--max_af \
				--hgvs \
				--no_escape \
				--canonical \
				--plugin MPC,{params.mpc_data} \
				--compress_output bgzip \
				--use_given_ref \
				--fasta {params.ref} \
				--dir_cache {params.cache} \
				--dir_plugins {params.plugins} \
				--input_file {input.vcf} \
				--output_file {output.vcf} \
				--fork {threads} 2>{log}
			"""

use rule r8_1_vep_germline_joint as r8_2_vep_somatic_paired with:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		vcf = rules.r6_5_filter_pass_exclude_normal_paired.output.vcf
	output:
		vcf = 'results/{run}/somatic/{patient}/somatic_annotated.vcf.gz'
	log: 
		'results/{run}/logs/somatic/{patient}/annotation.log'

use rule r8_1_vep_germline_joint as r8_3_vep_somatic_tmr_only with:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		vcf = rules.r7_3_filter_mutect_calls_tmr_only.output.vcf
	output:
		vcf = 'results/{run}/somatic/{patient}/somatic_annotated_tonly.vcf.gz'
	log: 
		'results/{run}/logs/somatic/{patient}/annotation.log'
