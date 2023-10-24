rule r9_vep_germline_joint:
	input:
		vcf = rules.r4_mergevcfs.output.vcf
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

# use rule r9_vep_germline_joint as r9_vep_germline_sample with:
# 	wildcard_constraints: 
# 		sample = "|".join(ngs.GRM_SAMPLES)
# 	input: 
# 		vcf = rules.r4_deepvariant.output.vcf
# 	output: 
# 		vcf = "results/{run}/germline/vcf/{sample}.annotated.vcf.gz"
# 	log: 
# 		'results/{run}/logs/germline/{sample}.annotation.log'


use rule r9_vep_germline_joint as r9_vep_somatic with:
	input: 
		vcf = rules.r6_filter_pass_exclude_normal_grm_vs_tmr.output.vcf if (ngs.GRM & ngs.TMR) else rules.r7_filter_mutect_calls_tmr_only.output.vcf
	output:
		vcf = 'results/{run}/somatic/{patient}/annotated.vcf.gz'
	log: 
		'results/{run}/logs/somatic/{patient}/annotation.log'

# use rule r9_vep_germline as r9_vep_sv_germline with:
# 	input: 
# 		vcf = 'results/{run}/germline/sv/results/variants/diploidSV.inv_converted.vcf.gz'
# 	output: 
# 		vcf = 'results/{run}/germline/sv/sv.annotated.vcf.gz'
# 	log: 
# 		'results/{run}/logs/germline/sv/annotation.log'
