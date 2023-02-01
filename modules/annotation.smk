rule vep_germline:
	wildcard_constraints: sample="|".join(GRM_SAMPLES)
	input: "results/{run}/germline/vcf/{sample}.vcf.gz"
	output: "results/{run}/germline/vcf/{sample}.annotated.vcf.gz"
	params: 
		assembly=config['assembly'],
		ref=config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		vep=config['tools']['vep']['path'],
		cache=config['tools']['vep']['cache'],
		plugins=config['tools']['vep']['plugins'],
		mpc_data=config['tools']['vep']['plugins_data']['MPC']
	log: 'results/{run}/logs/germline/{sample}.annotation.log'
	threads: workflow.cores
	shell: """singularity run -B /mnt:/mnt {params.vep} \
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
				--input_file {input} \
				--output_file {output} \
				--fork {threads} 2>{log}
			"""

use rule vep_germline as vep_germline_joint with:
	input:
		"results/{run}/germline/vcf/cohort.vcf.gz"
	output:
		'results/{run}/germline/vcf/cohort.annotated.vcf.gz'
	log:
		'results/{run}/logs/germline/annotation.log'

use rule vep_germline as vep_somatic with:
	input: 
		'results/{run}/somatic/{patient}/final.vcf.gz'
	output:
		'results/{run}/somatic/{patient}/annotated.vcf.gz'
	log: 
		'results/{run}/logs/somatic/{patient}/annotation.log'

use rule vep_germline as vep_sv_germline with:
	input: 
		'results/{run}/germline/sv/results/variants/diploidSV.inv_converted.vcf.gz'
	output: 
		'results/{run}/germline/sv/sv.annotated.vcf.gz'
	log: 
		'results/{run}/logs/germline/sv/annotation.log'
