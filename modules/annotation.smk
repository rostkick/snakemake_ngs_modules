rule vep_germline:
	wildcard_constraints: sample="|".join(GERMLINE_SAMPLES)
	input: "results/{run}/germline/vcf/{sample}.37.sorted.vcf.gz"
	output: "results/{run}/germline/vcf/{sample}.37.annotated.vcf.gz"
	params: 
		vep=config['tools']['vep']['path'],
		anno_assembly=config['anno_assembly'],
		reference_fasta=config['references']['genome_fa'],
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
				--assembly {params.anno_assembly} \
				--af \
				--af_gnomad \
				--max_af \
				--hgvs \
				--no_escape \
				--canonical \
				--plugin MPC,{params.mpc_data} \
				--compress_output bgzip \
				--use_given_ref \
				--fasta {params.reference_fasta} \
				--dir_cache {params.cache} \
				--dir_plugins {params.plugins} \
				--input_file {input} \
				--output_file {output} \
				--fork {threads} 2>{log}
			"""

use rule vep_germline as vep_germline_joint with:
	input:
		"results/{run}/germline/vcf/cohort.37.sorted.vcf.gz"
	output:
		'results/{run}/germline/vcf/cohort.37.annotated.vcf.gz'
	log:
		'results/{run}/logs/germline/annotation.log'

use rule vep_germline as vep_somatic with:
	input: 
		'results/{run}/somatic/{patient}/mutect2.37.sorted.vcf.gz'
	output: 
		'results/{run}/somatic/{patient}/mutect2.37.annotated.vcf.gz'
	log: 
		'results/{run}/logs/somatic/{patient}/annotation.log'

use rule vep_germline as vep_sv_germline with:
	input: 
		'results/{run}/germline/sv/results/variants/diploidSV.inv_converted.vcf.gz'
	output: 
		'results/{run}/germline/sv/sv.annotated.vcf.gz'
	log: 
		'results/{run}/logs/germline/sv/annotation.log'
