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

rule r8_5_parse_vcf_individual:
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input: 
		vcf = rules.r8_2_vep_germline_individual.output.vcf
	output:
		tsv = temp("results/{run}/germline/tsv/{sample}.unfiltered.tsv")
	params:
		bcftools = config['tools']['bcftools'],
	shell: """echo -e \
'Chr\\tRef\\tAlt\\tRefGene\\tExon\\tHGVS_description\\tZyg\\t\
rsID\\tGATK_FILTER\\tConsequence\\tConsequence_AA\\t\
gnomAD_exome_NFE\\tgnomAD_exome_Comb\\tgnomAD_genome_NFE\\tgnomAD_genome_Comb\\t\
CADD_PHRED\\tCADD_RAW\\tSIFT\\tPolyPhen\\tBIOTYPE\\tCANONICAL\\t\
am_class\\tam_pathogenicity\\tCoverageDepth\\tGenotypeQual\\tAlleleDepth' > {output.tsv}; \
{params.bcftools} +split-vep -s primary -d -f \
'%CHROM:%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%EXON\\t%SYMBOL;%HGVSg;%HGVSc;%HGVSp\\t[%GT]\\t\
%Existing_variation\\t%FILTER\\t%Consequence\\t%Amino_acids\\t\
%gnomADe_NFE_AF\\t%gnomADe_AF\\t%gnomADg_NFE_AF\\t%gnomADg_AF\\t\
%CADD_PHRED\\t%CADD_RAW\\t%SIFT\\t%PolyPhen\\t%BIOTYPE\\t%CANONICAL\\t\
%am_class\\t%am_pathogenicity\\t[%DP]\\t[%GQ]\\t[%AD]\\n' {input.vcf} >> {output.tsv}"""

rule r8_6_filter_tsv_individual:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_SAMPLES)
	input: 
		tsv = rules.r8_5_parse_vcf_individual.output.tsv
	output:
		tsv = "results/{run}/germline/tsv/{sample}.tsv"
	script: "scripts/filter_ind_table.py"

def get_input_files(wildcards):
    samples = ngs.GRM_SAMPLES
    return [f"data/{sample}.tsv" for sample in samples]

rule r8_6_merge_tsv_to_xlsx:
	# wildcard_constraints:
	# 	sample = "|".join(ngs.GRM_SAMPLES)
	input:
		tsv = expand("results/{run}/germline/tsv/{sample}.tsv", run=config['run'], sample=ngs.GRM_SAMPLES)
	output:
		xlsx = "results/{run}/germline/xlsx/individual.germline.results.xlsx"

	script:
		"scripts/collect_tsv_to_xml.py"
