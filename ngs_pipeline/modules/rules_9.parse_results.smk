rule r9_1_parse_vcf_individual:
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

rule r9_2_filter_tsv_individual:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_SAMPLES)
	input: 
		tsv = rules.r9_1_parse_vcf_individual.output.tsv
	output:
		tsv = "results/{run}/germline/tsv/{sample}.tsv"
	script: "scripts/filter_ind_table.py"

def get_input_files(wildcards):
	samples = ngs.GRM_SAMPLES
	return [f"data/{sample}.tsv" for sample in samples]

rule r9_3_merge_tsv_to_xlsx:
	input:
		tsv = expand("results/{run}/germline/tsv/{sample}.tsv", run=config['run'], sample=ngs.GRM_SAMPLES)
	output:
		xlsx = "results/{run}/germline/xlsx/individual.germline.results.xlsx"

	script:
		"scripts/collect_tsv_to_xml.py"