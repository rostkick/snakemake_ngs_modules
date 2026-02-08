rule r9_1_filter_panel_genes:
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input: 
		vcf = rules.r8_2_vep_germline_individual.output.vcf,
		bed = config['panel_capture']['target']
	output:
		vcf = "results/{run}/germline/vcf/{sample}.panel_filtered.vcf.gz"
	params:
		bcftools = config['tools']['bcftools']
	log:
		'results/{run}/logs/parse_results/{sample}.panel_filter.log'
	shell: """
		gene_file=$(mktemp)
		trap "rm -f $gene_file" EXIT
		
		grep -v "^browser" {input.bed} | \
			grep -v "^track" | \
			grep -v "^#" | \
			awk '{{print $4}}' | \
			sort -u > "$gene_file"
		
		{params.bcftools} +split-vep {input.vcf} \
			-c SYMBOL \
			-i "SYMBOL=@$gene_file" \
			-d -x \
			-O z -o {output.vcf} 2>>{log} || exit 1
	"""

rule r9_2_parse_vcf_individual:
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input: 
		vcf = rules.r9_1_filter_panel_genes.output.vcf if config.get('ngs_type') == 'panel' else rules.r8_2_vep_germline_individual.output.vcf
	output:
		tsv = "results/{run}/germline/tsv/{sample}.unfiltered.tsv"
	params:
		bcftools = config['tools']['bcftools']
	log:
		'results/{run}/logs/parse_results/{sample}.parse_vcf.log'
	shell: """
		echo -e 'Chr\\tRef\\tAlt\\tRefGene\\tExon\\tHGVS_description\\tZyg\\t\
rsID\\tGATK_FILTER\\tConsequence\\tConsequence_AA\\t\
gnomAD_exome_NFE\\tgnomAD_exome_Comb\\tgnomAD_genome_NFE\\tgnomAD_genome_Comb\\t\
CADD_PHRED\\tCADD_RAW\\tClinVar_CLNSIG\\tSIFT\\tPolyPhen\\tExAC_pLI\\tPUBMED\\t\
ClinVar_publications\\tClinVar_CLNDN\\tBIOTYPE\\tCANONICAL\\tPHENOTYPES\\t\
am_class\\tam_pathogenicity\\tSNPred_score\\tCoverageDepth\\tGenotypeQual\\tAlleleDepth' > {output.tsv}
		
		{params.bcftools} view -m2 -M2 {input.vcf} | \
		{params.bcftools} +split-vep -s primary -d -f \
'%CHROM:%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%EXON\\t%SYMBOL;%HGVSg;%HGVSc;%HGVSp\\t[%GT]\\t\
%Existing_variation\\t%FILTER\\t%Consequence\\t%Amino_acids\\t\
%gnomADe_NFE_AF\\t%gnomADe_AF\\t%gnomADg_NFE_AF\\t%gnomADg_AF\\t\
%CADD_PHRED\\t%CADD_RAW\\t%CLIN_SIG\\t%SIFT\\t%PolyPhen\\t%pLI_gene_value\\t%PUBMED\\t\
%ClinVar\\t%ClinVar_CLNDN\\t%BIOTYPE\\t%CANONICAL\\t%PHENOTYPES\\t\
%am_class\\t%am_pathogenicity\\t%SNPred_SNPred_score\\t[%DP]\\t[%GQ]\\t[%AD]\\n' >> {output.tsv} 2>>{log}
	"""

rule r9_3_filter_tsv_individual:
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input: 
		tsv = rules.r9_2_parse_vcf_individual.output.tsv,
		bed = config['panel_capture']['target']
	output:
		tsv = "results/{run}/germline/tsv/{sample}.tsv"
	params:
		mart = config['references']['vep_plugins_data']['custom']['mart'],
		rank = config['references']['vep_plugins_data']['custom']['rank'],
		gnomad_filter = config.get('filters', {}).get('gnomad', {}).get('enabled', False),
		gnomad_af_threshold = config.get('filters', {}).get('gnomad', {}).get('af_threshold', 0.01),
		gnomad_column = config.get('filters', {}).get('gnomad', {}).get('column', 'gnomAD_exome_NFE'),
		bed_file = config['panel_capture']['target'],
		ngs_type = config['ngs_type']
	log:
		'results/{run}/logs/parse_results/{sample}.filter_tsv.log'
	script: "scripts/postprocess.py"

rule r9_4_merge_tsv_to_xlsx:
	input:
		tsv = expand("results/{run}/germline/tsv/{sample}.tsv", run=config['run'], sample=ngs.GRM_SAMPLES)
	output:
		xlsx = "results/{run}/germline/xlsx/individual.{run}.germline.results.xlsx"
	params:
		tsv = "results/{run}/run_table.tsv"
	log:
		'results/{run}/logs/parse_results/tsv2xlsx.log'
	script:
		"scripts/collect_tsv_to_xml.py"