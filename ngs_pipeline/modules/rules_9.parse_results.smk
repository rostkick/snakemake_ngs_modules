# Build panel lookup from config for WES: {name -> bed}
# Empty dict for non-WES ngs_type.
_wes_panels = {
    p['name']: p['bed']
    for p in config['references'].get('wes_gene_panels', [])
} if config['ngs_type'] == 'WES' else {}

_panel_names = list(_wes_panels.keys())


rule r9_1_split_joint_vcf:
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		# vcf = rules.r8_3_restore_genotypes.output.vcf
		vcf = "results/{run}/germline/vcf/cohort.annotated.vcf.gz"
	output:
		vcf     = "results/{run}/germline/vcf/{sample}.joint.annotated.vcf.gz",
		vcf_tbi = "results/{run}/germline/vcf/{sample}.joint.annotated.vcf.gz.tbi"
	params:
		bcftools = config['tools']['bcftools']
	resources:
		mem_mb      = 2000,
		runtime_min = 60
	benchmark:
		'results/{run}/benchmarks/germline/vcf/{sample}.split_joint.bm'
	log:
		'results/{run}/logs/germline/{sample}.split_joint.log'
	shell: """
		set -euo pipefail
		{params.bcftools} view \
			-s {wildcards.sample} \
			-Oz -o {output.vcf} \
			{input.vcf} 2>{log}
		{params.bcftools} index -t {output.vcf} 2>>{log}
	"""


rule r9_2_filter_panel_genes:
	"""Filter per-sample VCF to genes in a WES gene panel BED."""
	wildcard_constraints:
		sample     = "|".join(ngs.GRM_SAMPLES),
		panel_name = "|".join(_panel_names) if _panel_names else "NOPANEL"
	input:
		vcf = rules.r9_1_split_joint_vcf.output.vcf \
			if config.get('calling_mode', 'joint') == 'joint' \
			else rules.r8_4_vep_germline_individual.output.vcf
	output:
		vcf = temp("results/{run}/germline/vcf/{sample}.{panel_name}.panel_filtered.vcf.gz")
	params:
		bcftools = config['tools']['bcftools'],
		bed      = lambda wc: _wes_panels[wc.panel_name]
	priority: 40
	resources:
		mem_mb      = 2000,
		runtime_min = 60
	log:
		'results/{run}/logs/parse_results/{sample}.{panel_name}.panel_filter.log'
	shell: """
		gene_file=$(mktemp)
		trap "rm -f $$gene_file" EXIT

		grep -v "^browser" {params.bed} | \
			grep -v "^track" | \
			grep -v "^#" | \
			awk '{{print $4}}' | \
			sort -u > "$$gene_file"

		{params.bcftools} +split-vep {input.vcf} \
			-c SYMBOL \
			-i "SYMBOL=@$$gene_file" \
			-d -x \
			-O z -o {output.vcf} 2>>{log} || exit 1
	"""


rule r9_3_filter_capture_genes:
	"""Filter per-sample VCF to genes in the capture BED. Used for panel ngs_type."""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		vcf = rules.r9_1_split_joint_vcf.output.vcf \
			if config.get('calling_mode', 'joint') == 'joint' \
			else rules.r8_4_vep_germline_individual.output.vcf,
		bed = config['panel_capture']['target']
	output:
		vcf = temp("results/{run}/germline/vcf/{sample}.capture_filtered.vcf.gz")
	params:
		bcftools = config['tools']['bcftools']
	priority: 40
	resources:
		mem_mb      = 2000,
		runtime_min = 60
	log:
		'results/{run}/logs/parse_results/{sample}.capture_filter.log'
	shell: """
		gene_file=$(mktemp)
		trap "rm -f $$gene_file" EXIT

		grep -v "^browser" {input.bed} | \
			grep -v "^track" | \
			grep -v "^#" | \
			awk '{{print $4}}' | \
			sort -u > "$$gene_file"

		{params.bcftools} +split-vep {input.vcf} \
			-c SYMBOL \
			-i "SYMBOL=@$$gene_file" \
			-d -x \
			-O z -o {output.vcf} 2>>{log} || exit 1
	"""


def _parse_vcf_input(wc):
	if config['ngs_type'] == 'panel':
		return rules.r9_3_filter_capture_genes.output.vcf
	elif config['ngs_type'] == 'WES':
		return rules.r9_2_filter_panel_genes.output.vcf
	else:
		return rules.r9_1_split_joint_vcf.output.vcf \
			if config.get('calling_mode', 'joint') == 'joint' \
			else rules.r8_4_vep_germline_individual.output.vcf


rule r9_4_parse_vcf:
	"""Parse panel-filtered or capture-filtered VCF to TSV."""
	wildcard_constraints:
		sample     = "|".join(ngs.GRM_SAMPLES),
		panel_name = "|".join(_panel_names) if _panel_names else "NOPANEL"
	input:
		vcf = _parse_vcf_input
	output:
		tsv = "results/{run}/germline/tsv/{sample}.{panel_name}.unfiltered.tsv"
	params:
		bcftools = config['tools']['bcftools']
	priority: 40
	log:
		'results/{run}/logs/parse_results/{sample}.{panel_name}.parse_vcf.log'
	shell: """
		echo -e 'Chr\\tRef\\tAlt\\tRefGene\\tExon\\tHGVS_description\\tZyg\\t\
rsID\\tGATK_FILTER\\tConsequence\\tConsequence_AA\\t\
gnomAD_exome_NFE\\tgnomAD_exome_Comb\\tgnomAD_genome_NFE\\tgnomAD_genome_Comb\\t\
CADD_PHRED\\tCADD_RAW\\tClinVar_CLNSIG\\tClinVar_CLNREVSTAT\\tSIFT\\tPolyPhen\\tExAC_pLI\\tPUBMED\\t\
ClinVar_publications\\tClinVar_CLNDN\\tBIOTYPE\\tCANONICAL\\tPHENOTYPES\\t\
am_class\\tam_pathogenicity\\tSNPred_score\\t\
REVEL\\tMPC\\tSpliceAI_DS_AG\\tSpliceAI_DS_AL\\tSpliceAI_DS_DG\\tSpliceAI_DS_DL\\tNMD\\t\
CoverageDepth\\tGenotypeQual\\tAlleleDepth' > {output.tsv}

		{params.bcftools} view -m2 -M2 {input.vcf} | \
		{params.bcftools} +split-vep -s primary -d -f \
'%CHROM:%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%EXON\\t%SYMBOL;%HGVSg;%HGVSc;%HGVSp\\t[%GT]\\t\
%Existing_variation\\t%FILTER\\t%Consequence\\t%Amino_acids\\t\
%gnomADe_NFE_AF\\t%gnomADe_AF\\t%gnomADg_NFE_AF\\t%gnomADg_AF\\t\
%CADD_PHRED\\t%CADD_RAW\\t%CLIN_SIG\\t%ClinVar_CLNREVSTAT\\t%SIFT\\t%PolyPhen\\t%pLI_gene_value\\t%PUBMED\\t\
%ClinVar\\t%ClinVar_CLNDN\\t%BIOTYPE\\t%CANONICAL\\t%PHENOTYPES\\t\
%am_class\\t%am_pathogenicity\\t%SNPred_SNPred_score\\t\
%REVEL\\t%MPC\\t%SpliceAI_pred_DS_AG\\t%SpliceAI_pred_DS_AL\\t%SpliceAI_pred_DS_DG\\t%SpliceAI_pred_DS_DL\\t%NMD\\t\
[%DP]\\t[%GQ]\\t[%AD]\\n' >> {output.tsv} 2>>{log}
	"""


rule r9_5_parse_vcf_wes_clinical:
	"""Parse full WES VCF (no panel gene filter) to TSV."""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		vcf = rules.r9_1_split_joint_vcf.output.vcf \
			if config.get('calling_mode', 'joint') == 'joint' \
			else rules.r8_4_vep_germline_individual.output.vcf
	output:
		tsv = "results/{run}/germline/tsv/{sample}.wes_clinical.unfiltered.tsv"
	params:
		bcftools = config['tools']['bcftools']
	priority: 40
	log:
		'results/{run}/logs/parse_results/{sample}.wes_clinical.parse_vcf.log'
	shell: """
		echo -e 'Chr\\tRef\\tAlt\\tRefGene\\tExon\\tHGVS_description\\tZyg\\t\
rsID\\tGATK_FILTER\\tConsequence\\tConsequence_AA\\t\
gnomAD_exome_NFE\\tgnomAD_exome_Comb\\tgnomAD_genome_NFE\\tgnomAD_genome_Comb\\t\
CADD_PHRED\\tCADD_RAW\\tClinVar_CLNSIG\\tClinVar_CLNREVSTAT\\tSIFT\\tPolyPhen\\tExAC_pLI\\tPUBMED\\t\
ClinVar_publications\\tClinVar_CLNDN\\tBIOTYPE\\tCANONICAL\\tPHENOTYPES\\t\
am_class\\tam_pathogenicity\\tSNPred_score\\t\
REVEL\\tMPC\\tSpliceAI_DS_AG\\tSpliceAI_DS_AL\\tSpliceAI_DS_DG\\tSpliceAI_DS_DL\\tNMD\\t\
CoverageDepth\\tGenotypeQual\\tAlleleDepth' > {output.tsv}

		{params.bcftools} view -m2 -M2 {input.vcf} | \
		{params.bcftools} +split-vep -s primary -d -f \
'%CHROM:%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%EXON\\t%SYMBOL;%HGVSg;%HGVSc;%HGVSp\\t[%GT]\\t\
%Existing_variation\\t%FILTER\\t%Consequence\\t%Amino_acids\\t\
%gnomADe_NFE_AF\\t%gnomADe_AF\\t%gnomADg_NFE_AF\\t%gnomADg_AF\\t\
%CADD_PHRED\\t%CADD_RAW\\t%CLIN_SIG\\t%ClinVar_CLNREVSTAT\\t%SIFT\\t%PolyPhen\\t%pLI_gene_value\\t%PUBMED\\t\
%ClinVar\\t%ClinVar_CLNDN\\t%BIOTYPE\\t%CANONICAL\\t%PHENOTYPES\\t\
%am_class\\t%am_pathogenicity\\t%SNPred_SNPred_score\\t\
%REVEL\\t%MPC\\t%SpliceAI_pred_DS_AG\\t%SpliceAI_pred_DS_AL\\t%SpliceAI_pred_DS_DG\\t%SpliceAI_pred_DS_DL\\t%NMD\\t\
[%DP]\\t[%GQ]\\t[%AD]\\n' >> {output.tsv} 2>>{log}
	"""


rule r9_6_filter_tsv:
	"""Postprocess panel-filtered TSV. gnomad filter disabled for gene panels."""
	wildcard_constraints:
		sample     = "|".join(ngs.GRM_SAMPLES),
		panel_name = "|".join(_panel_names) if _panel_names else "NOPANEL"
	input:
		tsv = rules.r9_4_parse_vcf.output.tsv,
		bed = lambda wc: _wes_panels.get(wc.panel_name, config['panel_capture']['target'])
	output:
		tsv = temp("results/{run}/germline/tsv/{sample}.{panel_name}.tsv")
	priority: 40
	params:
		mart                = config['references']['vep_plugins_data']['custom']['mart'],
		rank                = config['references']['vep_plugins_data']['custom']['rank'],
		gnomad_filter       = False,
		gnomad_af_threshold = config.get('filters', {}).get('gnomad', {}).get('af_threshold', 0.01),
		gnomad_column       = config.get('filters', {}).get('gnomad', {}).get('column', 'gnomAD_exome_NFE'),
		bed_file            = lambda wc: _wes_panels.get(wc.panel_name, config['panel_capture']['target']),
		ngs_type            = config['ngs_type'],
		consequence_filter  = True
	log:
		'results/{run}/logs/parse_results/{sample}.{panel_name}.filter_tsv.log'
	script: "scripts/postprocess.py"


rule r9_7_filter_tsv_wes_clinical:
	"""Postprocess full WES TSV. gnomad filter enabled at 0.01."""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		tsv = rules.r9_5_parse_vcf_wes_clinical.output.tsv,
		bed = config['panel_capture']['target']
	output:
		tsv = temp("results/{run}/germline/tsv/{sample}.wes_clinical.tsv")
	priority: 40
	params:
		mart                = config['references']['vep_plugins_data']['custom']['mart'],
		rank                = config['references']['vep_plugins_data']['custom']['rank'],
		gnomad_filter       = True,
		gnomad_af_threshold = 0.01,
		gnomad_column       = config.get('filters', {}).get('gnomad', {}).get('column', 'gnomAD_exome_NFE'),
		bed_file            = config['panel_capture']['target'],
		ngs_type            = config['ngs_type'],
		consequence_filter  = True
	log:
		'results/{run}/logs/parse_results/{sample}.wes_clinical.filter_tsv.log'
	script: "scripts/postprocess.py"


rule r9_8_merge_tsv_to_xlsx_panel:
	"""Merge all samples for one WES gene panel into a single xlsx with Legend."""
	wildcard_constraints:
		panel_name = "|".join(_panel_names) if _panel_names else "NOPANEL"
	input:
		tsv = expand(
			"results/{run}/germline/tsv/{sample}.{panel_name}.tsv",
			run=config['run'], sample=ngs.GRM_SAMPLES, allow_missing=True
		)
	output:
		xlsx = "results/{run}/germline/xlsx/panel.{panel_name}.{run}.germline.results.xlsx"
	params:
		tsv            = "results/{run}/run_table.tsv",
		hs_metrics_dir = f"results/{config['run']}/bam/hs_metrics",
		vcf_dir        = f"results/{config['run']}/germline/vcf",
		bcftools       = config['tools']['bcftools']
	log:
		'results/{run}/logs/parse_results/{panel_name}.tsv2xlsx.log'
	script:
		"scripts/collect_tsv_to_xml.py"


rule r9_9_merge_tsv_to_xlsx_wes_clinical:
	"""Merge all samples for full WES into a single xlsx with Legend."""
	input:
		tsv = expand(
			"results/{run}/germline/tsv/{sample}.wes_clinical.tsv",
			run=config['run'], sample=ngs.GRM_SAMPLES
		)
	output:
		xlsx = "results/{run}/germline/xlsx/wes_clinical.{run}.germline.results.xlsx"
	params:
		tsv            = "results/{run}/run_table.tsv",
		hs_metrics_dir = f"results/{config['run']}/bam/hs_metrics",
		vcf_dir        = f"results/{config['run']}/germline/vcf",
		bcftools       = config['tools']['bcftools']
	log:
		'results/{run}/logs/parse_results/wes_clinical.tsv2xlsx.log'
	script:
		"scripts/collect_tsv_to_xml.py"


rule r9_10_merge_tsv_to_xlsx:
	"""Merge all samples for panel ngs_type or WGS into a single xlsx with Legend."""
	input:
		tsv = expand("results/{run}/germline/tsv/{sample}.tsv", run=config['run'], sample=ngs.GRM_SAMPLES)
	output:
		xlsx = "results/{run}/germline/xlsx/individual.{run}.germline.results.xlsx"
	params:
		tsv            = "results/{run}/run_table.tsv",
		hs_metrics_dir = f"results/{config['run']}/bam/hs_metrics",
		vcf_dir        = f"results/{config['run']}/germline/vcf",
		bcftools       = config['tools']['bcftools']
	log:
		'results/{run}/logs/parse_results/tsv2xlsx.log'
	script:
		"scripts/collect_tsv_to_xml.py"


rule r9_11_parse_vcf_raw:
	"""Parse full annotated VCF to TSV without any filtering.
	Preserves all variants for archival purposes.
	"""
	wildcard_constraints:
		sample = "|".join(ngs.GRM_SAMPLES)
	input:
		vcf = rules.r9_1_split_joint_vcf.output.vcf \
			if config.get('calling_mode', 'joint') == 'joint' \
			else rules.r8_4_vep_germline_individual.output.vcf
	output:
		tsv = temp("results/{run}/germline/tsv/{sample}.raw.tsv")
	params:
		bcftools = config['tools']['bcftools']
	priority: 40
	resources:
		mem_mb      = 2000,
		runtime_min = 60
	log:
		'results/{run}/logs/parse_results/{sample}.raw.parse_vcf.log'
	shell: """
		echo -e 'Chr\tRef\tAlt\tRefGene\tExon\tHGVS_description\tZyg\trsID\tGATK_FILTER\tConsequence\tConsequence_AA\tgnomAD_exome_NFE\tgnomAD_exome_Comb\tgnomAD_genome_NFE\tgnomAD_genome_Comb\tCADD_PHRED\tCADD_RAW\tClinVar_CLNSIG\tClinVar_CLNREVSTAT\tSIFT\tPolyPhen\tExAC_pLI\tPUBMED\tClinVar_publications\tClinVar_CLNDN\tBIOTYPE\tCANONICAL\tPHENOTYPES\tam_class\tam_pathogenicity\tSNPred_score\tREVEL\tMPC\tSpliceAI_DS_AG\tSpliceAI_DS_AL\tSpliceAI_DS_DG\tSpliceAI_DS_DL\tNMD\tCoverageDepth\tGenotypeQual\tAlleleDepth' > {output.tsv}

		{params.bcftools} view -m2 -M2 {input.vcf} | \
		{params.bcftools} +split-vep -s primary -d -f '%CHROM:%POS\t%REF\t%ALT\t%SYMBOL\t%EXON\t%SYMBOL;%HGVSg;%HGVSc;%HGVSp\t[%GT]\t%Existing_variation\t%FILTER\t%Consequence\t%Amino_acids\t%gnomADe_NFE_AF\t%gnomADe_AF\t%gnomADg_NFE_AF\t%gnomADg_AF\t%CADD_PHRED\t%CADD_RAW\t%CLIN_SIG\t%ClinVar_CLNREVSTAT\t%SIFT\t%PolyPhen\t%pLI_gene_value\t%PUBMED\t%ClinVar\t%ClinVar_CLNDN\t%BIOTYPE\t%CANONICAL\t%PHENOTYPES\t%am_class\t%am_pathogenicity\t%SNPred_SNPred_score\t%REVEL\t%MPC\t%SpliceAI_pred_DS_AG\t%SpliceAI_pred_DS_AL\t%SpliceAI_pred_DS_DG\t%SpliceAI_pred_DS_DL\t%NMD\t[%DP]\t[%GQ]\t[%AD]\n' >> {output.tsv} 2>>{log}
	"""