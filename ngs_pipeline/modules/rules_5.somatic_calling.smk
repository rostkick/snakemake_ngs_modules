from collections import defaultdict


def get_somatic_input(wc, data):
	inp = defaultdict(str)
	sample_tumor = wc.patient + '_tmr'
	inp['tumor'] = f'results/{wc.run}/bam/{sample_tumor}.final.bam'
	if ngs.GRM:
		sample_germline = wc.patient + '_grm'
		inp['germline'] = f'results/{wc.run}/bam/{sample_germline}.final.bam'
	return inp

rule r5_1_collect_fir2_counts:
	input: 
		bam = lambda wc: get_somatic_input(wc, ngs.data)['tumor']
	output: 
		f1r2 = 'results/{run}/somatic/{patient}/f1r2.tsv'
	log: 
		'results/{run}/logs/somatic/{patient}/CollectF1R2Counts.log'
	params:
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
	shell:"""
		gatk CollectF1R2Counts \
				-R {params.ref} \
				-I {input.bam} \
				-O {output.f1r2} 2>{log}"""


rule r5_2_learn_read_orientation_model:
	input: 
		f1r2 = rules.r5_1_collect_fir2_counts.output.f1r2
	output: 
		rom = temp('results/{run}/somatic/{patient}/read-orientation-model.tar.gz')
	log: 
		'results/{run}/logs/somatic/{patient}/LearnReadOrientationModel.log'
	shell: """
			gatk LearnReadOrientationModel \
				-I {input.f1r2} \
				-O {output.rom} 2>{log}"""

rule r5_3_get_pileup_summaries_tmr:
	input: 
		bam = lambda wc: get_somatic_input(wc, ngs.data)['tumor']
	output: 
		getpileupsum = 'results/{run}/somatic/{patient}/getpileupsummaries_tmr.table'
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_tumor.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		exac = config['references38']['small_exac_common'] if config['assembly'] == 'GRCh38' else config['references37']['small_exac_common']
	shell: """
			{params.gatk} GetPileupSummaries \
				-I {input.bam} \
				-R {params.ref} \
				-V {params.exac} \
				-L {params.exac} \
				-O {output.getpileupsum} 2>{log}"""

# ===============================================
# ===== Somatic vs germline variant calling =====
# ===============================================
rule r5_4_mutect2_grm_vs_tmr:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		bam_tmr = lambda wc: get_somatic_input(wc, ngs.data)['tumor'],
		bam_grm = lambda wc: get_somatic_input(wc, ngs.data)['germline'],
		capture = rules.r2_8_bed_to_intervals.output.intervals if config['panel_capture']['target'][-4:] == 'bed' else config['panel_capture']['target']
	output: 
		vcf = temp('results/{run}/somatic/{patient}/raw.vcf.gz'),
		bam = temp('results/{run}/somatic/{patient}/raw.bam')
	log: 
		'results/{run}/logs/somatic/{patient}/Mutect2.log'
	params:
		gatk = config['tools']['gatk'],
		sample_name = lambda wc: wc.patient + '_grm',
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		grm_res = config['references38']['af_only_gnomad'] if config['assembly'] == 'GRCh38' else config['references37']['af_only_gnomad'],
		pon = config['references38']['snps'] if config['assembly'] == 'GRCh38' else config['references37']['snps']
	shell:
		"""{params.gatk} Mutect2 \
				-R {params.ref} \
				-I {input.bam_tmr} \
				-O {output.vcf} \
				-I {input.bam_grm} \
				-normal {params.sample_name} \
				-bamout {output.bam} \
				--germline-resource {params.grm_res} \
				--panel-of-normals {params.pon} 2>{log}"""

use rule r5_3_get_pileup_summaries_tmr as r5_5_get_pileup_summaries_grm with:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		bam = lambda wc: get_somatic_input(wc, ngs.data)['germline']
	output: 
		getpileupsum = temp('results/{run}/somatic/{patient}/getpileupsummaries_grm.table')
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_germline.log'


rule r5_6_calculate_contamination_grm_vs_tmr:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		getpileupsum_tmr = rules.r5_3_get_pileup_summaries_tmr.output.getpileupsum,
		getpileupsum_grm = rules.r5_5_get_pileup_summaries_grm.output.getpileupsum
	output: 
		contamination = 'results/{run}/somatic/{patient}/contamination.table'
	params:
		gatk = config['tools']['gatk']
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell:"""
		{params.gatk} CalculateContamination \
				-I {input.getpileupsum_tmr} \
				-O {output.contamination} \
				-matched {input.getpileupsum_grm} 2>{log}"""

rule r5_7_filter_mutect_calls_grm_vs_tmr:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		vcf = rules.r5_4_mutect2_grm_vs_tmr.output.vcf,
		contamination = rules.r5_6_calculate_contamination_grm_vs_tmr.output.contamination,
		rom = rules.r5_2_learn_read_orientation_model.output.rom
	output: 
		vcf = 'results/{run}/somatic/{patient}/filtered.vcf.gz'
	log: 
		'results/{run}/logs/somatic/{patient}/FilterMutectCalls.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa']
	shell:"""
			{params.gatk} FilterMutectCalls \
				-R {params.ref} \
				-V {input.vcf} \
				-O {output.vcf} \
				--contamination-table {input.contamination} \
				--orientation-bias-artifact-priors {input.rom} 2>{log}
				"""

rule r5_8_filter_pass_exclude_normal_grm_vs_tmr:
	wildcard_constraints:
		patient = "|".join(ngs.GRM_VS_TMR_PATIENTS)
	input: 
		vcf = rules.r5_7_filter_mutect_calls_grm_vs_tmr.output.vcf
	output: 
		vcf = 'results/{run}/somatic/{patient}/final.vcf.gz'
	params: 
		bcftools = config['tools']['bcftools']
	# threads: 
	# 	workflow.cores/max(len(ngs.GRM_VS_TMR_PATIENTS), 1)
	shell:"""
			tumor_sample=$(zcat {input.vcf} | \
					grep -oP '##tumor_sample=\K.+') &&\
					{params.bcftools} view -i "%FILTER='PASS'" \
					-s $tumor_sample \
					-Oz -o {output.vcf} {input.vcf}"""

# ===========================================
# ===== Tumor-only mode variant calling =====
# ===========================================
rule r5_9_mutect2_tumor_only:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		bam = lambda wc: get_somatic_input(wc, ngs.data)['tumor'],
		capture = rules.r2_8_bed_to_intervals.output.intervals if config['panel_capture']['target'][-4:] == 'bed' else config['panel_capture']['target']
	output: 
		vcf = 'results/{run}/somatic/{patient}/raw_tonly.vcf.gz',
		bam = 'results/{run}/somatic/{patient}/raw_tonly.bam'
	log: 
		'results/{run}/logs/somatic/{patient}/Mutect2.log'
	params:
		gatk = config['tools']['gatk'],
		ref = config['references38']['genome_fa'] if config['assembly'] == 'GRCh38' else config['references37']['genome_fa'],
		grm_res = config['references38']['af_only_gnomad'] if config['assembly'] == 'GRCh38' else config['references37']['af_only_gnomad'],
		pon = config['references38']['snps'] if config['assembly'] == 'GRCh38' else config['references37']['snps']
	shell: """
			{params.gatk} Mutect2 \
				-R {params.ref} \
				-I {input.bam} \
				-O {output.vcf} \
				-bamout {output.bam} \
				--germline-resource {params.grm_res} \
				--panel-of-normals {params.pon} 2>{log}"""

rule r5_10_calculate_contamination_tmr_only:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		getpileupsum_tmr = rules.r5_3_get_pileup_summaries_tmr.output.getpileupsum
	output: 
		contamination = 'results/{run}/somatic/{patient}/contamination_tonly.table'
	params:
		gatk = config['tools']['gatk']
	log: 
		'results/{run}/logs/somatic/{patient}/CalculateContamination.log'
	shell:"""
			{params.gatk} CalculateContamination \
				-I {input.getpileupsum_tmr} \
				-O {output.contamination} 2>{log}"""


use rule r5_7_filter_mutect_calls_grm_vs_tmr as r5_11_filter_mutect_calls_tmr_only with:
	wildcard_constraints:
		patient = "|".join(ngs.ONLY_TMR_PATIENTS)
	input: 
		vcf = rules.r5_9_mutect2_tumor_only.output.vcf,
		contamination = rules.r5_10_calculate_contamination_tmr_only.output.contamination,
		rom = rules.r5_2_learn_read_orientation_model.output.rom
	output:
		vcf = 'results/{run}/somatic/{patient}/final_tonly.vcf.gz'
