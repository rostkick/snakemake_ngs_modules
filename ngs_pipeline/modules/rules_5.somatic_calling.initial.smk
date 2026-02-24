from collections import defaultdict


def get_somatic_input(wc, data):
	inp = defaultdict(str)
	sample_tmr = wc.patient + '_tmr'
	inp['tumor'] = f'results/{wc.run}/bam/{sample_tmr}.final.bam'
	if ngs.GRM:
		sample_germline = wc.patient + '_grm'
		inp['germline'] = f'results/{wc.run}/bam/{sample_germline}.final.bam'
	return inp

rule r5_1_collect_fir2_counts:
	input: 
		bam = lambda wc: get_somatic_input(wc, ngs.data)['tumor']
	output: 
		f1r2 = temp('results/{run}/somatic/{patient}/f1r2.tsv')
	log: 
		'results/{run}/logs/somatic/{patient}/CollectF1R2Counts.log'
	priority: 35
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	resources:
		mem_mb={'panel': 2000, 'WES': 4000, 'WGS': 6000}.get(config['ngs_type'], 4000),
		runtime_min={'panel': 120, 'WES': 720, 'WGS': 1440}.get(config['ngs_type'], 720)
	shell:"""
		{params.gatk} --java-options "{params.java_opts}" CollectF1R2Counts \
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
	priority: 35
	params:
		gatk = config['tools']['gatk'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	resources:
		mem_mb={'panel': 2000, 'WES': 4000, 'WGS': 6000}.get(config['ngs_type'], 4000),
		runtime_min={'panel': 60, 'WES': 240, 'WGS': 480}.get(config['ngs_type'], 240)
	shell: """
			{params.gatk} --java-options "{params.java_opts}" LearnReadOrientationModel \
				-I {input.f1r2} \
				-O {output.rom} 2>{log}"""

rule r5_3_get_pileup_summaries_tmr:
	input: 
		bam = lambda wc: get_somatic_input(wc, ngs.data)['tumor']
	output: 
		getpileupsum = temp('results/{run}/somatic/{patient}/getpileupsummaries_tmr.table')
	log: 
		'results/{run}/logs/somatic/{patient}/GetPileupSummaries_tmr.log'
	priority: 35
	params:
		gatk = config['tools']['gatk'],
		ref = config['references']['genome_fa'],
		exac = config['references']['small_exac_common'],
		java_opts = lambda wc, resources: f"-Xmx{int(resources.mem_mb * 0.8)}m"
	resources:
		mem_mb={'panel': 4000, 'WES': 8000, 'WGS': 12000}.get(config['ngs_type'], 8000),
		runtime_min={'panel': 240, 'WES': 1440, 'WGS': 2880}.get(config['ngs_type'], 1440)
	shell: """
			{params.gatk} --java-options "{params.java_opts}" GetPileupSummaries \
				-I {input.bam} \
				-R {params.ref} \
				-V {params.exac} \
				-L {params.exac} \
				-O {output.getpileupsum} 2>{log}"""
