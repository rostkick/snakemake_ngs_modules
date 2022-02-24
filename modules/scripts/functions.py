from ast import literal_eval
from itertools import permutations


def get_min_substr(samples):
	series_list = samples.str.split(r'_|\.|-')
	
	i=0
	condition = '_'.join(series_list[0][:i]) == '_'.join(series_list[len(series_list)-1][:i])
	while condition:
		substr_set = set(series_list.str[:i].str.join('_').tolist())
		if len(substr_set) == len(set(samples.tolist())):
			break
		i+=1
	series = series_list.str[:i].str.join('_')
	return series


def apply_literal_eval(*l, data):
	"""Plain string to normal python list"""
	for i in l:
		print(data[i])
		data[i] = [literal_eval(i) for i in data[i]]


def get_somatic_input(wc, data):
	sample_tumor = data.loc[:, 'tumor_samples'][data.loc[:, 'patient']==wc.patient].values[0]
	sample_germline = data.loc[:, 'germline_samples'][data.loc[:, 'patient']==wc.patient].values[0]
	return {'tumor': f'results/{wc.run}/bam/{sample_tumor}.final.bam',
			'germline': f'results/{wc.run}/bam/{sample_germline}.final.bam'}


def get_inputs(run, mode, patients=None):
	if ('somatic' in mode) | (mode == 'all'):
		final_output_soma = [f'results/{run}/somatic/{patient}/plots/plot.pdf' for patient in patients] + \
							[f'results/{run}/somatic/{patient}/annotation/somatic.annotated.vcf.gz' for patient in patients]
	final_output_germ = [f'results/{run}/germline/vcf/germline.annotated.vcf.gz']
	final_output_sv = [f'results/{run}/germline/sv/sv.annotated.vcf.gz']


	possible_mode_spelling_all = set([','.join(i) for i in list(permutations(['germline', 'somatic', 'sv']))] + ['all'])
	possible_mode_spelling_soma_germ = set([','.join(i) for i in list(permutations(['germline', 'somatic']))])
	possible_mode_spelling_germ_sv = set([','.join(i) for i in list(permutations(['germline', 'sv']))])
	possible_mode_spelling_soma_sv = set([','.join(i) for i in list(permutations(['somatic', 'sv']))])

	if mode in possible_mode_spelling_all:
		return final_output_soma + final_output_germ + final_output_sv
	elif mode in possible_mode_spelling_soma_germ:
		return final_output_soma + final_output_germ
	elif mode in possible_mode_spelling_germ_sv:
		return final_output_germ + final_output_sv
	elif mode in possible_mode_spelling_soma_sv:
		return final_output_soma + final_output_sv
	elif mode == 'somatic':
		return final_output_soma
	elif mode == 'germline':
		return final_output_germ
	elif mode == 'sv':
		return final_output_sv
	else:
		error_message = "\x1b[0;31;40m" +"incorrect mode spelling!"+"\x1b[0m"
		assert False, error_message
