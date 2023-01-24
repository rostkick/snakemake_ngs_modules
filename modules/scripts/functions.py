from modules.scripts.params_builder import RUN, GERMLINE, SOMATIC


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


def get_somatic_input(wc, data):
	sample_tumor = data.loc[:, 'tmr_samples'][data.loc[:, 'patients']==wc.patient].values[0]
	sample_germline = data.loc[:, 'grm_samples'][data.loc[:, 'patients']==wc.patient].values[0]
	return {'tumor': f'results/{wc.run}/bam/{sample_tumor}.final.bam',
			'germline': f'results/{wc.run}/bam/{sample_germline}.final.bam'}

