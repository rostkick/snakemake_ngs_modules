import os
from itertools import product
import pandas as pd
import numpy as np
from glob import glob
from difflib import SequenceMatcher
from snakemake.workflow import config


def trim_sample_name(df):
	if len(df['sample']) > 1:
		series_list = df['sample'].str.split(r'_|\.|-')
		i=0
		condition = '_'.join(series_list[0][:i]) == '_'.join(series_list[len(series_list)-1][:i])
		while condition:
			substr_set = set(series_list.str[:i].str.join('_').tolist())
			if len(substr_set) == len(set(df['sample'].tolist())):
				break
			i+=1
		df['sample'] = series_list.str[:i].str.join('_')
	else:
		df['sample'] = df['sample'].str.extract('(.+?)[_\-\.\|].*')[0]
	return df


def get_fqs_table(dir_path, suffix):
	fqs = list(glob(dir_path+'/*fq.gz'))
	df = pd.DataFrame({'fastq': fqs})
	df['base_fastq'] = df['fastq'].apply(lambda x: os.path.basename(x))
	print(df)
	df_extracted = df['base_fastq'].str.extractall(r'(?P<sample>.*)[_\-\.](?P<lane>[lL]\d*)[_\-\.](?P<reads_orientation>[Rr][12])[_\-\.]?.*?').droplevel(1)
	df = pd.concat([df, df_extracted], axis=1)
	df = df.drop('base_fastq', axis=1)
	df = trim_sample_name(df)
	df['sample'] = df['sample'] + suffix
	return df

def estimate_distances(grm: list, tmr: list) -> list:
	distances = []
	for g, t in product(grm, tmr):
		distances.append((g, t, SequenceMatcher(None, g, t).ratio()))
	return distances

def map_names(mapping: list) -> pd.DataFrame:
	df = pd.DataFrame(mapping)
	df.columns = ['sample_grm', 'sample_tmr', 'distance']

	df = df.sort_values('distance', ascending=False)
	df_matched = df.loc[df['distance']>=0.5, :]
	
	df_grm_unmatched = pd.DataFrame()
	df_grm_unmatched['sample_grm'] = df.loc[df['distance']<0.5, 'sample_grm'].drop_duplicates().dropna()
	df_grm_unmatched['sample_tmr'] = np.nan
	df_grm_unmatched['distance'] = np.nan
	
	df_tmr_unmatched = pd.DataFrame()
	df_tmr_unmatched['sample_grm'] = np.nan
	df_tmr_unmatched['sample_tmr'] = df.loc[df['distance']<0.5, 'sample_tmr'].drop_duplicates().dropna()
	df_tmr_unmatched['distance'] = np.nan
	
	df = pd.concat([df_matched, df_grm_unmatched, df_tmr_unmatched])
	is_duplicated = df.apply(pd.Series.duplicated, axis=0)
	df = df.where(~is_duplicated, np.nan)
	df = df.drop('distance', axis=1)
	mapping = df.dropna(axis=0, how='all')
	return mapping

def get_long_table(grm, tmr, mapping) -> pd.DataFrame:
	wide_df = pd.merge(mapping, grm, how='outer', left_on='sample_grm', right_on='sample')
	wide_df = wide_df.drop('sample', axis=1)
	wide_df = pd.merge(wide_df, tmr, how='outer', left_on=['sample_tmr'], right_on=['sample'], suffixes=['_grm', '_tmr'])
	wide_df = wide_df.drop('sample', axis=1)
	long_df = pd.DataFrame({
	"sample": np.concatenate([wide_df.loc[:, 'sample_tmr'].values,
							wide_df.loc[:, 'sample_grm'].values]),
	"fastq": np.concatenate([wide_df.loc[:, "fastq_tmr"].values,
							wide_df.loc[:, "fastq_grm"].values]),
	"lane": np.concatenate([wide_df.loc[:, "lane_tmr"].values,
							wide_df.loc[:, "lane_grm"].values]),
	"reads_orientation": np.concatenate([wide_df.loc[:, "reads_orientation_tmr"].values,
										wide_df.loc[:, "reads_orientation_grm"].values])})
	long_df = long_df.drop_duplicates().dropna(axis=0, how='all')
	return long_df

def setup_run():
	grm_dir = config['grm_dir']
	tmr_dir = config['tmr_dir']
	grm = get_fqs_table(grm_dir, '_grm')
	tmr = get_fqs_table(tmr_dir, '_tmr')
	dist = estimate_distances(grm['sample'].unique(), tmr['sample'].unique())
	mapping = map_names(dist)
	long_df = get_long_table(grm, tmr, mapping)
	mapping['patient'] = mapping['sample_tmr'].str[:-4]
	return mapping, long_df
