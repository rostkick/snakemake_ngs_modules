import os
from glob import glob
from difflib import SequenceMatcher
import re
from itertools import product
import pandas as pd
import numpy as np
from snakemake.workflow import config


RUN = config['run']
GERMLINE = True if config['grm_dir'] is not None else False
SOMATIC = True if config['tmr_dir'] is not None else False
PAIR = True if config['reads_type'] == 'pair' else False

class FastqDirPath:

	def __init__(self, fastq_dir_path: str) -> None:
		self.fastq_dir_path = os.path.abspath(fastq_dir_path)


class FastqList(FastqDirPath):

	def __init__(self, fastq_dir_path: str) -> None:
		super().__init__(fastq_dir_path)
		self.fastqs = glob(self.fastq_dir_path+'/*.fastq.gz')
		self.fastq_f = [fq for fq in self.fastqs if re.match('.*_R1_.*', fq) is not None] if PAIR else self.fastqs # mb None is better
		self.fastq_r = [fq for fq in self.fastqs if re.match('.*_R2_.*', fq) is not None] if PAIR else None

class SampleSeries(FastqList):

	def __init__(self, fastq_dir_path: str) -> None:
		super().__init__(fastq_dir_path)
		self.sample = [os.path.basename(fastq) for fastq in self.fastq_f]
		self.sample = pd.Series(self.sample)
		self.__get_min_sample_name()

	def __get_min_sample_name(self):
		if len(self.sample) > 1:
			series_list = self.sample.str.split(r'_|\.|-')
			i=0
			condition = '_'.join(series_list[0][:i]) == '_'.join(series_list[len(series_list)-1][:i])
			while condition:
				substr_set = set(series_list.str[:i].str.join('_').tolist())
				if len(substr_set) == len(set(self.sample.tolist())):
					break
				i+=1
			self.sample = series_list.str[:i].str.join('_')
		else:
			self.sample = self.sample.str.extract('(.+?)[_\-\.\|].*')[0]


class SampleDataframe(SampleSeries):

	def __init__(self, fastq_dir_path: str) -> None:
		super().__init__(fastq_dir_path)
		self.sample_df = pd.DataFrame()
		self.__create_sample_data_frame()

	def __create_sample_data_frame(self):
		self.sample_df['sample'] = self.sample
		self.sample_df['fastq_f'] = self.fastq_f
		self.sample_df['fastq_r'] = self.fastq_r
		self.sample_df.dropna(axis=1, inplace=True)


class MapNames:

	def __init__(self, grm_dir_path: str|None, tmr_dir_path: str|None) -> None:
		self.grm = SampleDataframe(grm_dir_path) if GERMLINE else None
		self.tmr = SampleDataframe(tmr_dir_path) if SOMATIC else None
		self.distances = None
		self.map_df = None
		self.__estimate_distances()
		self.__map_names()

	def __estimate_distances(self):
		if GERMLINE and SOMATIC:
			distances = []
			for g, t in product(self.grm.sample.tolist(), self.tmr.sample.tolist()):
				distances.append((g, t, SequenceMatcher(None, g, t).ratio()))
			self.distances = distances

	def __map_names(self):
		if GERMLINE and SOMATIC is not None:
			map_df = pd.DataFrame(self.distances)
			map_df.columns = ['grm_samples', 'tmr_samples', 'distance']
			map_df = map_df.sort_values('distance', ascending=False)
			is_duplicated = map_df.apply(pd.Series.duplicated, axis=0)
			map_df = map_df.where(~is_duplicated, np.nan)
			map_df = map_df.drop('distance', axis=1)
			self.map_df = map_df.dropna(axis=0, how='all')


class MapData(MapNames):

	def __init__(self, grm_dir_path: str|None, tmr_dir_path: str|None) -> None:
		super().__init__(grm_dir_path, tmr_dir_path)
		self.wide_df = pd.DataFrame()
		self.__merge_data()
		self.__add_patient_ids()
		self.__add_suffixes_to_samples()
		self.__remove_suffixes_from_cols()

	def __merge_data(self):
		self.wide_df = pd.merge(self.map_df, self.grm.sample_df, how='outer', left_on='grm_samples', right_on='sample')
		self.wide_df = self.wide_df.drop('sample', axis=1)
		self.wide_df.columns = ['grm_' + str(col) if 'fastq' in col else col for col in self.wide_df.columns]
		self.wide_df = pd.merge(self.wide_df, self.tmr.sample_df, how='outer', left_on='tmr_samples', right_on='sample')
		self.wide_df = self.wide_df.drop('sample', axis=1)
		self.wide_df.columns = ['tmr_' + str(col) if col.startswith('fastq') else col for col in self.wide_df.columns]

	def __add_patient_ids(self):
		if self.wide_df['tmr_samples'].isna().sum() < self.wide_df['grm_samples'].isna().sum():
			self.wide_df['patients'] = self.wide_df['tmr_samples']
		else:
			self.wide_df['patients'] = self.wide_df['grm_samples']
		
	def __add_suffixes_to_samples(self):
		self.wide_df['grm_samples'] = self.wide_df['grm_samples'] + '_grm'
		self.wide_df['tmr_samples'] = self.wide_df['tmr_samples'] + '_tmr'

	def __remove_suffixes_from_cols(self):
		if PAIR is False:
			self.wide_df.columns = self.wide_df.columns.str.rstrip('_f')


class NGSData(MapData):
	def __init__(self, grm_dir_path: str|None, tmr_dir_path: str|None) -> None:
		super().__init__(grm_dir_path, tmr_dir_path)
		self.long_df = pd.DataFrame()
		self.__create_long_dataframe()

	def __create_long_dataframe(self):
		if PAIR:
			self.long_df = pd.DataFrame(
				{
				"samples": np.concatenate([self.wide_df.loc[:, 'tmr_samples'].values,
										self.wide_df.loc[:, 'grm_samples'].values]),
				"fastq_forward": np.concatenate([self.wide_df.loc[:, "tmr_fastq_f"].values,
												self.wide_df.loc[:, "grm_fastq_f"].values]),
				"fastq_reverse": np.concatenate([self.wide_df.loc[:, "tmr_fastq_r"].values,
												self.wide_df.loc[:, "grm_fastq_r"].values])
				}
										)
			self.long_df = self.long_df.dropna(axis=0, how='all')
		else:
			self.long_df = pd.DataFrame(
				{
				"samples": np.concatenate([self.wide_df.loc[:, 'tmr_samples'].values,
										self.wide_df.loc[:, 'grm_samples'].values]),
				"fastq_forward": np.concatenate([self.wide_df.loc[:, "tmr_fastq"].values,
										self.wide_df.loc[:, "grm_fastq"].values])})
