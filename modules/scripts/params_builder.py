import os
import re
from glob import glob
from abc import ABCMeta, abstractmethod
import pandas as pd
import numpy as np
from itertools import chain, product
from difflib import SequenceMatcher
from dataclasses import dataclass
from snakemake.workflow import config

#####################################
# family of sequencing data classes #
#####################################
class SeqDir:
	def __init__(self, dir_path: str) -> list:
		dir_path = self.abs_seq_path(dir_path)
		self.seqs = self.get_seq_list(dir_path)

	def abs_seq_path(self, dir_path: str) -> str:
		return os.path.abspath(dir_path)

	def get_seq_list(self, dir_path: str) -> list:
		return list(chain(*[glob(dir_path + '/' + ext) for ext in ['*.fastq', '*.fastq.gz', '*.fq', '*.fq.gz']]))

class DataProcessor(metaclass=ABCMeta):
	@abstractmethod
	def __init__(self, path:str) -> pd.DataFrame:
		pass
	@abstractmethod
	def extract_list(self, path: str) -> list:
		pass
	@abstractmethod
	def convert_to_dataframe(self, sample: pd.Series, fastq: list) -> pd.Series:
		pass
	def trim_sample_name(self, df=None) -> pd.DataFrame:
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

class DataProcessorSingle(DataProcessor):
	def __init__(self, path:str) -> pd.DataFrame:
		fastq = self.extract_list(path)
		df = self.convert_to_dataframe(fastq)
		self.df = super().trim_sample_name(df)

	def extract_list(self, path: str) -> list:
		seq = SeqDir(path)
		return seq.seqs

	def convert_to_dataframe(self, fastq):
		df = pd.DataFrame({'fastq': fastq})
		df['base_fastq'] = df['fastq'].apply(lambda x: os.path.basename(x))
		df_extracted = df['base_fastq'].str.extractall(r'(?P<sample>.*)[_\-\.](?P<lane>[lL]\d*)[_\-\.]?.*?').droplevel(1)
		df = pd.concat([df, df_extracted], axis=1)
		df = df.drop('base_fastq', axis=1)
		return df

class DataProcessorPaired(DataProcessorSingle):
	def __init__(self, path: str) -> pd.DataFrame:
		fastq = super().extract_list(path)
		df = self.convert_to_dataframe(fastq)
		self.df = super().trim_sample_name(df)

	def convert_to_dataframe(self, fastq):
		df = pd.DataFrame({'fastq': fastq})
		df['base_fastq'] = df['fastq'].apply(lambda x: os.path.basename(x))
		df_extracted = df['base_fastq'].str.extractall(r'(?P<sample>.*)[_\-\.](?P<lane>[lL]\d*)[_\-\.](?P<reads_orientation>[Rr][12])[_\-\.]?.*?').droplevel(1)
		df = pd.concat([df, df_extracted], axis=1)
		df = df.drop('base_fastq', axis=1)
		return df

class Mapping:
	@staticmethod
	def map_names(df, suffix) -> pd.DataFrame:
		mapping = pd.DataFrame()
		mapping[f'sample{suffix}'] = df['sample'] + suffix
		mapping.loc[:, 'patient'] = df['sample']
		mapping = mapping.drop_duplicates()
		return mapping

class MappedGermTumor:
	DISTANCE_THRESHOLD = 0.5
	def __init__(self, grm: DataProcessor, tmr: DataProcessor):
		distances = self.estimate_distances(grm, tmr)
		self.mapping = self.map_names(distances)

	def estimate_distances(self, grm: DataProcessor, tmr: DataProcessor) -> list:
		distances = []
		for g, t in product(grm.df['sample'].tolist(), tmr.df['sample'].tolist()):
			distances.append((g, t, SequenceMatcher(None, g, t).ratio()))
		return distances

	def map_names(self, mapping: list) -> pd.DataFrame:
		df = pd.DataFrame(mapping)
		df.columns = ['sample_grm', 'sample_tmr', 'distance']

		df = df.sort_values('distance', ascending=False)
		df_matched = df.loc[df['distance']>=self.DISTANCE_THRESHOLD, :]

		df_grm_unmatched = pd.DataFrame()
		df_grm_unmatched.loc[:, 'sample_grm'] = df.loc[df['distance']<self.DISTANCE_THRESHOLD, 'sample_grm'].drop_duplicates().dropna()
		df_grm_unmatched.loc[:, 'sample_tmr'] = np.nan
		df_grm_unmatched.loc[:, 'distance'] = np.nan

		df_tmr_unmatched = pd.DataFrame()
		df_tmr_unmatched['sample_grm'] = np.nan
		df_tmr_unmatched.loc[:, 'sample_tmr'] = df.loc[df['distance']<self.DISTANCE_THRESHOLD, 'sample_tmr'].drop_duplicates().dropna()
		df_tmr_unmatched.loc[:, 'distance'] = np.nan

		df = pd.concat([df_matched, df_grm_unmatched, df_tmr_unmatched])
		is_duplicated = df.apply(pd.Series.duplicated, axis=0)
		df = df.where(~is_duplicated, np.nan)
		df = df.drop('distance', axis=1)
		mapping = df.dropna(axis=0, how='all').reset_index()
		mapping.loc[:, 'sample_grm'] = mapping['sample_grm'] + '_grm'
		mapping.loc[:, 'sample_tmr'] = mapping['sample_tmr'] + '_tmr'
		mapping.loc[:, 'patient'] = mapping['sample_tmr'].str[:-4]
		return mapping

##################################
# family of wide-format ngs data #
##################################
class NGS(metaclass=ABCMeta):
	@abstractmethod
	def __init__(self):
		pass
	@abstractmethod
	def get_data_processor(self, path) -> DataProcessor:
		pass

class NGSIndividualPaired(NGS):
	def __init__(self, path):
		self.mapping = pd.DataFrame()
		self.data = self.get_data_processor(path)

	def get_data_processor(self, path: str) -> DataProcessorPaired:
		return DataProcessorPaired(path).df


class NGSIndividualSingle(NGSIndividualPaired):
	def get_data_processor(self, path: str) -> DataProcessorSingle:
		return DataProcessorSingle(path).df

class NGSJointPaired(NGS):
	def __init__(self, grm_path, tmr_path):
		grm, tmr = self.get_data_processor(grm_path, tmr_path)
		self.mapping = self.get_mapping(grm, tmr)
		wide_df = self.get_wide_df(self.mapping, grm, tmr)
		self.data = self.get_long_df(wide_df)

	def get_data_processor(self, grm_path: str, tmr_path: str) -> DataProcessorPaired:
		return DataProcessorPaired(grm_path), DataProcessorPaired(tmr_path)

	def get_mapping(self, grm: DataProcessor, tmr: DataProcessor) -> pd.DataFrame:
		return MappedGermTumor(grm, tmr).mapping

	def get_wide_df(self, mapping, grm, tmr):
		wide_df = pd.merge(mapping, grm.df, how='outer', left_on='sample_grm', right_on='sample')
		wide_df = wide_df.drop('sample', axis=1)
		wide_df = pd.merge(wide_df, tmr.df, how='outer', left_on=['sample_tmr'], right_on=['sample'], suffixes=['_grm', '_tmr'])
		wide_df = wide_df.drop('sample', axis=1)
		return wide_df

	def get_long_df(self, wide_df):
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

class NGSJointSingle(NGSJointPaired):
	def get_long_df(self, wide_df):
		long_df = pd.DataFrame({
		"sample": np.concatenate([wide_df.loc[:, 'sample_tmr'].values,
								wide_df.loc[:, 'sample_grm'].values]),
		"fastq": np.concatenate([wide_df.loc[:, "fastq_tmr"].values,
								wide_df.loc[:, "fastq_grm"].values]),
		"lane": np.concatenate([wide_df.loc[:, "lane_tmr"].values,
								wide_df.loc[:, "lane_grm"].values])})
		long_df = long_df.drop_duplicates().dropna(axis=0, how='all')
		return long_df

####################################
# family of high-level ngs classes #
####################################
class Germline:
	def __init__(self, data: pd.DataFrame):
		self.SAMPLES = data['sample'].dropna().unique().tolist()
		self.GRM_SAMPLES = [i for i in self.SAMPLES if '_grm' in i]
		self.GRM_SAMPLES = self.SAMPLES if len(self.GRM_SAMPLES) == 0 else self.GRM_SAMPLES

class Tumor:
	def __init__(self, data: pd.DataFrame):
		self.SAMPLES = data['sample'].dropna().unique().tolist()
		self.TMR_SAMPLES = [i for i in self.SAMPLES if '_tmr' in i]
		self.TMR_SAMPLES = self.SAMPLES if len(self.TMR_SAMPLES) == 0 else self.TMR_SAMPLES
		self.ALL_PATIENTS = self.TMR_SAMPLES.copy()
		self.ONLY_TMR_PATIENTS = self.TMR_SAMPLES.copy()

class GermlineAndTumor:
	def __init__(self, mapping: pd.DataFrame):
		self.ALL_PATIENTS = mapping['patient'].dropna().to_list()
		self.GRM_VS_TMR_PATIENTS = mapping.query("~(sample_tmr.isnull() | sample_grm.isnull())")['patient'].to_list()
		self.ONLY_TMR_PATIENTS =  mapping.query("sample_grm.isnull()")['patient'].to_list()

class NGSSetup(Germline, Tumor, GermlineAndTumor):
	def __init__(self):

		self.GRM=config['grm_dir'] != ''
		self.TMR=config['tmr_dir'] != ''
		self.PAIR=config['reads_type'] == 'pair'

		self.GRM_SAMPLES = []
		self.ALL_PATIENTS = []
		self.GRM_VS_TMR_PATIENTS = []
		self.ONLY_TMR_PATIENTS = []

		ngs = self.create_params()
		self.data = ngs.data
		self.mapping = self.check_mapping(ngs)

		self.LANES = ngs.data['lane'].dropna().unique().tolist()
		if self.GRM:
			Germline.__init__(self, self.data)
		if self.TMR:
			Tumor.__init__(self, self.data)
		if self.GRM & self.TMR:
			GermlineAndTumor.__init__(self, self.mapping)

	def create_params(self) -> NGS:
		if self.GRM and self.TMR and self.PAIR:
			ngs = NGSJointPaired(config['grm_dir'], config['tmr_dir'])
		elif (self.GRM and self.TMR == False) and self.PAIR:
			ngs = NGSIndividualPaired(config['grm_dir'])
		elif (self.GRM == False and self.TMR) and self.PAIR:
			ngs = NGSIndividualPaired(config['tmr_dir'])
		elif self.GRM and self.TMR and self.PAIR == False:
			ngs = NGSJointSingle(config['grm_dir'], config['tmr_dir'])
		elif (self.GRM == False or self.TMR == False) and self.PAIR == False:
			ngs = NGSIndividualSingle(config['grm_dir'])
		return ngs
	
	def check_mapping(self, ngs) -> pd.DataFrame:
		if self.GRM and self.TMR == False:
			return Mapping.map_names(ngs.data, '_grm')
		elif self.GRM == False and self.TMR:
			return Mapping.map_names(ngs.data, '_tmr')
		else:
			return ngs.mapping

