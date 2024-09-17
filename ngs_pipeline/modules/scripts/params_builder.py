import os
import re
import subprocess as sp
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
	"""
	The DataProcessor class is an abstract base class, meaning it cannot be instantiated directly.
	It is meant to be subclassed and its abstract methods implemented.
	"""
	@abstractmethod
	def __init__(self, path:str) -> pd.DataFrame:
		pass
	@abstractmethod
	def extract_list(self, path: str) -> list:
		pass
	@abstractmethod
	def extract_params(self, sample: pd.Series, fastq: list) -> pd.Series:
		pass
	@abstractmethod
	def convert_to_dataframe(self, sample: pd.Series, fastq: list) -> pd.Series:
		pass
	@abstractmethod
	def create_link(self, df: pd.Series) -> pd.DataFrame:
		pass

	def convert_to_dataframe(self, fastq: list) -> pd.DataFrame:
		df = pd.DataFrame({'fastq': fastq})
		df['base_fastq'] = df['fastq'].apply(lambda x: os.path.basename(x))
		df['dir_fastq'] = df['fastq'].apply(lambda x: os.path.dirname(x)+'/')
		df['link_dir'] = df['dir_fastq']+ "links/"
		return df

	def trim_sample_name(self, df=None) -> pd.DataFrame:
		if len(df['patient']) > 1:
			series_list = df['patient'].str.split(r'_|\.|-')
			i=0
			condition = '_'.join(series_list[0][:i]) == '_'.join(series_list[len(series_list)-1][:i])
			while condition:
				substr_set = set(series_list.str[:i].str.join('_').tolist())
				if len(substr_set) == len(set(df['patient'].tolist())):
					break
				i+=1
			df['patient'] = series_list.str[:i].str.join('_')
		else:
			df['patient'] = df['patient'].str.extract('(.+?)[_\-\.\|].*')[0]
		return df

	def write_link(self, df: pd.DataFrame) -> pd.DataFrame:
		fastq = df['fastq'].tolist()
		link_dir = df['link_dir'].tolist()
		links = df['link_name'].tolist()
		link_dir = df['link_dir']
	
		for fq, link_d, link in zip(fastq, link_dir, links):
			sp.run(f'mkdir -p {link_d}', shell=True)
			sp.run(f'ln -fs {fq} {link}', shell=True)

		df = df.drop(['fastq', 'base_fastq', 'dir_fastq', 'link_dir'], axis=1)
		df = df.rename({'link_name': 'fastq'}, axis=1)
		return df

class DataProcessorSingle(DataProcessor):
	def __init__(self, path:str) -> pd.DataFrame:
		fastq = self.extract_list(path)
		df = super().convert_to_dataframe(fastq)
		df = self.extract_params(df)
		df = super().trim_sample_name(df)
		df = self.create_link(df)
		self.df = super().write_link(df)

	def extract_list(self, path: str) -> list:
		seq = SeqDir(path)
		return seq.seqs

	def extract_params(self, df: pd.DataFrame) -> pd.DataFrame:
		df_extracted = df['base_fastq'].str.extractall(r'(?P<patient>.*)[_\-\.](?P<lane>[lL][\d]*)[_\-\.]?.*?').droplevel(1)
		df_extracted['lane'] = df_extracted['lane'].fillna('L001')
		df = pd.concat([df, df_extracted], axis=1)
		return df

	def create_link(self, df: pd.DataFrame) -> pd.DataFrame:
		df['link_name'] = df['link_dir']+df['patient'] + '_' + df['lane'] + '_' + '.fastq'
		return df

class DataProcessorPaired(DataProcessorSingle):
	def __init__(self, path: str) -> pd.DataFrame:
		fastq = super().extract_list(path)
		df = super().convert_to_dataframe(fastq)
		df = self.extract_params(df)
		df = super().trim_sample_name(df)
		df = self.create_link(df)
		self.df = super().write_link(df)

	def extract_params(self, df: pd.DataFrame) -> pd.DataFrame:
		df_extracted = df['base_fastq'].str.extractall(r'(?P<patient>.+?)[_\-\.](?P<lane>[lL][0-9M]*)?[_\-\.]?(?P<reads_orientation>[Rr][12])[_\-\.]?.*?').droplevel(1)
		df_extracted['lane'] = df_extracted['lane'].fillna('L001')
		df = pd.concat([df, df_extracted], axis=1)
		return df
	
	def create_link(self, df: pd.DataFrame) -> pd.DataFrame:
		df['link_name'] = df['link_dir']+df['patient'] + '_' + df['lane'] + '_' + df['reads_orientation'] + '.fastq'
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
	def __init__(self, grm: DataProcessor, tmr: DataProcessor):
		self.DISTANCE_THRESHOLD = 0.6
		self.distances = self.estimate_distances(grm, tmr)
		self.mapping = self.map_names(self.distances)

	def estimate_distances(self, grm: DataProcessor, tmr: DataProcessor) -> list:
		distances = []
		for g, t in product(grm.df['patient'].tolist(), tmr.df['patient'].tolist()):
			distances.append((g, t, SequenceMatcher(None, g, t).ratio()))
		return distances

	def map_names(self, mapping: list) -> pd.DataFrame:
		df = pd.DataFrame(mapping)
		df.columns = ['sample_grm', 'sample_tmr', 'distance']

		df = df.sort_values('distance', ascending=False)
		max_dist_current = df['distance'].max()
		max_dist = max_dist_current if max_dist_current > self.DISTANCE_THRESHOLD else self.DISTANCE_THRESHOLD
		df_matched = df.loc[df['distance']==max_dist, :].drop_duplicates()
		df_grm_unmatched = pd.DataFrame()
		df_grm_unmatched.loc[:, 'sample_grm'] = df.loc[df['distance']<max_dist, 'sample_grm'].drop_duplicates().dropna()
		grm_mask = ~df_grm_unmatched.loc[:, 'sample_grm'].isin(df_matched.loc[:, 'sample_grm'])
		df_grm_unmatched = df_grm_unmatched.loc[grm_mask, :].dropna()

		if df_grm_unmatched.size != 0:
			df_grm_unmatched.loc[:, 'sample_tmr'] = np.nan
			df_grm_unmatched.loc[:, 'distance'] = np.nan

		df_tmr_unmatched = pd.DataFrame()
		df_tmr_unmatched.loc[:, 'sample_tmr'] = df.loc[df['distance']<=max_dist, 'sample_tmr'].drop_duplicates().dropna()
		tmr_mask = ~df_tmr_unmatched.loc[:, 'sample_tmr'].isin(df_matched.loc[:, 'sample_tmr'])
		df_tmr_unmatched = df_tmr_unmatched.loc[tmr_mask, :].dropna()

		if df_tmr_unmatched.size != 0:
			df_tmr_unmatched.loc[:, 'distance'] = np.nan
			df_tmr_unmatched['sample_grm'] = np.nan

		df = pd.concat([df_matched, df_grm_unmatched, df_tmr_unmatched])
		is_duplicated = df.apply(pd.Series.duplicated, axis=0)
		df = df.where(~is_duplicated, np.nan)
		df = df.drop('distance', axis=1)
		mapping = df.dropna(axis=0, how='all').reset_index()
		mapping.loc[:, 'sample_grm'] = mapping['sample_grm'] + '_grm'
		mapping.loc[:, 'sample_tmr'] = mapping['sample_tmr'] + '_tmr'
		mapping.loc[:, 'patient'] = mapping['sample_tmr'].str[:-4]
		mapping['patient'].fillna(mapping['sample_grm'].str[:-4], inplace=True)
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
	def __init__(self, path: str, suffix: str):
		self.mapping = pd.DataFrame()
		self.data = self.get_data_processor(path)
		self.data['sample'] = self.data['patient'] + suffix

	def get_data_processor(self, path: str) -> DataProcessorPaired:
		return DataProcessorPaired(path).df

class NGSIndividualSingle(NGSIndividualPaired):
	def get_data_processor(self, path: str, suffix: str) -> DataProcessorSingle:
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
		mapping['sample'] = mapping['sample_grm'].str[:-4]
		wide_df = pd.merge(mapping, grm.df, how='outer', left_on='patient', right_on='patient')
		wide_df = pd.merge(wide_df, tmr.df, how='outer', left_on=['patient'],
														right_on=['patient'], suffixes=['_grm', '_tmr'])
		wide_df['patient'] = wide_df.loc[:, 'sample'].str[:-4]
		wide_df = wide_df.drop('sample', axis=1)
		wide_df = wide_df.drop('index', axis=1)
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
		long_df['patient'] = long_df['sample'].str[:-4]
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
		long_df['patient'] = long_df['sample'].str[:-4]
		return long_df

####################################
# family of high-level ngs classes #
####################################
class Germline:
	def __init__(self):
		self.SAMPLES = self.data.loc[:, 'sample'].dropna().unique().tolist()
		self.GRM_SAMPLES = [i for i in self.SAMPLES if '_grm' in i]
		self.GRM_SAMPLES = self.SAMPLES if len(self.GRM_SAMPLES) == 0 else self.GRM_SAMPLES

class Tumor:
	def __init__(self):
		self.SAMPLES = self.data.loc[:, 'sample'].dropna().unique().tolist()
		self.TMR_SAMPLES = [i for i in self.SAMPLES if '_tmr' in i]
		self.TMR_SAMPLES = self.SAMPLES if len(self.TMR_SAMPLES) == 0 else self.TMR_SAMPLES
		self.TMR_PATIENTS = [tmr[:-4] for tmr in self.TMR_SAMPLES] # mapping['patient'].dropna().unique().tolist()
		self.ONLY_TMR_PATIENTS = self.TMR_PATIENTS.copy()

class GermlineAndTumor:
	def __init__(self):
		self.SAMPLES = self.mapping['sample_grm'].dropna().unique().tolist() + self.mapping['sample_tmr'].dropna().unique().tolist()
		self.TMR_PATIENTS = self.mapping['patient'].dropna().to_list()
		self.TMR_PATIENTS = self.mapping.loc[~self.mapping['sample_tmr'].isna(), 'patient'].to_list()
		self.GRM_VS_TMR_PATIENTS = self.mapping.query("~(sample_tmr.isnull() | sample_grm.isnull())")['patient'].to_list()
		self.ONLY_TMR_PATIENTS =  self.mapping.query("sample_grm.isnull()")['patient'].to_list()

class NGSSetup(Germline):#, Tumor, GermlineAndTumor):
	def __init__(self):

		self.GRM=config['grm_dir'] != ''
		self.TMR=config['tmr_dir'] != ''

		self.check_files()

		self.PAIR=config['reads_type'] == 'pair'

		self.GRM_SAMPLES = []
		self.TMR_SAMPLES = []
		self.TMR_PATIENTS = []
		self.GRM_VS_TMR_PATIENTS = []
		self.ONLY_TMR_PATIENTS = []

		ngs = self.create_params()
		self.data = ngs.data
		self.mapping = self.check_mapping(ngs)

		self.LANES = ngs.data['lane'].dropna().unique().tolist()
		if self.GRM:
			Germline.__init__(self)
		if self.TMR:
			Tumor.__init__(self)
		if self.GRM & self.TMR:
			GermlineAndTumor.__init__(self)
	
	def check_files(self):
		grm_dir = config['grm_dir']
		if self.GRM:
			files = list(chain(*[glob(grm_dir + '/' + ext) for ext in ['*.fastq', '*.fastq.gz', '*.fq', '*.fq.gz']]))
			if len(files) == 0:
				print(f"No reads files were found in germline directory {grm_dir}!")
				self.GRM = False

		tmr_dir = config['tmr_dir']
		if self.TMR:
			files = list(chain(*[glob(tmr_dir + '/' + ext) for ext in ['*.fastq', '*.fastq.gz', '*.fq', '*.fq.gz']]))
			if len(files) == 0:
				print(f"No reads files were found in tumor directory {tmr_dir}!")
				self.TMR = False
		if not (self.GRM | self.TMR):
			raise FileNotFoundError(f"No files matching the pattern were found in both germline and tumor directories.")

	def create_params(self) -> NGS:
		if self.GRM and self.TMR and self.PAIR:
			ngs = NGSJointPaired(config['grm_dir'], config['tmr_dir'])
		elif (self.GRM and self.TMR == False) and self.PAIR:
			ngs = NGSIndividualPaired(config['grm_dir'], '_grm')
		elif (self.GRM == False and self.TMR) and self.PAIR:
			ngs = NGSIndividualPaired(config['tmr_dir'], '_tmr')
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

