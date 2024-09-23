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

from typing import List, Tuple


pd.set_option('display.max_columns', None)

#####################################
# family of sequencing data classes #
#####################################
class SeqDir:
	def __init__(self, seq_dir: str):
		self.seq_dir = os.path.abspath(seq_dir)
		self.seqs = self.get_seq_list(self.seq_dir)

	def get_seq_list(self, directory_path: str) -> List[str]:
		"""Return a list of sequence files in the given directory."""
		sequence_extensions = ('*.fastq', '*.fastq.gz', '*.fq', '*.fq.gz')
		sequence_files = (glob(os.path.join(directory_path, extension)) for extension in sequence_extensions)
		return list(chain.from_iterable(sequence_files))

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

	def convert_to_dataframe(self, fastq_list: List[str]) -> pd.DataFrame:
		"""Convert a list of fastq files to a pandas DataFrame."""
		df = pd.DataFrame({'fastq': fastq_list})
		df['base_fastq'] = df['fastq'].apply(os.path.basename)
		df['dir_fastq'] = df['fastq'].apply(os.path.dirname) + '/'
		df['link_dir'] = df['dir_fastq'] + 'links/'
		return df

	def trim_sample_name(self, df: pd.DataFrame) -> pd.DataFrame:
		df['full_name'] = df['patient']
		"""Trim the sample names to a common prefix for all samples."""
		if len(df["patient"].unique()) > 1:
			prefix_lengths = []
			for prefix_length in range(1, len(df["patient"].iloc[0]) + 1):
				prefixes = df["patient"].str[:prefix_length]
				if len(prefixes.unique()) == len(df["patient"].unique()):
					prefix_lengths.append(prefix_length)
					break
			df["patient"] = df["patient"].str.replace("[\-\.\|]", "_", regex=True)
			df["patient"] = df["patient"].str[:prefix_lengths[-1]]
		return df

	def write_links(self, df: pd.DataFrame) -> pd.DataFrame:
		"""Write links from input fastq files to a new location."""
		for _, row in df.iterrows():
			source_file = row['fastq']
			link_directory = row['link_dir']
			link_filename = row['link_name']

			# Create the link directory if it doesn't exist
			os.makedirs(link_directory, exist_ok=True)
			if not os.path.exists(link_filename):
				# Create the link
				os.symlink(source_file, link_filename)

		# Drop columns that are no longer needed
		df = df.drop(['fastq', 'base_fastq', 'dir_fastq', 'link_dir'], axis=1)

		# Rename the link_name column to fastq
		df = df.rename({'link_name': 'fastq'}, axis=1)

		return df

class DataProcessorSingle(DataProcessor):
	def __init__(self, path: str) -> pd.DataFrame:
		"""
		Initialize DataProcessorSingle with a path to a fastq file.

		Args:
			path (str): Path to the fastq file.

		Returns:
			pd.DataFrame: A pandas DataFrame containing the sample information.
		"""
		fastqs = self.extract_list(path)
		df = super().convert_to_dataframe(fastqs)
		df = self.extract_params(df)
		df = super().trim_sample_name(df)
		df = self.create_link(df)
		self.df = super().write_links(df)

	def extract_list(self, path: str) -> List[str]:
		"""Extract a list of fastq files from a given directory."""
		seq_dir = SeqDir(path)
		return seq_dir.seqs

	def extract_params(self, df: pd.DataFrame) -> pd.DataFrame:
		"""Extract patient and lane information from the base_fastq column."""
		extracted = (
			df['base_fastq']
			.str.extractall(r'(?P<patient>.+?)[_\-.](?P<lane>[lL][\d]*)[_\-\.]?.*?')
			.droplevel(1)
		)

		extracted['lane'] = extracted['lane'].fillna('L001')
		return pd.concat([df, extracted], axis=1)

	def create_link(self, df: pd.DataFrame) -> pd.DataFrame:
		"""Create a 'link_name' column with the path to the link for each fastq file."""
		df["link_name"] = df["link_dir"] + df["patient"] + "_" + df["lane"] + ".fastq"
		return df

class DataProcessorPaired(DataProcessorSingle):
	def __init__(self, path: str) -> pd.DataFrame:
		"""
		Initialize DataProcessorPaired with a path to a directory containing fastq files.

		Args:
			path (str): Path to the directory.

		Returns:
			pd.DataFrame: A pandas DataFrame containing sample information.
		"""
		fastqs = self.extract_list(path)
		df = super().convert_to_dataframe(fastqs)
		df = self.extract_params(df)
		df = super().trim_sample_name(df)
		df = self.create_link(df)
		self.df = super().write_links(df)

	def extract_params(self, df: pd.DataFrame) -> pd.DataFrame:
		"""
		Extract patient and lane information from the base_fastq column.
		"""
		extracted = (
			df['base_fastq']
			.str.extractall(r'(?P<patient>.+?)[_\-\.](?P<lane>[lL][0-9M]*)?[_\-\.]?(?P<reads_orientation>[Rr][12])[_\-\.]?.*?')
			.droplevel(1)
		)

		extracted['lane'] = extracted['lane'].fillna('L001')
		return pd.concat([df, extracted], axis=1)
	
	def create_link(self, df: pd.DataFrame) -> pd.DataFrame:
		"""Create a 'link_name' column with the path to the link for each fastq file."""
		df['link_name'] = df['link_dir'] + df['patient'] + '_' + df['lane'] + '_' + df['reads_orientation'] + '.fastq'
		return df

class Mapping:
	@staticmethod
	def map_names(df, suffix) -> pd.DataFrame:
		"""Map a dataframe to a new dataframe with the sample names	appended with the given suffix."""
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

	def estimate_distances(self, grm: DataProcessor, tmr: DataProcessor) -> List[Tuple[float, float]]:
		"""Estimates the distances between patients in two datasets, returning a list of tuples containing patient pairs and their sequence similarity ratios."""
		distances = []
		germ_patients = grm.df['patient'].tolist()
		tmr_patients = tmr.df['patient'].tolist()
		for gp in germ_patients:
			for tp in tmr_patients:
				distances.append((gp, tp, SequenceMatcher(None, gp, tp).ratio()))
		return distances

	def map_names(self, mapping: list) -> pd.DataFrame:
		"""
		Maps the distances between patients in two datasets, returning a dataframe with the matched and unmatched samples.
		
		Parameters:
		mapping (list): A list of tuples containing patient pairs and their sequence similarity ratios.
		
		Returns:
		pd.DataFrame: A dataframe with the matched and unmatched samples, including their sample names and patient IDs.
		"""
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
