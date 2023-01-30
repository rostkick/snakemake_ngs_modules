import os
import re
from glob import glob
from abc import ABCMeta, abstractmethod
import pandas as pd
import numpy as np
from itertools import product
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
		return glob(dir_path+'/*.fastq.gz')

class DataProcessor(metaclass=ABCMeta):
	@abstractmethod
	def __init__(self, path:str) -> pd.DataFrame:
		pass
	@abstractmethod
	def extract_list(self, path: str) -> list:
		pass
	@abstractmethod
	def get_sample_name(self, fastq: list) -> pd.Series:
		pass
	@abstractmethod
	def convert_to_df(self, sample: pd.Series, fastq: list) -> pd.Series:
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
		sample = self.get_sample_name(fastq)
		df = self.convert_to_df(sample, fastq)
		self.df = self.trim_sample_name(df)

	def extract_list(self, path: str) -> list:
		seq = SeqDir(path)
		return seq.seqs

	def get_sample_name(self, fastq) -> pd.Series:
		sample = [os.path.basename(fq) for fq in fastq]
		return pd.Series(sample)

	def convert_to_df(self, sample: pd.Series, fastq: list) -> pd.DataFrame:
		df = pd.DataFrame()
		df['sample'] = sample
		df['fastq'] = fastq
		return df

class DataProcessorPaired(DataProcessorSingle):
	def __init__(self, path: str) -> pd.DataFrame:
		fastq = super().extract_list(path)
		fastq_f, fastq_r = self.extract_paired_reads(fastq)
		sample = super().get_sample_name(fastq_f)
		df = self.convert_to_df(sample, fastq_f, fastq_r)
		self.df = super().trim_sample_name(df)

	def extract_paired_reads(self, fastq: list) -> list:
		fastq_f = [fq for fq in fastq if re.match('.*_R1_.*', fq) is not None]
		fastq_r = [fq for fq in fastq if re.match('.*_R2_.*', fq) is not None]
		return fastq_f, fastq_r

	def convert_to_df(self, sample: pd.Series, fastq_f: list, fastq_r: list) -> pd.DataFrame:
		df = pd.DataFrame()
		df['sample'] = sample
		df['fastq_f'] = fastq_f
		df['fastq_r'] = fastq_r
		return df

class MappedGermTumor:
	DISTANCE_THRESHOLD = 0.5
	def __init__(self, grm: DataProcessor, tmr: DataProcessor):
		distances = self.estimate_distances(grm, tmr)
		self.df = self.map_names(distances)


	def estimate_distances(self, grm: DataProcessor, tmr: DataProcessor) -> list:
		distances = []
		for g, t in product(grm.df['sample'].tolist(), tmr.df['sample'].tolist()):
			distances.append((g, t, SequenceMatcher(None, g, t).ratio()))
		return distances

	def map_names(self, distances: list) -> pd.DataFrame:
		df = pd.DataFrame(distances)
		df.columns = ['grm_samples', 'tmr_samples', 'distance']

		df = df.sort_values('distance', ascending=False)
		df_matched = df.loc[df['distance']>=self.DISTANCE_THRESHOLD, :]
		
		df_grm_unmatched = pd.DataFrame()
		df_grm_unmatched['grm_samples'] = df.loc[df['distance']<self.DISTANCE_THRESHOLD, 'grm_samples'].drop_duplicates().dropna()
		df_grm_unmatched['tmr_samples'] = np.nan
		df_grm_unmatched['distance'] = np.nan
		
		df_tmr_unmatched = pd.DataFrame()
		df_tmr_unmatched['grm_samples'] = np.nan
		df_tmr_unmatched['tmr_samples'] = df.loc[df['distance']<self.DISTANCE_THRESHOLD, 'tmr_samples'].drop_duplicates().dropna()
		df_tmr_unmatched['distance'] = np.nan
		
		df = pd.concat([df_matched, df_grm_unmatched, df_tmr_unmatched])
		is_duplicated = df.apply(pd.Series.duplicated, axis=0)
		df = df.where(~is_duplicated, np.nan)
		df = df.drop('distance', axis=1)
		df = df.dropna(axis=0, how='all')
		return df

##################################
# family of wide-format ngs data #
##################################
class NGSWide(metaclass=ABCMeta):
	@abstractmethod
	def __init__(self):
		pass
	@abstractmethod
	def get_data_processor(self, path) -> DataProcessor:
		pass
	@abstractmethod
	def get_wide_table(self) -> pd.DataFrame:
		pass
	@abstractmethod
	def add_patient_ids(self) -> pd.DataFrame:
		pass

class NGSWideJointPaired(NGSWide):
	def __init__(self, grm_path: str, tmr_path: str):
		grm, tmr = self.get_data_processor(grm_path, tmr_path)
		mapping = self.get_mapping(grm, tmr)
		wide_df = self.get_wide_table(grm, tmr, mapping)
		wide_df = self.add_patient_ids(wide_df)
		self.wide_df = self.add_suffixes(wide_df)

	def get_data_processor(self, grm_path: str, tmr_path: str) -> DataProcessorPaired:
		return DataProcessorPaired(grm_path), DataProcessorPaired(tmr_path)

	def get_mapping(self, grm: DataProcessor, tmr: DataProcessor) -> MappedGermTumor:
		return MappedGermTumor(grm, tmr)
	
	def get_wide_table(self, grm: DataProcessor, tmr: DataProcessor, mapping: MappedGermTumor) -> pd.DataFrame:
		wide_df = pd.merge(mapping.df, grm.df, how='outer', left_on='grm_samples', right_on='sample')
		wide_df = wide_df.drop('sample', axis=1)
		wide_df.columns = ['grm_' + str(col) if 'fastq' in col else col for col in wide_df.columns]
		wide_df = pd.merge(wide_df, tmr.df, how='outer', left_on='tmr_samples', right_on='sample')
		wide_df = wide_df.drop('sample', axis=1)
		wide_df.columns = ['tmr_' + str(col) if col.startswith('fastq') else col for col in wide_df.columns]
		return wide_df

	def add_patient_ids(self, wide_df: pd.DataFrame) -> pd.DataFrame:
		if wide_df['tmr_samples'].isna().sum() < wide_df['grm_samples'].isna().sum():
			wide_df['patients'] = wide_df['tmr_samples']
		else:
			wide_df['patients'] = wide_df['grm_samples']
		return wide_df
		
	def add_suffixes(self, wide_df: pd.DataFrame) -> pd.DataFrame:
		wide_df['grm_samples'] = wide_df['grm_samples'] + '_grm'
		wide_df['tmr_samples'] = wide_df['tmr_samples'] + '_tmr'
		return wide_df

class NGSWideIndividualPaired(NGSWide):
	def __init__(self, path):
		ngs = self.get_data_processor(path)
		wide_df = self.get_wide_table(ngs)
		self.wide_df = self.add_patient_ids(wide_df)

	def get_data_processor(self, path: str) -> DataProcessorPaired:
		return DataProcessorPaired(path)

	def get_wide_table(self, ngs: DataProcessorPaired) -> pd.DataFrame:
		return ngs.df

	def add_patient_ids(self, wide_df: pd.DataFrame) -> pd.DataFrame:
		wide_df['patients'] = wide_df['sample']
		return wide_df

class NGSWideJointSingle(NGSWideJointPaired):
	def __init__(self, grm_path, tmr_path):
		grm, tmr = self.get_data_processor(grm_path, tmr_path)
		mapping = super().get_mapping(grm, tmr)
		wide_df = super().get_wide_table(grm, tmr, mapping)
		wide_df = super().add_patient_ids(wide_df)
		self.wide_df = super().add_suffixes(wide_df)

	def get_data_processor(self, grm_path: str, tmr_path: str) -> DataProcessorSingle:
		return DataProcessorSingle(grm_path), DataProcessorSingle(tmr_path)


class NGSWideIndividualSingle(NGSWideIndividualPaired):
	def __init__(self, path: str):
		ngs = self.get_data_processor(path)
		wide_df = super().get_wide_table(ngs)
		self.wide_df = super().add_patient_ids(wide_df)

	def get_data_processor(self, path: str) -> DataProcessorSingle:
		return DataProcessorSingle(path)


##################################
# family of long-format ngs data #
##################################
class NGSLong(metaclass=ABCMeta):
	@abstractmethod
	def get_long_table(self):
		pass

class NGSLongJointPaired(NGSWideJointPaired, NGSLong):
	def __init__(self, grm_path: str, tmr_path: str):
		super().__init__(grm_path, tmr_path)
		self.long_df = self.get_long_table(self.wide_df)
	def get_long_table(self, wide_df) -> pd.DataFrame:
		long_df = pd.DataFrame(
			{
			"samples": np.concatenate([wide_df.loc[:, 'tmr_samples'].values,
									wide_df.loc[:, 'grm_samples'].values]),
			"fastq_forward": np.concatenate([wide_df.loc[:, "tmr_fastq_f"].values,
											wide_df.loc[:, "grm_fastq_f"].values]),
			"fastq_reverse": np.concatenate([wide_df.loc[:, "tmr_fastq_r"].values,
											wide_df.loc[:, "grm_fastq_r"].values])
			}
									)
		long_df = long_df.dropna(axis=0, how='all')
		return long_df

class NGSLongIndividualPaired(NGSWideIndividualPaired, NGSLong):
	def __init__(self, path: str):
		super().__init__(path)
		self.long_df = self.get_long_table(self.wide_df)
	def get_long_table(self, wide_df) -> pd.DataFrame:
		long_df = wide_df.drop('patients', axis=1)
		long_df = long_df.rename({'sample':'samples', 'fastq_f':'fastq_forward', 'fastq_r':'fastq_reverse'}, axis=1)
		long_df = long_df.dropna(axis=0, how='all')
		return long_df


class NGSLongJointSingle(NGSWideJointSingle, NGSLong):
	def __init__(self, grm_path: str, tmr_path: str):
		super().__init__(grm_path, tmr_path)
		self.long_df = self.get_long_table(self.wide_df)
	def get_long_table(self, wide_df) -> pd.DataFrame:
		long_df = pd.DataFrame(
			{
			"samples": np.concatenate([wide_df.loc[:, 'tmr_samples'].values,
									wide_df.loc[:, 'grm_samples'].values]),
			"fastq": np.concatenate([wide_df.loc[:, "tmr_fastq"].values,
									wide_df.loc[:, "grm_fastq"].values])})
		long_df = long_df.dropna(axis=0, how='all')
		return long_df

class NGSLongIndividualSingle(NGSWideIndividualSingle, NGSLong):
	def __init__(self, path: str):
		super().__init__(path)
		self.long_df = self.get_long_table(self.wide_df)
	def get_long_table(self, wide_df) -> pd.DataFrame:
		long_df = wide_df.drop('patients', axis=1)
		long_df = long_df.rename({'sample':'samples', 'fastq':'fastq'}, axis=1)
		long_df = long_df.dropna(axis=0, how='all')
		return long_df



####################################
# family of high-level ngs classes #
####################################
class NGS(metaclass=ABCMeta):
	@abstractmethod
	def __init__(self):
		pass

class NGSJointPaired(NGSLongJointPaired, NGS):
	def __init__(self, grm_path: str, tmr_path: str):
		super().__init__(grm_path, tmr_path)

class NGSIndividualPaired(NGSLongIndividualPaired, NGS):
	def __init__(self, path: str):
		super().__init__(path)

class NGSJointSingle(NGSLongJointSingle, NGS):
	def __init__(self, grm_path: str, tmr_path: str):
		super().__init__(grm_path, tmr_path)

class NGSIndividualSingle(NGSLongIndividualSingle, NGS):
	def __init__(self, path: str):
		super().__init__(path)

class NGSSetup:
	def __init__(self):
		self.GRM=config['grm_dir'] != ''
		self.TMR=config['tmr_dir'] != ''
		self.PAIR=config['reads_type'] == 'pair'

	def create_params(self) -> NGS:
		if self.GRM and self.TMR and self.PAIR:
			data = NGSJointPaired(config['grm_dir'], config['tmr_dir'])
		elif (self.GRM == False or self.TMR == False) and self.PAIR:
			data = NGSIndividualPaired(config['grm_dir'])
		elif self.GRM and self.TMR and self.PAIR == False:
			data = NGSJointSingle(config['grm_dir'], config['tmr_dir'])
		elif (self.GRM == False or self.TMR == False) and self.PAIR == False:
			data = NGSIndividualSingle(config['grm_dir'])
		return data
