import os
import re
import subprocess as sp
from glob import glob
from abc import ABCMeta, abstractmethod
import pandas as pd
from snakemake.workflow import config


#####################################
# family of sequencing data classes #
#####################################
class SeqDir:
	def __init__(self, fastq_pattern: str) -> list:
		self.seqs = self.get_seq_list(fastq_pattern)

	def get_seq_list(self, fastq_pattern: str) -> list:
		paths = glob(fastq_pattern)
		dname = os.path.dirname(os.path.abspath(paths[0]))
		paths = [os.path.join(dname, os.path.basename(path)) for path in paths] 
		return paths

class DataProcessor:
	def __init__(self, path: str) -> pd.DataFrame:
		fastq = self.extract_list(path)
		df = self.convert_to_dataframe(fastq)
		df = self.extract_params(df)
		df = self.trim_sample_name(df)
		df = self.create_link(df)
		self.df = self.write_link(df)
	
	def extract_list(self, path: str) -> list:
		seq = SeqDir(path)
		return seq.seqs

	def convert_to_dataframe(self, fastq: list) -> pd.DataFrame:
		df = pd.DataFrame({'fastq': fastq})
		df['base_fastq'] = df['fastq'].apply(lambda x: os.path.basename(x))
		df['dir_fastq'] = df['fastq'].apply(lambda x: os.path.dirname(x)+'/')
		df['link_dir'] = df['dir_fastq']+ "links/"
		return df

	def extract_params(self, df: pd.DataFrame) -> pd.DataFrame:
		df_extracted = df['base_fastq'].str.extractall(r'(?P<sample>.+?)[_\-\.](?P<lane>[lL][0-9M]*)?[_\-\.]?(?P<reads_orientation>[Rr][12])[_\-\.]?.*?').droplevel(1)
		df_extracted['lane'] = df_extracted['lane'].fillna('L001')
		df = pd.concat([df, df_extracted], axis=1)
		return df

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
	
	def create_link(self, df: pd.DataFrame) -> pd.DataFrame:
		df['link_name'] = df['link_dir']+df['sample'] + '_' + df['lane'] + '_' + df['reads_orientation'] + '.fastq'
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


# class DataProcessor(DataProcessor):



##################################
# family of wide-format ngs data #
##################################
class NGS:
	def __init__(self, path: str):
		self.data = self.get_data_processor(path)

	def get_data_processor(self, path: str) -> DataProcessor:
		return DataProcessor(path).df

class NGSVariables:
	def __init__(self):
		self.SAMPLES = self.data.loc[:, 'sample'].dropna().unique().tolist()

class NGSSetup(NGSVariables):
	def __init__(self):

		if 	config['reads'] == '':
			raise Exception('Указан некорректный путь!')

		self.SAMPLES = []
		ngs = NGS(config['reads'])
		self.data = ngs.data
		self.LANES = ngs.data['lane'].dropna().unique().tolist()

		NGSVariables.__init__(self)
