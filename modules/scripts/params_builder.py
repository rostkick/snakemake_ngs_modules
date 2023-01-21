import re
import os
import pandas as pd


class ParamsTable:
	pass


class ReadWorker:
	@staticmethod
	def sort_fastq_by_pattern(fastq_list: list, pattern: str) -> list:
		return sorted(fastq_list, key = lambda x: int(re.search(pattern, x).group(1)))

	@staticmethod
	def extract_forward(fastq_list: list) -> list:
		return [i for i in fastq_list if re.match('.*_R1_.*', i) is not None]

	@staticmethod
	def extract_reverse(fastq_list: list) -> list:
		return [i for i in fastq_list if re.match('.*_R2_.*', i) is not None]

class ReadDFWorker:
	@staticmethod
	def get_reads_table(forward_reads: list, reverse_reads: list) -> pd.core.frame.DataFrame:
		return 	pd.DataFrame({'forward_reads': forward_reads, 'reverse_reads': reverse_reads})

	@staticmethod
	def extract_sample_name(fastq_df: pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
		fastq_df.loc[:, 'samples'] = \
			fastq_df.loc[:, 'forward_reads'] \
				.apply(lambda x: os.path.basename(x)) \
				.str.extract('(?P<samples>.+)_R1.*')

class FastqList:
	def __init__(self, fastq_list: list) -> None:
		self.fastq_list = fastq_list

class FastqListSortedByPattern(FastqList):
	def __init__(self, fastq_list: list, sorting_pattern: str) -> None:
		super().__init__(fastq_list)
		self.sorting_pattern = sorting_pattern
		self.fastq_list = ReadWorker.sort_fastq_by_pattern(self.fastq_list, self.sorting_pattern)

class PairedReads(FastqListSortedByPattern):
	def __init__(self, fastq_list: list, sorting_pattern: str) -> None:
		super().__init__(fastq_list, sorting_pattern)
		self.forward_reads = ReadWorker.extract_forward(self.fastq_list)
		self.reverse_reads = ReadWorker.extract_reverse(self.fastq_list)

class PairedReadsDF(PairedReads):
	def __init__(self, fastq_list: list, sorting_pattern: str) -> None:
		super().__init__(fastq_list, sorting_pattern)
		self.fastq_df = ReadDFWorker.get_reads_table(self.forward_reads, self.reverse_reads)
