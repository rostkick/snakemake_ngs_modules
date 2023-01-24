import os
from itertools import chain
import pandas as pd
import numpy as np
from shutil import copyfile
import snakemake

from modules.scripts.functions import *


configfile: "configure.yml"

RUN = config['run_id']
MODE = config['mode']

data = pd.read_table(config['params_table'])

GERMLINE_SAMPLES = data.loc[:, 'germline_samples'].to_list()
if ('somatic' in MODE) | (MODE == 'all'):
	long_data = pd.DataFrame(
		{
			"samples": np.concatenate([data.loc[:, 'tumor_samples'].values,
									data.loc[:, 'germline_samples'].values]),
			"fastq_forward": np.concatenate([data.loc[:, "tumor_forward"].values,
											data.loc[:, "germline_forward"].values]),
			"fastq_reverse": np.concatenate([data.loc[:, "tumor_reverse"].values,
											data.loc[:, "germline_reverse"].values])
		}
	)
	PATIENTS = data.loc[:, 'patient'].to_list()
	TUMOR_SAMPLES = data.loc[:, 'tumor_samples'].to_list()	
	SAMPLES = long_data.loc[:, "samples"].tolist()
else:
	SAMPLES = GERMLINE_SAMPLES
	PATIENTS = ''
	data.rename(columns={'germline_samples': 'samples',
						'germline_forward': 'fastq_forward',
						'germline_reverse': 'fastq_reverse'}, inplace=True)
	long_data = data.copy()

print(SAMPLES)
wildcard_constraints:
	patient="|".join(PATIENTS),
	sample="|".join(SAMPLES)

rule all:
	input: get_inputs(RUN, MODE, PATIENTS)

include: config["snakemake_modules"] + "aligning.smk"
include: config["snakemake_modules"] + "preprocessing.smk"
# include: config["snakemake_modules"] + "germline_calling.gatk.smk"
include: config["snakemake_modules"] + "germline_calling.deepvariant.smk"
include: config["snakemake_modules"] + "somatic_calling.smk"
include: config["snakemake_modules"] + "sv_calling.smk"
include: config["snakemake_modules"] + "liftover.smk"
include: config["snakemake_modules"] + "annotation.smk"
