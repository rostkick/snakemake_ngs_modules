import os
import argparse
import pandas as pd
import subprocess as sp
from tempfile import NamedTemporaryFile
import re

from modules.scripts.functions import *


def run_snake(args):

	germ = sorted(args.fastq_germline, key = lambda x: int(re.search(args.sample_pattern, x).group(1)))
	germline_forward = [i for i in germ if re.match('.*_R1_.*', i) is not None]
	germline_reverse = [i for i in germ if re.match('.*_R2_.*', i) is not None]

	df = pd.DataFrame({'germline_forward': germline_forward,
					'germline_reverse': germline_reverse})

	df.loc[:, 'germline_samples'] = df.loc[:, 'germline_forward'].apply(lambda x: os.path.basename(x)).str.extract('(?P<germline_samples>.+)_R1.*')
	df.loc[:, 'germline_samples'] = get_min_substr(df.loc[:, 'germline_samples'])

	if ('somatic' in args.mode) | (args.mode == 'all'):

		tumor = sorted(args.fastq_tumor, key = lambda x: int(re.search(args.sample_pattern, x).group(1)))

		tumor_forward = [i for i in tumor if re.match('.*_R1_.*', i) is not None]
		tumor_reverse = [i for i in tumor if re.match('.*_R2_.*', i) is not None]

		if args.patients is None:
			patient = [f'patient_{i}' for i in range(1, len(germline_forward)+1)]
		else:
			patient = args.patients.split(',')

		df['tumor_forward'] = tumor_forward
		df['tumor_reverse'] = tumor_reverse
		df.loc[:, 'tumor_samples'] = df.loc[:, 'tumor_forward'].apply(lambda x: os.path.basename(x)).str.extract('(?P<tumor_samples>.+)_R1.*')
		df['tumor_samples'] = get_min_substr(df.loc[:, 'tumor_samples'])

		df['patient'] = patient
		df.loc[:, 'patient'] = df.loc[:, 'patient'].astype(str).str.replace('\.', '_')
	
	# print(df.sample(n=5).loc[:, df.columns.map(lambda x: x.contains("samples"))])
	print(df.sample(n=5).loc[:, df.columns.str.contains('samples')])
	print('\n')

	df.to_csv(args.output_params, sep='\t', index=False)
	wd_path = os.path.dirname(os.path.realpath(__file__))

	slash = '\\'
	nl = '\n'
	snake_command = \
				(f"snakemake {slash}{nl}"
				f"	--nolock {slash}{nl}"
				f"	-s {wd_path}/snakefile.smk {slash}{nl}"
				f"	--configfile {wd_path}/configure.yml {slash}{nl}"
				f"	--jobs {args.jobs} {slash}{nl}"
				f"	--cores {args.cores} {slash}{nl}"
				f"	--config run_id='{args.run_id}' {slash}{nl}"
				f"			mode='{args.mode}' {slash}{nl}"
				f"			params_table='{args.output_params}' {slash}{nl}"
				f"{'	-p ' + slash + nl if args.printshellcmds else ''}"
				f"{'	-F ' + slash + nl if args.force else ''}"
				f"{'	-q ' + slash + nl if args.quite else ''}"
				f"{'	-n ' + slash + nl if args.dryrun else ''}"
				f"{'	--use-conda ' + slash + nl if args.use_conda else ''}"
				f"{'	--ri ' + slash + nl if args.rerun_incomplete else ''}"
				f"{'	--dag' + slash + nl if args.dag else ''}"
				f"{'	--rulegraph' + slash + nl if args.rulegraph else ''}"
				f"{'	-R ' + args.forcerun + slash + nl if args.forcerun is not None else ''}"
				f"{'	-U ' + args.until + slash + nl if args.until is not None else ''}").rstrip('\\\n')
	if args.cluster & ~(args.dag | args.rulegraph | args.dryrun):
		sbatch_body_begin = \
					(f"sbatch << ENDINPUT{nl}"
					f"#!/bin/sh{nl}"
					f"{nl}"
					f"#SBATCH --job-name={args.run_id}{nl}"
					f"#SBATCH --cpus-per-task={args.jobs}{nl}"
					f"{'#SBATCH --mem=' + args.max_mem + nl if args.max_mem else ''}"
					f"{'#SBATCH --qos=' + args.qos + nl if args.qos else ''}"
					f"{'#SBATCH --time=' + args.time + nl if args.time else ''}"
					f"{'#SBATCH --output=' + args.log + nl if args.log else ''}"
					f"#SBATCH--priority=1{nl}"
					f"{nl}")
		sbatch_body_end = \
					(f"{nl}"
					f"{nl}"
					f"ENDINPUT"
					f"{nl}")
		
		sbatch_task = sbatch_body_begin + snake_command + sbatch_body_end
		
		tmp = NamedTemporaryFile()
		with open(tmp.name, 'w') as w:
			w.write(sbatch_task)
		print(sbatch_task)
		sp.run(f"bash {tmp.name}", shell=True)
	else:
		print(snake_command)
		sp.run(snake_command, shell=True)

	

	


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Process some integers.')

	parser.add_argument('--mode', '-m', type=str, required=True,
						help='')
	parser.add_argument('--run_id', '-ri', type=str, required=True,
						help='')
	parser.add_argument('--patients', '-pa', type=str,
						help='')
	parser.add_argument('--fastq_germline', '-fg', type=str, nargs='+', required=True,
						help='')
	parser.add_argument('--fastq_tumor', '-ft', type=str, nargs='+',
						help='')
	parser.add_argument('--output_params', '-o', type=str, required=True,
						help='')
	parser.add_argument('--sample_pattern', '-sp', type=str, required=True,
						help='')
	parser.add_argument('--jobs', '-j', type=int, required=True,
						help='')
	parser.add_argument('--cores', '-c', type=int, required=True,
						help='')
	parser.add_argument('--force', '-F', action='store_true',
						help='')
	parser.add_argument('--forcerun', '-R', default=None,
						help='')
	parser.add_argument('--until', '-U', default=None,
						help='')
	parser.add_argument('--use_conda', action='store_true',
						help='')
	parser.add_argument('--quite', '-q', action='store_true',
						help='')
	parser.add_argument('--printshellcmds', '-p', action='store_true',
						help='')
	parser.add_argument('--dryrun', '-n', action='store_true',
						help='')
	parser.add_argument('--rerun_incomplete', action='store_true',
						help='')
	parser.add_argument('--dag', action='store_true',
						help='')
	parser.add_argument('--rulegraph', action='store_true',
						help='')
	parser.add_argument('--cluster', action='store_true',
						help='')
	parser.add_argument('--max_mem', type=str,
						help='')
	parser.add_argument('--qos', type=str,
						help='')
	parser.add_argument('--time', type=str,
						help='')
	parser.add_argument('--log', type=str,
						help='')

	args = parser.parse_args()

	run_snake(args)
