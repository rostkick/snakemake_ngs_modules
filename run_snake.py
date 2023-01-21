import os
import argparse
import pandas as pd
import subprocess as sp
from tempfile import NamedTemporaryFile
import re

from modules.scripts.functions import get_min_substr
from modules.scripts.command_builder import SlurmDecorator, UserInput, CommandBuilder


def run_snake(args):

	germ = sorted(args.fastq_germline, key = lambda x: int(re.search(args.sample_pattern, x).group(1)))
	germline_forward = [i for i in germ if re.match('.*_[\w]1_.*', i) is not None]
	germline_reverse = [i for i in germ if re.match('.*_[\w]2_.*', i) is not None]

	df = pd.DataFrame({'germline_forward': germline_forward,
					'germline_reverse': germline_reverse})

	df.loc[:, 'germline_samples'] = df.loc[:, 'germline_forward'].apply(lambda x: os.path.basename(x)).str.extract('(?P<germline_samples>.+)_R1.*')
	print(df.loc[:, 'germline_samples'])
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
	print(df)
	# print(df.sample(n=3).loc[:, df.columns.str.contains('samples')])
	print('\n')

	df.to_csv(args.output_params, sep='\t', index=False)

	wd_path = os.path.dirname(os.path.realpath(__file__))
	smk_path = os.path.join(wd_path, 'snakefile.smk')

	if args.configfile:
		config_path = args.configfile
	else:
		config_path = os.path.join(wd_path, 'configure.yml')

	
	smk_command_builder = CommandBuilder(UserInput('snakemake', True))

	smk_command_elements = []
	smk_command_elements.append(UserInput('-s', smk_path))
	smk_command_elements.append(UserInput('--nolock', True))
	smk_command_elements.append(UserInput('--configfile', config_path))
	smk_command_elements.append(UserInput('--jobs', args.jobs))
	smk_command_elements.append(UserInput('--cores', args.cores))
	smk_command_elements.append(UserInput('--config', {
													'run_id': args.run_id, 
													'mode': args.mode, 
													'params_table': args.output_params
													}))
	smk_command_elements.append(UserInput('-p', args.printshellcmds))
	smk_command_elements.append(UserInput('-F', args.force))
	smk_command_elements.append(UserInput('-q', args.quite))
	smk_command_elements.append(UserInput('-n', args.dryrun))
	smk_command_elements.append(UserInput('--use-conda', args.use_conda))
	smk_command_elements.append(UserInput('--ri', args.rerun_incomplete))
	smk_command_elements.append(UserInput('--dag', args.dag))
	smk_command_elements.append(UserInput('--rulegraph', args.rulegraph))
	smk_command_elements.append(UserInput('-R', args.forcerun))
	smk_command_elements.append(UserInput('-U', args.until))

	for element_x in smk_command_elements:
		smk_command_builder.add_cmd_element(element_x)

	snake_command = str(smk_command_builder)

	if args.cluster & ~(args.dag | args.rulegraph | args.dryrun):

		slurm_command_builder = CommandBuilder(UserInput('#!/bin/sh', True))
		slurm_command_builder.set_str_end_root('')

		slurm_command_elements = []
		slurm_command_elements.append(UserInput('sbatch << ENDINPUT', True))
		slurm_command_elements.append(UserInput('#!/bin/sh', True))
		slurm_command_elements.append(UserInput("#SBATCH --job-name", args.run_id))
		slurm_command_elements.append(UserInput("#SBATCH --cpus-per-task", args.jobs))
		slurm_command_elements.append(UserInput("#SBATCH --mem", args.max_mem))
		slurm_command_elements.append(UserInput("#SBATCH --qos", args.qos))
		slurm_command_elements.append(UserInput("#SBATCH --time", args.time))
		slurm_command_elements.append(UserInput("#SBATCH --output", args.log))
		slurm_command_elements.append(UserInput("#SBATCH --priority", "1"))

		for element_x in slurm_command_elements:
			slurm_command_builder.add_cmd_element(element_x)

		slurm_command_builder.set_str_end('')
		slurm_command_builder.set_indent_size(0)

		sbatch_task = SlurmDecorator(smk_command_builder, slurm_command_builder)

		tmp = NamedTemporaryFile()
		with open(tmp.name, 'w') as w:
			w.write(str(sbatch_task))
		print(sbatch_task)
		sp.run(f"bash {tmp.name}", shell=True)
	else:
		print(snake_command)
		sp.run(snake_command, shell=True)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Process some integers.')
	parser.add_argument('--configfile', '-cf', type=str,
						help='')
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
	parser.add_argument('--jobs', '-j', required=True,
						help='')
	parser.add_argument('--cores', '-c', required=True,
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
