"""
Run IgDiscover

Calls Snakemake to produce all the output files.
"""
import os
import logging
import subprocess
import pkg_resources
from snakemake import snakemake
from sqt.utils import available_cpu_count

logger = logging.getLogger(__name__)

def add_arguments(parser):
	add = parser.add_argument
	add('--dryrun', '-n', default=False, action='store_true',
		help='Do not execute anything')
	add('--cores', '--jobs', '-j', metavar='N', type=int, default=available_cpu_count(),
		help='Run on at most N CPU cores in parallel. '
		'Default: Use as many cores as available (%(default)s)')
	add('--keepgoing', '-k', default=False, action='store_true',
		help='If one job fails, finish the others.')
	add('--printshellcmds', '-p', default=False, action='store_true',
		help='Print out the shell commands that will be executed.')


def subprocess_snakemake(snakefile_path, dryrun, cores, keepgoing, printshellcmds):
	"""
	Call snakemake. There are some bugs in its Python API, so this routine
	just spawns a subprocess.
	"""
	arguments = ['snakemake', '-s', snakefile_path, '-j', str(cores)]
	if dryrun:
		arguments += ['-n']
	if printshellcmds:
		arguments += ['-p']
	subprocess.call(arguments)


def main(args):
	# snakemake sets up its own logging and this cannot be easily changed
	# (setting keep_logger=True crashes), so remove our own log handler
	# for now
	logger.root.handlers = []
	snakefile_path = pkg_resources.resource_filename('igdiscover', 'Snakefile')
	snakemake(snakefile_path,
		dryrun=args.dryrun,
		cores=args.cores,
		keepgoing=args.keepgoing,
		printshellcmds=args.printshellcmds,
	)
