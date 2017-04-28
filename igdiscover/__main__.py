#!/usr/bin/env python3
"""
IgDiscover computes V/D/J gene usage profiles and discovers novel V genes

- Run IgBLAST in parallel (wrapper inspired by igblastwrp).
- Parse IgBLAST output into a tab-separated table
- Group sequences by barcode
- Plot V gene usage
- Discover new V genes given more than one dataset
"""
import sys
import logging
import importlib
from sqt import HelpfulArgumentParser
import matplotlib as mpl
import warnings
import resource
from . import __version__

__author__ = "Marcel Martin"

mpl.use('Agg')
warnings.filterwarnings('ignore', 'axes.color_cycle is deprecated and replaced with axes.prop_cycle')
warnings.filterwarnings('ignore', 'The `IPython.html` package')

# List of all subcommands. A module of the given name must exist and define
# add_arguments() and main() functions. Documentation is taken from the first
# line of the module’s docstring.
COMMANDS = [
	'init',
	'run',
	'config',
	'commonv',
	'igblast',
	'parse',
	'filter',
	'count',
	'group',
	#'multidiscover',
	'germlinefilter',
	'discover',
	'discoverj',
	'clusterplot',
	'errorplot',
	'upstream',
	'dendrogram',
	'rename',
	'union',
	'clonotypes',
	'plotalleles',
]

logger = logging.getLogger(__name__)


def format_duration(seconds):
	h = int(seconds // 3600)
	seconds -= h * 3600
	m = int(seconds // 60)
	seconds -= m * 60
	return '{:02d}:{:02d}:{:04.1f}'.format(h, m, seconds)


def main(arguments=None):
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	parser = HelpfulArgumentParser(description=__doc__, prog='igdiscover')
	parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

	subparsers = parser.add_subparsers()
	for command_name in COMMANDS:
		module = importlib.import_module('.' + command_name, 'igdiscover')
		subparser = subparsers.add_parser(command_name,
			help=module.__doc__.split('\n')[1], description=module.__doc__)
		subparser.set_defaults(func=module.main)
		module.add_arguments(subparser)

	args = parser.parse_args(arguments)
	if not hasattr(args, 'func'):
		parser.error('Please provide the name of a subcommand to run')
	else:
		args.func(args)
	if sys.platform == 'linux':
		rself = resource.getrusage(resource.RUSAGE_SELF)
		rchildren = resource.getrusage(resource.RUSAGE_CHILDREN)
		memory_kb = rself.ru_maxrss + rchildren.ru_maxrss
		cpu_time = rself.ru_utime + rself.ru_stime + rchildren.ru_utime + rchildren.ru_stime
		cpu_time_s = format_duration(cpu_time)
		logger.info('CPU time {}. Maximum memory usage {:.3f} GB'.format(
			cpu_time_s, memory_kb / 1E6))


if __name__ == '__main__':
	main()
