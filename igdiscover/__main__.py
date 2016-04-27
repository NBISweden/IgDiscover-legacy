#!/usr/bin/env python3
"""
IgDiscover computes V/D/J gene usage profiles and discovers novel V genes

- Run IgBLAST in parallel (wrapper inspired by igblastwrp).
- Parse IgBLAST output into a tab-separated table
- Group sequences by barcode
- Plot V gene usage
- Discover new V genes given more than one dataset
"""
__author__ = "Marcel Martin"

import logging
import importlib
from sqt import HelpfulArgumentParser
import matplotlib as mpl
import warnings
from . import __version__

mpl.use('Agg')
warnings.filterwarnings('ignore', 'axes.color_cycle is deprecated and replaced with axes.prop_cycle')
warnings.filterwarnings('ignore', 'The `IPython.html` package')

# List of all subcommands. A module of the given name must exist and define
# add_arguments() and main() functions. Documentation is taken from the first
# line of the module’s docstring.
COMMANDS = [
	'init',
	'run',
	'commonv',
	'igblast',
	'parse',
	'filter',
	'count',
	'group',
	#'multidiscover',
	'compose',
	'discover',
	'clusterplot',
	'errorplot',
	'upstream',
	'dendrogram',
	'rename',
	'union',
]

logger = logging.getLogger(__name__)

def main():
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

	args = parser.parse_args()
	if not hasattr(args, 'func'):
		parser.error('Please provide the name of a subcommand to run')
	else:
		args.func(args)


if __name__ == '__main__':
	main()
