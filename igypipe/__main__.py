#!/usr/bin/env python3
"""
Antibody pipeline helper script. It can:

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

from . import __version__

COMMANDS = [
	'commonv',
	'igblast',
	'parse',
	'filter',
	'count',
	'group',
	'multidiscover',
	'compose',
	'discover',
	'init',
	'clusterplot',
	'errorplot',
	'upstream'
]

logger = logging.getLogger(__name__)

def main():
	logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
	parser = HelpfulArgumentParser(description=__doc__, prog='igypipe')
	parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

	subparsers = parser.add_subparsers()
	for command_name in COMMANDS:
		module = importlib.import_module('.' + command_name, 'igypipe')
		subparser = subparsers.add_parser(command_name,
			help=module.__doc__.split('\n')[1], description=module.__doc__)
		subparser.set_defaults(func=module.main)
		module.add_arguments(subparser)

	args = parser.parse_args()
	if not hasattr(args, 'func'):
		parser.error("Please provide a command")
	else:
		args.func(args)


if __name__ == '__main__':
	main()
