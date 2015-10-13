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
from sqt import HelpfulArgumentParser

from . import (commonv, igblast, parse, filter, count, group, multidiscover, compose,
	discover, init, clusterplot, errorplot, __version__)

logger = logging.getLogger(__name__)

def main():
	logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
	parser = HelpfulArgumentParser(description=__doc__, prog='igypipe')
	parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

	subparsers = parser.add_subparsers()
	init.add_subcommand(subparsers)
	igblast.add_subcommand(subparsers)
	parse.add_subcommand(subparsers)
	filter.add_subcommand(subparsers)
	count.add_subcommand(subparsers)
	group.add_subcommand(subparsers)
	multidiscover.add_subcommand(subparsers)
	commonv.add_subcommand(subparsers)
	errorplot.add_subcommand(subparsers)
	clusterplot.add_subcommand(subparsers)
	discover.add_subcommand(subparsers)
	compose.add_subcommand(subparsers)

	args = parser.parse_args()
	if not hasattr(args, 'func'):
		parser.error("Please provide a command")
	else:
		args.func(args)


if __name__ == '__main__':
	main()
