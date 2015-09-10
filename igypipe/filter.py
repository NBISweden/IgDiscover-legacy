"""
Filter table with parsed IgBLAST results

The filtered table is printed to standard output.
"""
import logging
from itertools import islice

from .table import read_table, filtered_table

logger = logging.getLogger(__name__)

def add_subcommand(subparsers):
	subparser = subparsers.add_parser('filter', help=__doc__.split('\n')[1], description=__doc__)
	subparser.set_defaults(func=filter_command)
	subparser.add_argument('table', help='Table with filtered IgBLAST results.')
	return subparser


def filter_command(args):
	d = read_table(args.table, log=True)
	d = filtered_table(d, log=True)
	print(d.to_csv(sep='\t', index=False))
	logger.info('%d rows written', len(d))
