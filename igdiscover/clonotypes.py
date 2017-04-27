"""
Create table of clonotypes (preliminary)
"""
import logging

from .table import read_table

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('table', help='Table with parsed and filtered IgBLAST results')


def main(args):
	table = read_table(args.table)
	table['CDR3_length'] = table['CDR3_nt'].apply(len)
#	table.sort_values(['V_gene', 'J_gene', 'CDR3_nt'])

	if False:
		print(table.head(0).to_csv(sep='\t', header=True, index=False), end='')
		for name, group in table.groupby(('V_gene', 'J_gene', 'CDR3_nt')):
			# We get an intentional empty line between groups since
			# to_csv() already includes a line break
			print(group.to_csv(sep='\t', header=False, index=False))

	columns = ['V_gene', 'J_gene', 'CDR3_length', 'CDR3_nt']
	print('group_size', *columns, sep='\t')
	for name, group in table.groupby(columns):
		print(group['count'].sum(), *name, sep='\t')
