"""
Filter table with parsed IgBLAST results

Discard the following rows in the table:
- no J assigned
- stop codon found
- V gene coverage less than 90
- J gene coverage less than 60
- V gene E-value greater than 1E-3

The filtered table is printed to standard output.
"""
import logging
from itertools import islice
import pandas as pd

from .table import read_table

logger = logging.getLogger(__name__)

def add_arguments(parser):
	parser.add_argument('table', help='Table with filtered IgBLAST results.')


def filtered_table(table,
		v_gene_coverage=90,  # at least
		j_gene_coverage=60,  # at least
		v_gene_evalue=1E-3,  # at most
	):
	"""
	Discard the following rows in the table (read in by read_table):
	- no J assigned
	- stop codon found
	- V gene coverage less than v_gene_coverage
	- J gene coverage less than j_gene_coverage
	- V gene E-value greater than v_gene_evalue

	Return the filtered table.
	"""
	# Both V and J must be assigned
	# (Note V_gene and J_gene columns use empty strings instead of NA)
	filtered = table[(table['V_gene'] != '') & (table['J_gene'] != '')][:]
	logger.info('%s rows have both V and J assignment', len(filtered))
	filtered['V_gene'] = pd.Categorical(filtered['V_gene'])

	# Filter out sequences that have a stop codon
	filtered = filtered[filtered.stop == 'no']
	logger.info('%s of those do not have a stop codon', len(filtered))

	# Filter out sequences with a too low V gene hit E-value
	filtered = filtered[filtered.V_evalue <= v_gene_evalue]
	logger.info('%s of those have an E-value of at most %s', len(filtered), v_gene_evalue)

	# Filter out sequences with too low V gene coverage
	filtered = filtered[filtered.V_covered >= v_gene_coverage]
	logger.info('%s of those cover the V gene by at least %s%%', len(filtered), v_gene_coverage)

	# Filter out sequences with too low J gene coverage
	filtered = filtered[filtered.J_covered >= j_gene_coverage]
	logger.info('%s of those cover the J gene by at least %s%%', len(filtered), j_gene_coverage)

	return filtered


def main(args):
	d = read_table(args.table, log=True)
	d = filtered_table(d)
	print(d.to_csv(sep='\t', index=False), end='')
	logger.info('%d rows written', len(d))
