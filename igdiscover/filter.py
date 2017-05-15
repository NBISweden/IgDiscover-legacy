"""
Filter table with parsed IgBLAST results

Discard the following rows in the table:
- no J assigned
- stop codon found
- V gene coverage less than 90%
- J gene coverage less than 60%
- V gene E-value greater than 1E-3

The filtered table is printed to standard output.
"""
import logging
import pandas as pd

from .table import fix_columns

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--v-coverage', type=float, default=90, metavar='PERCENT',
		help='Require that the sequence covers at least PERCENT of the V gene. '
		'Default: %(default)s')
	arg('--v-evalue', type=float, default=1E-3, metavar='EVALUE',
		help='Require that the E-value for the V gene match is at most EVALUE. '
		'Default: %(default)s')
	arg('--j-coverage', type=float, default=60, metavar='PERCENT',
		help='Require that the sequence covers at least PERCENT of the J gene. '
		'Default: %(default)s')
	arg('table', help='Table with filtered IgBLAST results.')


class FilteringStatistics:
	__slots__ = ('n', 'vjassigned', 'stop', 'v_evalue', 'v_coverage', 'j_coverage')

	def __init__(self):
		self.n = 0
		self.vjassigned = 0
		self.stop = 0
		self.v_evalue = 0
		self.v_coverage = 0
		self.j_coverage = 0

	def __iadd__(self, other):
		for att in self.__slots__:
			v = getattr(self, att)
			setattr(self, att, v + getattr(other, att))
		return self


def filtered_table(table,
		v_gene_coverage,  # at least
		j_gene_coverage,  # at least
		v_gene_evalue,  # at most
	):
	"""
	Discard the following rows in the table:
	- no J assigned
	- stop codon found
	- V gene coverage less than v_gene_coverage
	- J gene coverage less than j_gene_coverage
	- V gene E-value greater than v_gene_evalue

	Return the filtered table.
	"""
	stats = FilteringStatistics()
	stats.n = len(table)
	# Both V and J must be assigned
	# (Note V_gene and J_gene columns use empty strings instead of NA)
	filtered = table[(table['V_gene'] != '') & (table['J_gene'] != '')][:]
	stats.vjassigned = len(filtered)
	filtered['V_gene'] = pd.Categorical(filtered['V_gene'])

	# Filter out sequences that have a stop codon
	filtered = filtered[filtered.stop == 'no']
	stats.stop = len(filtered)

	# Filter out sequences with a too low V gene hit E-value
	filtered = filtered[filtered.V_evalue <= v_gene_evalue]
	stats.v_evalue = len(filtered)

	# Filter out sequences with too low V gene coverage
	filtered = filtered[filtered.V_covered >= v_gene_coverage]
	stats.v_coverage = len(filtered)

	# Filter out sequences with too low J gene coverage
	filtered = filtered[filtered.J_covered >= j_gene_coverage]
	stats.j_coverage = len(filtered)

	return filtered, stats


def main(args):
	n = 0
	first = True
	written = 0
	stats = FilteringStatistics()
	for chunk in pd.read_csv(args.table, chunksize=10000, sep='\t'):
		fix_columns(chunk)
		n += len(chunk)
		filtered, chunk_stats = filtered_table(chunk, v_gene_coverage=args.v_coverage,
			j_gene_coverage=args.j_coverage, v_gene_evalue=args.v_evalue)
		stats += chunk_stats
		print(filtered.to_csv(sep='\t', index=False, header=first), end='')
		first = False
		written += len(filtered)

	logger.info('%s rows in input table', stats.n)
	logger.info('%s rows have both V and J assignment', stats.vjassigned)
	logger.info('%s of those do not have a stop codon', stats.stop)
	logger.info('%s of those have an E-value of at most %s', stats.v_evalue, args.v_evalue)
	logger.info('%s of those cover the V gene by at least %s%%', stats.v_coverage, args.v_coverage)
	logger.info('%s of those cover the J gene by at least %s%%', stats.j_coverage, args.j_coverage)
	logger.info('%d rows written', written)
