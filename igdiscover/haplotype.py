"""
Haplotype
"""
import sys
import logging
from collections import defaultdict
import pandas as pd
from sqt import SequenceReader
from .table import read_table

logger = logging.getLogger(__name__)


# The second-most expressed allele of a gene must be expressed at at least this
# fraction of the highest-expressed allele in order for the gene to be considered
# heterozygous.
HETEROZYGOUS_THRESHOLD = 0.1


def add_arguments(parser):
	arg = parser.add_argument
	arg('--d-evalue', type=float, default=1E-4,
		help='Maximal allowed E-value for D gene match. Default: %(default)s')
	arg('--d-coverage', '--D-coverage', type=float, default=65,
		help='Minimum D coverage (in percent). Default: %(default)s%%)')
	arg('--database', metavar='FASTA',
		help='Restrict plotting to the sequences named in the FASTA file. '
		'Only the sequence names are used!')
	arg('--x', choices=('V', 'D', 'J'), default='V',
		help='Type of gene on x axis. Default: %(default)s')
	arg('--gene', choices=('V', 'D', 'J'), default='J',
		help='Type of gene on y axis. Default: %(default)s')
	arg('table', help='Table with parsed and filtered IgBLAST results')


def expression_counts(table, gene_type):
	"""
	Return a DataFrame with gene and allele as the row index and columns 'name' and
	'count'. When 'name' is VH1-1*01, gene would be 'VH1-1' and allele
	would be '01'.
	"""
	counts = table.groupby(gene_type + '_gene').size()
	names, _, alleles = zip(*[s.partition('*') for s in counts.index])
	expressions = pd.DataFrame(
		{'gene': names, 'allele': alleles, 'count': counts, 'name': counts.index},
		columns=['gene', 'allele', 'name', 'count']).set_index(['gene', 'allele'])
	return expressions


def filter_alleles(table):
	"""
	Remove alleles that have too low expression relative to the highest-expressed allele
	"""
	max_expression = table['count'].max()
	return table[table['count'] >= HETEROZYGOUS_THRESHOLD * max_expression]


def cooccurrences(table, gene_type1, gene_type2, groups1, groups2):
	logger.info('Analyzing heterozygous %s genes and comparing to %s genes', gene_type1, gene_type2)
	coexpressions = table.groupby(
		(gene_type1 + '_gene', gene_type2 + '_gene')).size().to_frame()
	coexpressions.columns = ['count']
	hom = []
	for alleles1 in groups1:
		if len(alleles1) == 1:
			hom.append(alleles1.iloc[0]['name'])
			continue
		assert len(alleles1) == 2  # TODO
		print('# Haplotypes from {}:'.format('/'.join(alleles1['name'])))
		print('haplotype1', 'haplotype2', 'type', sep='\t')
		haplotypes = []
		for alleles2 in groups2:
			is_expressed_list = []
			names = []
			for name2, _ in alleles2.itertuples(index=False):
				ex = []
				for name1, _ in alleles1.itertuples(index=False):
					try:
						e = coexpressions.loc[(name1, name2), 'count']
					except KeyError:
						e = 0
					ex.append(e)
				ex_total = sum(ex) + 1  # +1 avoids division by zero
				ratios = [x / ex_total for x in ex]
				is_expressed = [ratio >= 0.1 for ratio in ratios]
				is_expressed_list.append(is_expressed)
				names.append(name2)
			if len(is_expressed_list) == 1:
				is_expressed = is_expressed_list[0]
				if is_expressed == [True, False]:
					haplotypes.append((names[0], '', 'deletion'))
				elif is_expressed == [False, True]:
					haplotypes.append(('', names[0], 'deletion'))
				elif is_expressed == [True, True]:
					haplotypes.append((names[0], names[0], 'homozygous'))
				else:
					assert is_expressed == [False, False]
			elif is_expressed_list == [[True, False], [False, True]]:
				haplotypes.append((names[0], names[1], 'heterozygous'))
			elif is_expressed_list == [[False, True], [True, False]]:
				haplotypes.append((names[1], names[0], 'heterozygous'))
			else:
				for is_expressed, name in zip(is_expressed_list, names):
					if is_expressed == [False, False]:
						continue
					haplotypes.append((
						name if is_expressed[0] else '',
						name if is_expressed[1] else '', ''))
		for h1, h2, typ in haplotypes:
			print(h1, h2, typ, sep='\t')
	logger.info('homozygous: %s', ', '.join(hom))


def main(args):
	usecols = ['V_gene', 'D_gene', 'J_gene', 'V_errors', 'D_errors', 'J_errors', 'D_covered',
		'D_evalue']
	# Support reading a table without D_errors
	try:
		table = read_table(args.table, usecols=usecols)
	except ValueError:
		usecols.remove('D_errors')
		table = read_table(args.table, usecols=usecols)
	logger.info('Table with %s rows read', len(table))

	table = table[table.V_errors == 0]
	logger.info('%s rows remain after requiring V errors = 0', len(table))
	table = table[table.J_errors == 0]
	logger.info('%s rows remain after requiring J errors = 0', len(table))
	table = table[table.D_evalue <= args.d_evalue]
	logger.info('%s rows remain after requiring D E-value <= %s', len(table), args.d_evalue)
	table = table[table.D_covered >= args.d_coverage]
	logger.info('%s rows remain after requiring D coverage >= %s', len(table), args.d_coverage)
	if 'D_errors' in table.columns:
		table = table[table.D_errors == 0]
		logger.info('%s rows remain after requiring D errors = 0', len(table))

	# Pre-compute expression levels
	expressions = {gt: expression_counts(table, gt) for gt in ('V', 'D', 'J')}

	heterozygous = defaultdict(list)
	for gene_type in ('V', 'D', 'J'):
		logger.info('Heterozygous %s genes:', gene_type)
		for _, group in expressions[gene_type].groupby(level='gene'):
			group = filter_alleles(group)
			if len(group) >= 2:
				logger.info('%s with alleles %s -- Counts: %s',
					group.index[0][0],
					', '.join(group['name']),
					', '.join(str(x) for x in group['count'])
				)

	def count_heterozygous(groups):
		return sum(1 for g in groups if len(g) > 1)

	def filtered_group(gene_type):
		"""Group on gene level, filter out low-expressed genes in each group"""
		result = []
		for _, group in expressions[gene_type].groupby(level='gene'):
			result.append(filter_alleles(group))
		return result

	for gene_type1, gene_type2 in [('J', 'V')]:  # TODO [('V', 'D'), ('D', 'J'), ('V', 'J')]:
		# groups1 and groups2 are expressions grouped by gene (one row is an allele)
		groups1 = filtered_group(gene_type1)
		groups2 = filtered_group(gene_type2)

		cooccurrences(table, gene_type1, gene_type2, groups1, groups2)
		# TODO cooccurrences(table, gene_type2, gene_type1, groups2, groups1)
