"""
Haplotype
"""
import sys
import logging
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

	for gene_type in ('V', 'D', 'J'):
		print('Heterozygous', gene_type, 'genes:')
		for _, group in expressions[gene_type].groupby(level='gene'):
			group = filter_alleles(group)
			if len(group) >= 2:
				print(group.index[0][0], 'with alleles', ', '.join(group['name']), '-- Counts:',
					', '.join(str(x) for x in group['count'])
				)
		print()

	return

	def count_heterozygous(groups):
		return sum(1 for g in groups if len(g) > 1)

	def filtered_group(expressions):
		result = []
		for _, group in expressions.groupby(level='gene'):
			result.append(filter_alleles(group))
		return result

	for gene_type1, gene_type2 in [('V', 'D'), ('D', 'J'), ('V', 'J')]:
		groups1 = filtered_group(expressions[gene_type1])
		groups2 = filtered_group(expressions[gene_type2])
		n_het1 = count_heterozygous(groups1)
		n_het2 = count_heterozygous(groups2)

		# Swap if necessary to make the first group the smaller one
		if n_het1 > n_het2:
			groups1, groups2 = groups2, groups1
			n_het1, n_het2 = n_het2, n_het1
			gene_type1, gene_type2 = gene_type2, gene_type1

		expressions1 = expressions[gene_type1]
		expressions2 = expressions[gene_type2]
		print('=== Inspecting', n_het1, 'heterozygous', gene_type1, 'genes and', n_het2,
			'heterozygous', gene_type2, 'genes')
		print()
		coexpressions = table.groupby(
			(gene_type1 + '_gene', gene_type2 + '_gene')).size().to_frame()
		coexpressions.columns = ['count']
		hom = []
		for alleles1 in groups1:
			# print('a1', alleles1)
			if len(alleles1) == 1:
				hom.append(alleles1.iloc[0]['name'])  # , alleles1.iloc[0]['count'])
				continue
			print(alleles1.index[0][0], 'is heterozygous')
			# heterozygous gene, look at its alleles
			for name1, _ in alleles1.itertuples(index=False):
				print('Allele', name1, 'occurs with:')
				# iterate over the genes of other type
				for alleles2 in groups2:
					if len(alleles2) < 2:
						continue
					ex = []
					for name2, _ in alleles2.itertuples(index=False):
						try:
							e = coexpressions.loc[(name1, name2), 'count']
						except KeyError:
							e = 0
						ex.append(e + 1)  # +1 for pseudo-count
					# print('co-expression with', name2, e)
					ex_total = sum(ex)
					ratios = [x / ex_total for x in ex]
					# print(ex)
					# print(ratios)
					for (name, _), ratio, e in zip(alleles2.itertuples(index=False), ratios, ex):
						if ratio > 0.25:
							print(name, e)
						else:
							print('NOT', name, e)
					print()
		print('homozygous:', ', '.join(hom))
