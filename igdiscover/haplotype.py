"""
Determine haplotypes based on co-occurrences of alleles
"""
import sys
import logging
from typing import List, Tuple, Iterator
import pandas as pd
from argparse import ArgumentParser
from sqt import SequenceReader
from .table import read_table

logger = logging.getLogger(__name__)


# The second-most expressed allele of a gene must be expressed at at least this
# fraction of the highest-expressed allele in order for the gene to be considered
# heterozygous.
HETEROZYGOUS_THRESHOLD = 0.1

#
EXPRESSED_RATIO = 0.1


def add_arguments(parser: ArgumentParser):
	arg = parser.add_argument
	arg('--d-evalue', type=float, default=1E-4,
		help='Maximal allowed E-value for D gene match. Default: %(default)s')
	arg('--d-coverage', '--D-coverage', type=float, default=65,
		help='Minimum D coverage (in percent). Default: %(default)s%%)')
	arg('--het', metavar='GENETYPE', choices=('V', 'D', 'J'), default='J',
		help='Use heterozygous GENETYPE genes as anchoring point for haplotyping. '
			'Default: %(default)s')
	arg('--gene', metavar='GENETYPE', choices=('V', 'D', 'J'),
		help='Haplotype genes of this GENETYPE. Default is to haplotype all except '
			'the --het gene')
	arg('--order', metavar='FASTA', default=None,
		help='Sort the output according to the order of the records in '
			'this FASTA file.')
	arg('table', help='Table with parsed and filtered IgBLAST results')


def expression_counts(table: pd.DataFrame, gene_type: str) -> Iterator[pd.DataFrame]:
	"""
	Yield DataFrames for each gene with gene and allele as the row index and columns 'name' and
	'count'. For example, when 'name' is VH1-1*01, gene would be 'VH1-1' and allele
	would be '01'.
	"""
	counts = table.groupby(gene_type + '_gene').size()
	names, _, alleles = zip(*[s.partition('*') for s in counts.index])
	expressions = pd.DataFrame(
		{'gene': names, 'allele': alleles, 'count': counts, 'name': counts.index},
		columns=['gene', 'allele', 'name', 'count']).set_index(['gene', 'allele'])
	del alleles

	# Example expressions at this point:
	#
	#                             name  count
	# gene       allele
	# IGHV1-18   01        IGHV1-18*01    166
	#            03        IGHV1-18*03      1
	# IGHV1-2    02         IGHV1-2*02     85
	#            04         IGHV1-2*04     16
	# IGHV1-24   01        IGHV1-24*01      5

	logger.info('Heterozygous %s genes:', gene_type)
	for _, alleles in expressions.groupby(level='gene'):
		# Remove alleles that have too low expression relative to the highest-expressed allele
		max_expression = alleles['count'].max()
		alleles = alleles[alleles['count'] >= HETEROZYGOUS_THRESHOLD * max_expression]
		if len(alleles) >= 2:
			logger.info('%s with alleles %s -- Counts: %s',
				alleles.index[0][0],
				', '.join(alleles['name']),
				', '.join(str(x) for x in alleles['count'])
			)
		yield alleles


class GeneMissing(Exception):
	pass


def sorted_haplotype(haplotypes, gene_order):
	gene_order = {name: i for i, name in enumerate(gene_order)}

	def keyfunc(hap):
		name = hap[0] if hap[0] else hap[1]
		gene, _, allele = name.partition('*')
		try:
			allele = int(allele)
		except ValueError:
			allele = 999
		try:
			index = gene_order[gene]
		except KeyError:
			logger.warning('Gene %s not found in gene order file, placing it at the end',
				gene)
			index = 1000000
		return index * 1000 + allele
	return sorted(haplotypes, key=keyfunc)


class HeterozygousGene:
	def __init__(self, name: str, alleles: List[str]):
		"""
		name -- name of this gene, such as 'VH4-4'
		alleles -- list of its alleles, such as ['VH4-4*01', 'VH4-4*02']
		"""
		self.name = name
		self.alleles = alleles


def compute_coexpressions(table: pd.DataFrame, gene_type1: str, gene_type2: str):
	coexpressions = table.groupby(
		(gene_type1 + '_gene', gene_type2 + '_gene')).size().to_frame()
	coexpressions.columns = ['count']
	return coexpressions


def cooccurrences(coexpressions, het_alleles: Tuple[str, str], target_groups):
	"""
	het_alleles -- a pair of alleles of a heterozygous gene,
	such as ('IGHJ6*02', 'IGHJ6*03').
	"""
	assert len(het_alleles) == 2

	haplotype = []
	for target_alleles in target_groups:
		is_expressed_list = []
		names = []
		counts = []
		for target_allele, _ in target_alleles.itertuples(index=False):
			ex = []
			for het_allele in het_alleles:
				try:
					e = coexpressions.loc[(het_allele, target_allele), 'count']
				except KeyError:
					e = 0
				ex.append(e)
			ex_total = sum(ex) + 1  # +1 avoids division by zero
			ratios = [x / ex_total for x in ex]
			is_expressed = [ratio >= EXPRESSED_RATIO for ratio in ratios]
			if is_expressed != [False, False]:
				is_expressed_list.append(is_expressed)
				names.append(target_allele)
				counts.append(ex)
		if len(is_expressed_list) == 1:
			is_expressed = is_expressed_list[0]
			if is_expressed == [True, False]:
				haplotype.append((names[0], '', 'deletion', counts[0]))
			elif is_expressed == [False, True]:
				haplotype.append(('', names[0], 'deletion', counts[0]))
			elif is_expressed == [True, True]:
				haplotype.append((names[0], names[0], 'homozygous', counts[0]))
			else:
				assert False
		elif is_expressed_list == [[True, False], [False, True]]:
			haplotype.append((names[0], names[1], 'heterozygous', (counts[0][0], counts[1][1])))
		elif is_expressed_list == [[False, True], [True, False]]:
			haplotype.append((names[1], names[0], 'heterozygous', (counts[0][1], counts[1][0])))
		else:
			for is_expressed, name, count in zip(is_expressed_list, names, counts):
				haplotype.append((
					name if is_expressed[0] else '',
					name if is_expressed[1] else '',
					'',
					count,
				))
	return haplotype


def read_and_filter(path: str, d_evalue: float, d_coverage: float):
	usecols = ['V_gene', 'D_gene', 'J_gene', 'V_errors', 'D_errors', 'J_errors', 'D_covered',
		'D_evalue']
	# Support reading a table without D_errors
	try:
		table = read_table(path, usecols=usecols)
	except ValueError:
		usecols.remove('D_errors')
		table = read_table(path, usecols=usecols)
	logger.info('Table with %s rows read', len(table))

	table = table[table.V_errors == 0]
	logger.info('%s rows remain after requiring V errors = 0', len(table))
	table = table[table.J_errors == 0]
	logger.info('%s rows remain after requiring J errors = 0', len(table))
	table = table[table.D_evalue <= d_evalue]
	logger.info('%s rows remain after requiring D E-value <= %s', len(table), d_evalue)
	table = table[table.D_covered >= d_coverage]
	logger.info('%s rows remain after requiring D coverage >= %s', len(table), d_coverage)
	if 'D_errors' in table.columns:
		table = table[table.D_errors == 0]
		logger.info('%s rows remain after requiring D errors = 0', len(table))
	return table


def main(args):
	if args.het == args.gene:
		logger.error('Gene types given for --het and --gene must not be the same')
		sys.exit(1)
	if args.order is not None:
		with SequenceReader(args.order) as sr:
			gene_order = [r.name for r in sr]
	else:
		gene_order = None
	het_gene_type = args.het
	if args.gene:
		target_gene_types = [args.gene]
	else:
		target_gene_types = {'V': ['D', 'J'], 'D': ['V', 'J'], 'J': ['V', 'D']}[het_gene_type]

	table = read_and_filter(args.table, args.d_evalue, args.d_coverage)
	het_expressions = {
		gene_type: list(expression_counts(table, gene_type)) for gene_type in ['V', 'D', 'J']}

	coexpressions_dict = {
		gene_type: compute_coexpressions(table, het_gene_type, gene_type)
		for gene_type in ['V', 'D', 'J']}

	for het_alleles in het_expressions[het_gene_type]:
		het_alleles = tuple(het_alleles['name'])
		if len(het_alleles) == 1:
			continue
		if len(het_alleles) != 2:
			logger.warning('More than two alleles found: %s. Skipping!',
				', '.join(het_alleles))
			continue

		logger.info('Using heterozygous %s gene alleles %s and %s to haplotype genes',
			het_gene_type, het_alleles[0], het_alleles[1])
		haplotype = []
		for target_gene_type in target_gene_types:
			coexpressions = coexpressions_dict[target_gene_type]
			target_groups = het_expressions[target_gene_type]
			haplotype.extend(cooccurrences(coexpressions, het_alleles, target_groups))

		if gene_order is not None:
			haplotype = sorted_haplotype(haplotype, gene_order)

		print('# {}:'.format(' vs '.join(het_alleles)))
		print('haplotype1', 'haplotype2', 'type', 'count1', 'count2', sep='\t')
		for h1, h2, typ, count in haplotype:
			print(h1, h2, typ, count[0], count[1], sep='\t')
		print()
