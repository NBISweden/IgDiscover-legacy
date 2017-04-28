"""
Plot allele usage
"""
import sys
import logging
import pandas as pd
from sqt import SequenceReader
from .table import read_table

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	arg('--d-evalue', type=float, default=1E-4,
		help='Maximal allowed E-value for D gene match. Default: %(default)s')
	arg('--d-coverage', '--D-coverage', type=float, default=65,
		help='Minimum D coverage (in percent). Default: %(default)s%%)')
	arg('--vdatabase', metavar='FASTA',
		help='Restrict plotting to the V sequences named in the FASTA file. '
		'Only the sequence names are used!')
	arg('--gene', choices=('D', 'J'), help='Type of gene on y axis. Default: %(default)s', default='J')
	arg('alleles', help='List of alleles to plot on y axis, separated by comma')
	arg('table', help='Table with parsed and filtered IgBLAST results')
	arg('plot', help='Path to output PDF or PNG')


def main(args):
	usecols = ['V_gene', 'D_gene', 'J_gene', 'V_errors', 'D_errors', 'J_errors', 'D_covered',
		'D_evalue']
	# Support reading a table without D_errors
	try:
		table = read_table(args.table, usecols=usecols)
	except ValueError:
		usecols = [col for col in usecols if col != 'D_errors']
		table = read_table(args.table, usecols=usecols)
	logger.info('Table with %s rows read', len(table))
	table = table[table.V_errors == 0]
	logger.info('%s rows remain after requiring V errors = 0', len(table))
	if args.gene == 'J':
		table = table[table.J_errors == 0]
		logger.info('%s rows remain after requiring J errors = 0', len(table))
	elif args.gene == 'D':
		table = table[table.D_evalue <= args.d_evalue]
		logger.info('%s rows remain after requiring D E-value <= %s', len(table), args.d_evalue)
		table = table[table.D_covered >= args.d_coverage]
		logger.info('%s rows remain after requiring D coverage >= %s', len(table), args.d_coverage)
		if 'D_errors' in table.columns:
			table = table[table.D_errors == 0]
			logger.info('%s rows remain after requiring D errors = 0', len(table))

	gene1 = 'V_gene'
	gene2 = args.gene + '_gene'
	expression_counts = table.groupby((gene1, gene2)).size().to_frame().reset_index()
	matrix = pd.DataFrame(
		expression_counts.pivot(index=gene1, columns=gene2, values=0).fillna(0), dtype=int)
	# matrix[v_gene,d_gene] gives co-occurrences of v_gene and d_gene
	print('#\n# Expressed genes with counts\n#')
	# The .sum() is along axis=0, that is, the V gene counts are summed up,
	# resulting in counts for each D/J gene
	for g, count in matrix.sum().iteritems():
		print(g, '{:8}'.format(count))

	alleles = args.alleles.split(',')
	for allele in alleles:
		if allele not in matrix.columns:
			logger.error('Allele %s not expressed in this dataset', allele)
			sys.exit(1)

	if args.vdatabase:
		with SequenceReader(args.vdatabase) as f:
			v_names = [record.name for record in f if record.name in matrix.index]
		if not v_names:
			logger.error('None of the sequence names in %r were found in the input table',
				args.vdatabase)
			sys.exit(1)
		allele_expressions = matrix.loc[v_names, alleles]
	else:
		allele_expressions = matrix.loc[:, alleles]
	print('#\n# Allele-specific expression\n#')
	print(allele_expressions)

	if len(alleles) == 2:
		allele_expressions.loc[:, alleles[1]] *= -1

	# remove all-zero rows
	# allele_expressions = allele_expressions[(allele_expressions > 0.001).any(axis=1)]
	ax = allele_expressions.plot(kind='bar', stacked=True, figsize=(12, 6))
	ax.legend(title=None)
	ax.set_title('Allele-specific expression counts')
	ax.set_xlabel('V gene')
	ax.set_ylabel('Count')
	ax.figure.set_tight_layout(True)
	# ax.legend(bbox_to_anchor=(1.15, 0.5))
	ax.figure.savefig(args.plot)
	logger.info('Plotted %r', args.plot)
