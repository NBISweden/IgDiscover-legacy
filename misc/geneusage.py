#!/usr/bin/env python3
"""
Compute and plot V/J gene usage from a clonotypes.tab file
"""
import sys
from argparse import ArgumentParser
import seaborn
import pandas as pd


def add_arguments(parser):
	arg = parser.add_argument
	arg('--mincount', metavar='N', default=1, type=int,
		help='Filter out clonotypes with less than N members. Default: %(default)s')
	arg('--gene', default='V', choices=('V', 'J'),
		help='Which gene to plot. Choose V or J. Default: %(default)s')
	arg('--plot', metavar='FILE',
		help='Plot results into FILE (PNG, PDF)')
	arg('table', help='clonotypes.tab table created with the "clonotypes" subcommand')


def main(args):
	# Read in input table
	table = pd.read_csv(args.table, sep='\t')

	# Extract the 'V_gene' or 'J_gene' column
	counts_column = table[args.gene + '_gene']

	# Count how often each value (gene) appears
	gene_frequencies = counts_column.value_counts()

	# Keep only those with minimum count
	gene_frequencies = gene_frequencies[gene_frequencies >= args.mincount]

	# Sort by gene name
	gene_frequencies = gene_frequencies.sort_index()

	# Print out the frequencies
	print('gene', 'frequency', sep='\t')
	for gene_name, frequency in gene_frequencies.items():
		print(gene_name, frequency, sep='\t')

	# If requested, create a plot
	if args.plot:
		ax = gene_frequencies.plot(kind='bar', figsize=(20, 5))
		ax.set_title('Gene expression counts')
		ax.set_xlabel('{} gene'.format(args.gene))
		ax.set_ylabel('Clonotype count')
		ax.figure.tight_layout()
		ax.figure.savefig(args.plot)



if __name__ == '__main__':
	parser = ArgumentParser(description=__doc__)
	add_arguments(parser)
	args = parser.parse_args()
	main(args)
