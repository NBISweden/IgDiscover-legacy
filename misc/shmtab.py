#!/usr/bin/env python3
"""
Extract V_SHM values from a filtered.tab.gz and group them by V gene

Create a table that has a column for each V gene.
"""
import logging
import pandas as pd


def add_arguments(parser):
	arg = parser.add_argument
	arg('--field', default='V_SHM', help='Which column to extract. Default: %(default)s')
	arg('table', help='filtered.tab.gz file')


def main(args):
	table = pd.read_csv(args.table, sep='\t', usecols=(args.field, 'V_gene'), sep='\t')
	df = pd.DataFrame(index=range(table.groupby('V_gene').size().max()))
	for gene, group in table.groupby('V_gene'):
		df[gene] = group[args.field].reset_index(drop=True)
	print(df.to_csv(sep='\t', index=False))


if __name__ == '__main__':
	from argparse import ArgumentParser
	parser = ArgumentParser()
	add_arguments(parser)
	args = parser.parse_args()
	main(args)
