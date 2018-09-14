#!/usr/bin/env python3
"""
Extract V_SHM values from multiple clonoquery tables
"""
import logging
import pandas as pd
import io
import os
import sys


def add_arguments(parser):
	arg = parser.add_argument
	arg('--field', default='V_SHM', help='Which column to extract. Default: %(default)s')
	arg('-n', default=5, type=int, help='Require at least N values in each file for a query to be listed')
	arg('paths', nargs='+')


def read_one_file(path, field_name):
	with open(path) as f:
		header = f.readline()
		sections = f.read().strip().split('\n\n')
	column_names = header.strip().split('\t')
	sio = io.StringIO(sections[0])
	sio.readline()

	queries = dict()
	for section in sections:
		sio = io.StringIO(section)
		query = sio.readline()
		if not query.startswith('#'):
			print(f'query is: {query!r}')
		assert query.startswith('# Query: '), query
		query_name = query.split('\t')[0].partition('# Query: ')[2]
		table = pd.read_table(sio, header=None, names=column_names, usecols=(field_name,))
		queries[query_name] = table[field_name]
	return queries


def main(args):
	files = {os.path.basename(path): read_one_file(path, args.field) for path in args.paths}
	print('read {} files'.format(len(files)), file=sys.stderr)

	query_names = dict()
	for queries in files.values():
		for query_name in queries.keys():
			query_names[query_name] = None
	query_names = list(query_names)

	maxlen = 0
	for queries in files.values():
		maxlen = max(maxlen, max(len(q) for q in queries.values()))

	df = pd.DataFrame(index=range(maxlen))
	for query_name in query_names:
		long_enough_in_all_files = all(len(files[file].get(query_name, [])) >= args.n for file in files)
		if not long_enough_in_all_files:
			continue
		for file in files:
			column_name = '{file}-{query_name}'.format(file=file, query_name=query_name)
			df[column_name] = files[file].get(query_name, pd.Series([]))
	print(df.to_csv(sep='\t', index=False))


if __name__ == '__main__':
	from argparse import ArgumentParser
	parser = ArgumentParser()
	add_arguments(parser)
	args = parser.parse_args()
	main(args)
