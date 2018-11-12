"""
Compare two germline gene databases given as FASTA files
"""
import sys
import logging
from sqt import FastaReader
from sqt.align import hamming_distance
#from .utils import natural_sort_key

logger = logging.getLogger(__name__)


def add_arguments(parser):
	arg = parser.add_argument
	# arg('--no-sort', dest='sort', action='store_false', default=True,
	# 	help='Do not sort sequences by name')
	# arg('--not-found', metavar='TEXT', default=' (not found)',
	# 	help='Append this text to the record name when the sequence was not found '
	# 	'in the template. Default: %(default)r')
	# arg('--rename-from', metavar='TEMPLATE',
	# 	help='FASTA template file with correctly named sequences. If a sequence in '
	# 		'the target file is identical to one in the template, it is assigned the '
	# 		'name of the sequence in the template.')
	# arg('--order-by', metavar='TEMPLATE',
	# 	help='FASTA template that contains genes in the desired order. '
	# 		'If a name contains a "*" (asterisk), only the preceding part '
	# 		'is used. Thus, VH4-4*01 and VH4-4*02 are equivalent.')
	arg('a', help='FASTA file with expected sequences')
	arg('b', help='FASTA file with actual sequences')


class PrefixDict:
	"""
	A dict that maps strings to values, but where a prefix of a key is enough
	to retrieve the value.
	"""
	def __init__(self, items):
		self._items = []
		for k, v in items:
			self.add(k, v)

	def add(self, k, v):
		try:
			self[k]
			raise ValueError('Key {!r} already exists'.format(k))
		except KeyError:
			self._items.append((k, v))

	def __getitem__(self, key):
		found = None
		for seq, value in self._items:
			if seq.startswith(key) or key.startswith(seq):
				if found is not None:
					# TODO don't use keyerror here
					raise KeyError('Key {!r} is ambiguous'.format(key))
				found = value
		if found is None:
			raise KeyError(key)
		return found

	def get(self, key, default=None):
		try:
			v = self[key]
		except KeyError:
			return default
		else:
			return v

	def __len__(self):
		return len(self._items)


class GeneMissing(Exception):
	pass


def gene_name(record):
	record_name, _, record_comment = record.name.partition(' ')
	return record_name.partition('*')[0]


def sorted_by_gene(records, gene_order):
	d = {name: i for i, name in enumerate(gene_order)}

	def keyfunc(record):
		gene = gene_name(record)
		try:
			return d[gene]
		except KeyError:
			raise GeneMissing(gene)

	return sorted(records, key=keyfunc)


def check_duplicate_names(records):
	names = set()
	for record in records:
		if record.name in names:
			yield record.name
		names.add(record.name)


def check_exact_duplicate_sequences(records):
	sequences = dict()
	for record in records:
		if record.sequence in sequences:
			yield record.name, sequences[record.sequence]
		else:
			sequences[record.sequence] = record.name


def pair_up(a_records, b_records):
	b_records = b_records[:]
	only_a = []
	pairs = []
	# Try to find a partner (either identical or similar) for every A sequence
	for a in a_records:
		best_score = -10000
		best = None
		for b in b_records:
			l = min(len(a.sequence), len(b.sequence))
			length_diff = max(len(a.sequence), len(b.sequence)) - l
			dist_prefixes = hamming_distance(a.sequence[:l], b.sequence[:l])
			dist_suffixes = hamming_distance(a.sequence[-l:], b.sequence[-l:])
			score = l - 5 * min(dist_prefixes, dist_suffixes) - length_diff
			if score >= l - 19 and score > best_score:
				best_score = score
				best = b
		if best is not None:
			b_records.remove(best)
			pairs.append((a, best))
		else:
			only_a.append(a)
	only_b = b_records
	identical = [(a, b) for a, b in pairs if a.sequence == b.sequence]
	similar = [(a, b) for a, b in pairs if a.sequence != b.sequence]

	return only_a, only_b, identical, similar


def main(args):
	with FastaReader(args.a) as f:
		a_records = list(f)
	with FastaReader(args.b) as f:
		b_records = list(f)

	for records, path in ((a_records, args.a), (b_records, args.b)):
		dups = list(check_duplicate_names(records))
		if dups:
			print('Duplicate record names found in', path)
			for name in dups:
				print('-', name)

	for record, path in ((a_records, args.a), (b_records, args.b)):
		dups = list(check_exact_duplicate_sequences(records))
		if dups:
			print('Duplicate sequences found in', path)
			for name, name_orig in dups:
				print('-', name, 'is identical to earlier record', name_orig)

	only_a, only_b, identical, similar = pair_up(a_records, b_records)

	# One-line summary
	print('A:', len(a_records), 'records')
	print('B:', len(b_records), 'records')
	print('{} lost. {} gained. {} identical. {} similar.'.format(len(only_a), len(only_b), len(identical), len(similar)))
	print()
	print('## Only in A')
	for record in only_a:
		print('-', record.name)
	print()
	print('## Only in B')
	for record in only_b:
		print('+', record.name)
	print()
	print('## Identical')
	for a, b in identical:
		print('=', a.name, '--', b.name)
	print()
	print('## Similar')
	for a, b in similar:
		l = min(len(a.sequence), len(b.sequence))
		length_diff = max(len(a.sequence), len(b.sequence)) - l
		dist_prefixes = hamming_distance(a.sequence[:l], b.sequence[:l])
		dist_suffixes = hamming_distance(a.sequence[-l:], b.sequence[-l:])
		print('~', a.name, '--', b.name, 'length_diff', length_diff, 'hamming distance', min(dist_prefixes, dist_suffixes))

