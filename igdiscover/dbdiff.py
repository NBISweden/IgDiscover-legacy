"""
Compare two germline gene databases given as FASTA files
"""
import sys
import logging
from sqt import FastaReader
from sqt.align import hamming_distance
#from .utils import natural_sort_key

logger = logging.getLogger(__name__)

do_not_show_cpustats = 1


def add_arguments(parser):
	arg = parser.add_argument
	arg('a', help='FASTA file with expected sequences')
	arg('b', help='FASTA file with actual sequences')


RED = "\x1b[0;31m"
GREEN = "\x1b[0;32m"
RESET = "\x1b[0m"


def red(s):
	return RED + s + RESET


def green(s):
	return GREEN + s + RESET


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


def print_similar(a, b):
	l = min(len(a.sequence), len(b.sequence))
	length_diff = max(len(a.sequence), len(b.sequence)) - l
	dist_prefixes = hamming_distance(a.sequence[:l], b.sequence[:l])
	dist_suffixes = hamming_distance(a.sequence[-l:], b.sequence[-l:])
	if dist_prefixes <= dist_suffixes:
		a_prefix = ''
		b_prefix = ''
		a_common = a.sequence[:l]
		b_common = b.sequence[:l]
		a_suffix = a.sequence[l:]
		b_suffix = b.sequence[l:]
	else:
		a_prefix = a.sequence[:-l]
		b_prefix = b.sequence[:-l]
		a_common = a.sequence[-l:]
		b_common = b.sequence[-l:]
		a_suffix = ''
		b_suffix = ''

	def format_indel(a, b):
		if len(a) > len(b):
			assert len(b) == 0
			return red('{-' + a + '}')
		elif len(b) > len(a):
			assert len(a) == 0
			return green('{+' + b + '}')
		else:
			return ''

	s = format_indel(a_prefix, b_prefix)
	edits = []
	for i, (ac, bc) in enumerate(zip(a_common, b_common)):
		if ac != bc:
			edits.append('{' + red(ac) + ' → ' + green(bc) + '}')
		else:
			edits.append(ac)
	s += ''.join(edits)

	s += format_indel(a_suffix, b_suffix)

	print('~', a.name, '--', b.name)
	print(s)
	print()


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
		print_similar(a, b)
