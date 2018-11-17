"""
Some utility functions that work on sequences and lists of sequences.
"""
import os
import sys
import re
import random
import hashlib
import resource
from collections import OrderedDict, Counter, defaultdict
from itertools import groupby
from typing import List, Tuple

import numpy as np
from sqt.align import edit_distance, multialign, consensus, globalalign
from sqt.dna import GENETIC_CODE, nt_to_aa as _nt_to_aa
from sqt import SequenceReader
from cutadapt.align import Aligner


def downsampled(population, size):
	"""
	Return a random subsample of the population.

	Uses reservoir sampling. See https://en.wikipedia.org/wiki/Reservoir_sampling
	"""
	sample = population[:size]
	for index in range(size, len(population)):
		r = random.randint(0, index)
		if r < size:
			sample[r] = population[index]
	return sample


def distances(sequences, band=0.2):
	"""
	Compute all pairwise edit distances and return a square matrix.

	Entry [i,j] in the matrix is the edit distance between sequences[i]
	and sequences[j].
	"""
	# Pre-compute distances between unique sequences
	unique_sequences = list(set(sequences))
	unique_distances = dict()  # maps (seq1, seq2) tuples to edit distance
	maxdiff = max((int(len(s) * band) for s in sequences), default=0)  # TODO double-check this
	for i, s in enumerate(unique_sequences):
		for j, t in enumerate(unique_sequences):
			if i < j:
				dist = min(maxdiff+1, edit_distance(s, t, maxdiff=maxdiff))
				unique_distances[(t, s)] = dist
				unique_distances[(s, t)] = dist

	# Fill the result matrix
	m = np.zeros((len(sequences), len(sequences)), dtype=float)
	for i, s in enumerate(sequences):
		for j, t in enumerate(sequences):
			if i < j:
				d = 0 if s == t else unique_distances[(s, t)]
				m[j, i] = m[i, j] = d
	return m


# copied from sqt and modified
def consensus(aligned, threshold=0.7, ambiguous='N', keep_gaps=False):
	"""
	Compute a consensus from multialign() output, allowing degraded sequences
	in the 3' end.

	aligned -- a dict mapping names to sequences or a list of sequences
	keep_gaps -- whether the returned sequence contains gaps (-)
	"""
	n = len(aligned)
	result = []
	if hasattr(aligned, 'values'):
		sequences = aligned.values()
	else:
		sequences = aligned

	active = int(len(aligned) * 0.05)
	for i, chars in enumerate(reversed(list(zip(*sequences)))):
		counter = Counter(chars)
		active = max(n - counter['-'], active)
		assert counter['-'] >= n - active
		counter['-'] -= n - active
		char, freq = counter.most_common(1)[0]
		if i >= 10:  # TODO hard-coded
			active = n
		if freq / active >= threshold:
			if keep_gaps or char != '-':
				result.append(char)
		else:
			result.append(ambiguous)
	return ''.join(result[::-1])


def iterative_consensus(sequences, program='muscle-medium', threshold=0.6,
		subsample_size=200, maximum_subsample_size=1600):
	"""
	Compute a consensus sequence of the given sequences, but do not use all
	sequences if there are many: First, try to compute the consensus from a
	small subsample. If there are 'N' bases, increase the subsample size and
	repeat until either there are no more 'N' bases, all available sequences
	have been used or maximum_subsample_size is reached.
	"""
	while True:
		sample = downsampled(sequences, subsample_size)
		aligned = multialign(OrderedDict(enumerate(sample)), program=program)

		cons = consensus(aligned, threshold=threshold).strip('N')
		if 'N' not in cons:
			# This consensus is good enough
			break
		if len(sequences) <= subsample_size:
			# We have already used all the sequences that are available
			break
		subsample_size *= 2
		if subsample_size > maximum_subsample_size:
			break
	return cons


def sequence_hash(s, digits=4):
	"""
	For a string, return a 'fingerprint' that looks like 'S1234' (the character
	'S' is fixed). The idea is that this allows one to quickly see whether two
	sequences are not identical.
	"""
	h = int(hashlib.md5(s.encode()).hexdigest()[-4:], base=16)
	return 'S' + str(h % 10**digits).rjust(digits, '0')


def unique_name(name, sequence):
	"""Create a unique name based on the current name and the sequence

	The returned name looks like name_S1234. If the current name contains
	already a S_.... suffix, it is removed before the new suffix is appended.

	name -- current name
	"""
	return '{}_{}'.format(name.rsplit('_S', 1)[0], sequence_hash(sequence))


class SerialPool:
	"""
	An alternative to multiprocessing.Pool that runs things in serial for
	easier debugging
	"""
	def __init__(self, *args, **kwargs):
		pass

	def __enter__(self):
		return self

	def __exit__(self, *args):
		pass

	def imap(self, func, iterable, chunksize):
		for i in iterable:
			yield func(i)


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
	"""
	Use this function as sorting key to sort in 'natural' order.

	>>> names = ['file10', 'file1.5', 'file1']
	>>> sorted(names, key=natural_sort_key)
	['file1', 'file1.5', 'file10']

	Source: http://stackoverflow.com/a/16090640/715090
	"""
	return [int(text) if text.isdigit() else text.lower()
		for text in re.split(_nsre, s)]


class UniqueNamer:
	"""
	Assign unique names by appending letters to already seen names.
	"""
	def __init__(self):
		self._names = set()

	def __call__(self, name):
		ext = 'A'
		new_name = name
		while new_name in self._names:
			if ext == '[':
				raise ValueError('Too many duplicate names')
			new_name = name + ext
			ext = chr(ord(ext) + 1)
		self._names.add(new_name)
		return new_name


class Merger:
	"""
	Merge similar items into one. To specify what "similar" means, implement
	the merged() method in a subclass.
	"""
	def __init__(self):
		self._items = []

	def add(self, item):
		# This method could possibly be made simpler if the graph structure
		# was made explicit.

		items = []
		for existing_item in self._items:
			m = self.merged(existing_item, item)
			if m is None:
				items.append(existing_item)
			else:
				item = m
		items.append(item)
		self._items = items

	def extend(self, iterable):
		for i in iterable:
			self.add(i)

	def __iter__(self):
		if self._items and hasattr(self._items, 'name'):
			yield from sorted(self._items, key=lambda x: x.name)
		else:
			yield from self._items

	def __len__(self):
		return len(self._items)

	def merged(self, existing_item, item):
		"""
		If existing_item and item can be returned, this method must return
		a new item that represents a merged version of both. If they cannot
		be merged, it must return None.
		"""
		raise NotImplementedError("not implemented")


def relative_symlink(src, dst, force=False):
	"""
	Create a symbolic link in any directory.

	force -- if True, then overwrite an existing file/symlink
	"""
	if force:
		try:
			os.remove(dst)
		except FileNotFoundError:
			pass
	target = os.path.relpath(os.path.abspath(src), start=os.path.dirname(dst))
	os.symlink(target, dst)


def nt_to_aa(s):
	"""Translate a nucleotide sequence to an amino acid sequence"""
	try:
		# try fast version first
		return _nt_to_aa(s)
	except ValueError:
		# failure because there was an unknown nucleotide
		return ''.join(GENETIC_CODE.get(s[i:i+3], '*') for i in range(0, len(s), 3))


def has_stop(sequence):
	"""
	Return a boolean indicating whether the sequence has an internal stop codon.
	An incomplete codon at the end is allowed.

	>>> has_stop('GGG')
	False
	>>> has_stop('TAA')
	True
	>>> has_stop('GGGAC')
	False
	"""
	s = sequence[:len(sequence) // 3 * 3]
	return '*' in nt_to_aa(s)


def plural_s(n):
	return 's' if n != 1 else ''


class FastaValidationError(Exception):
	pass


def validate_fasta(path):
	"""
	Ensure that the FASTA file is suitable for use with makeblastdb.
	Raise a FastaValidationError if any of the following are true:

	- a record is empty
	- a record name occurs more than once
	- a sequence occurs more than once
	"""
	with SequenceReader(path) as sr:
		records = list(sr)

	names = set()
	sequences = dict()
	for r in records:
		if len(r.sequence) == 0:
			raise FastaValidationError("Record {!r} is empty".format(r.name))
		if r.name in names:
			raise FastaValidationError("Record name {!r} occurs more than once".format(r.name))
		s = r.sequence.upper()
		if s in sequences:
			raise FastaValidationError("Records {!r} and {!r} contain the same sequence".format(
				r.name, sequences[s]))
		sequences[s] = r.name
		names.add(r.name)


def find_overlap(s, t, min_overlap=1):
	"""
	Detect if s and t overlap.

	Returns:

	None if no overlap was detected.
	0 if s is a prefix of t or t is a prefix of s.
	Positive int gives index where t starts within s.
	Negative int gives -index where s starts within t.

	>>> find_overlap('ABCDE', 'CDE')
	2
	>>> find_overlap('CDE', 'ABCDEFG')
	-2
	>>> find_overlap('ABC', 'X') is None
	True
	"""
	aligner = Aligner(s, max_error_rate=0)
	aligner.min_overlap = min_overlap
	result = aligner.locate(t)
	if result is None:
		return None
	s_start, _, t_start, _, _, _ = result
	return s_start - t_start


def merge_overlapping(s, t):
	"""
	Return merged sequences or None if they do not overlap
	The minimum overlap is 50% of the length of the shorter sequence.
	"""
	i = find_overlap(s, t, min_overlap=max(1, min(len(s), len(t)) // 2))
	if i is None:
		return None

	if i >= 0:
		# positive: index of t in s
		if i + len(t) < len(s):
			# t is in s
			return s
		return s[:i] + t
	if -i + len(s) < len(t):
		# s is in t
		return t
	return t[:-i] + s


def get_cpu_time():
	"""Return CPU time used by process and children"""
	if sys.platform != 'linux':
		return None
	rs = resource.getrusage(resource.RUSAGE_SELF)
	rc = resource.getrusage(resource.RUSAGE_CHILDREN)
	return rs.ru_utime + rs.ru_stime + rc.ru_utime + rc.ru_stime


def slice_arg(s):
	"""
	Parse a string that describes a slice with start and end.

	>>> slice_arg('2:-3')
	slice(2, -3, None)

	>> slice_arg(':-3')
	slice(None, -3, None)

	>> slice_arg('2:')
	slice(2, None, None)
	"""

	start, end = s.split(':')
	start = None if start == '' else int(start)
	end = None if end == '' else int(end)
	return slice(start, end)


def is_same_gene(name1: str, name2: str):
	"""
	Compare gene names to find out whether they are alleles of each other.
	Both names must have a '*' in them
	"""
	return '*' in name1 and '*' in name2 and name1.split('*')[0] == name2.split('*')[0]


def describe_nt_change(s: str, t: str):
	"""
	Describe changes between two nucleotide sequences

	>>> describe_nt_change('AAA', 'AGA')
	'2A>G'
	>>> describe_nt_change('AAGG', 'AATTGG')
	'2_3insTT'
	>>> describe_nt_change('AATTGG', 'AAGG')
	'3_4delTT'
	>>> describe_nt_change('AATTGGCG', 'AAGGTG')
	'3_4delTT; 7C>T'
	"""
	row1, row2, start1, stop1, start2, stop2, errors = \
		globalalign(s.encode('ascii'), t.encode('ascii'), flags=0, match=0)
	row1 = row1.decode('ascii').replace('\0', '-')
	row2 = row2.decode('ascii').replace('\0', '-')
	changes = []

	def grouper(c):
		c1, c2 = c
		if c1 == c2:
			return 'MATCH'
		elif c1 == '-':
			return 'INS'
		elif c2 == '-':
			return 'DEL'
		else:
			return 'SUBST'

	index = 1
	for event, group in groupby(zip(row1, row2), grouper):
		if event == 'MATCH':
			index += len(list(group))
		elif event == 'SUBST':
			# ungroup
			for c1, c2 in group:
				change = '{}{}>{}'.format(index, c1, c2)
				changes.append(change)
				index += 1
		elif event == 'INS':
			inserted = ''.join(c[1] for c in group)
			change = '{}_{}ins{}'.format(index-1, index, inserted)
			changes.append(change)
		elif event == 'DEL':
			deleted = ''.join(c[0] for c in group)
			change = '{}_{}del{}'.format(index, index + len(deleted) - 1, deleted)
			changes.append(change)
			index += len(deleted)
	return '; '.join(changes)


class ChimeraFinder:
	def __init__(self, sequences: List[str], min_length: int=10):
		self._sequences = sequences
		self._min_length = min_length
		self._build_index()

	def _build_index(self):
		min_length = self._min_length
		# Create two dictionaries that map all prefixes and suffixes to indexes of all
		# sequences they occur in
		prefixes = defaultdict(list)
		suffixes = defaultdict(list)

		for i, sequence in enumerate(self._sequences):
			if len(sequence) < min_length:
				continue
			for stop in range(min_length, len(sequence) + 1):
				prefix = sequence[:stop]
				prefixes[prefix].append(i)
			for start in range(0, len(sequence) + 1 - min_length):
				suffix = sequence[start:]
				suffixes[suffix].append(i)
		self._prefixes = prefixes
		self._suffixes = suffixes

	def find_exact(self, query: str):
		"""
		Find out whether the query string can be explained as a concatenation of
		a prefix of one of the strings plus a suffix of one of the strings in
		the given list of sequences. Both the prefix and the suffix must have a
		length of at least min_length.

		If the answer is yes, a tuple (prefix_length, prefix_indices, suffix_indices)
		is returned where prefix_length is the length of the query prefix,
		prefix_indices is a list of int indices into the sequences of all possible
		prefix sequences that match, and suffix_indices is the same for the suffix.

		The prefix_length returned is the first that yields a result. More are
		possible.

		If the answer is no, None is returned.
		"""
		min_length = self._min_length
		for split_index in range(min_length, len(query) + 1 - min_length):
			prefix = query[:split_index]
			suffix = query[split_index:]
			if prefix in self._prefixes and suffix in self._suffixes:
				return (split_index, self._prefixes[prefix], self._suffixes[suffix])

		return None


def available_cpu_count():
	return len(os.sched_getaffinity(0))
