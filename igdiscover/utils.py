"""
Some utility functions that work on sequences and lists of sequences.
"""
import os
import sys
import re
import random
import hashlib
import resource
from collections import OrderedDict

import numpy as np
from sqt.align import edit_distance, multialign, consensus
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
	m = np.zeros((len(sequences), len(sequences)), dtype=float)
	maxdiff = max((int(len(s) * band) for s in sequences), default=0)
	for i, s in enumerate(sequences):
		for j, t in enumerate(sequences):
			if i < j:
				m[j,i] = m[i,j] = min(maxdiff+1, edit_distance(s, t, maxdiff=maxdiff))
	return m


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
