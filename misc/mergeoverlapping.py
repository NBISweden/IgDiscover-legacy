#!/usr/bin/env python3
"""
Merge overlapping sequences
"""
import sys
from argparse import ArgumentParser
from sqt.io.fasta import FastaReader
from cutadapt.align import Aligner


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


class SequenceInfo:
	__slots__ = ('name', 'sequence', 'count')

	def __init__(self, name, sequence):
		self.name = name
		self.sequence = sequence
		# self.count = count


class OverlappingSequenceMerger(Merger):
	"""
	Merge sequences that overlap
	"""
	def merged(self, s, t):
		"""
		Merge two sequences if they overlap. Return None if they should not be merged.
		"""
		m = merge_overlapping(s.sequence, t.sequence)
		if m is not None:
			# overlap found
			return SequenceInfo(s.name + ';' + t.name, m)

		# no overlap found
		return None


def add_arguments(parser):
	arg = parser.add_argument
	arg('fasta', help='File with sequences')


def main(args):
	merger = OverlappingSequenceMerger()
	with FastaReader(args.fasta) as f:
		for record in f:
			merger.add(SequenceInfo(record.name, record.sequence))

	for info in merger:
		print('>{}\n{}'.format(info.name, info.sequence))


if __name__ == '__main__':
	parser = ArgumentParser(description=__doc__)
	add_arguments(parser)
	args = parser.parse_args()
	main(args)
