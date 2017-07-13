#!/usr/bin/env python3
"""
Merge overlapping sequences
"""
from igdiscover.utils import merge_overlapping, Merger
import sys
from argparse import ArgumentParser
from sqt.io.fasta import FastaReader


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
