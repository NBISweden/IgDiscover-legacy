"""
Group sequences that share a barcode (molecular identifier, MID)

Since the same barcode can sometimes be used by different sequences, a 'pseudo
CDR3' sequence (defined to be 80 to 60 bases upstream of the 3' end) is used
to further distinguish sequences.

The barcode is assumed to be in the 5' end of the sequence.

Use --trim-g to remove initial runs of G (at 5' end, artifact from RACE
protocol).

For all the found groups, one sequence is output to standard output (in FASTA
format). Which sequence that is depends on the group size:
- If the group consists of a single sequence, that sequence is output
- If the group consists of two sequences, one of them is picked randomly
- If the group has at least three sequences, a consensus is computed. The
  consensus is output if it contains no ambiguous bases. Otherwise, also here a
  random sequence is chosen.
"""
"""
Things to keep in mind
- There are some indels in homopolymers (probably PCR problem)
- Different lengths of the initial G run cannot be used to distinguish sequences
  since they can come from polymerase problems in homopolymers.
- There are also regular sequencing errors in the initial run of G nucleotides.
- Some paired reads aren’t correctly merged into single reads. They end up being
  too long.
- When grouping by barcode and pseudo CDR3, sequence lengths vary within groups.
  However, this affects only ~1% of sequences, so it is not necessary to compute
  a multiple alignment. Just taking the consensus will drown the incorrect
  sequences, at least if the group size is large.
- It does not hurt to reduce the minimimum number of sequences per group for
  taking a consensus to 2, but it also does not help much (results in 0.5% more
  sequences). (The consensus is only successful if both sequences are identical.)
  However, since this is also a simple way to deal with exact duplicates, we do
  it anyway and can then skip the separate duplicate removal step (VSEARCH).

TODO
- Use pandas.DataFrame
"""

import csv
import sys
import logging
from collections import Counter, defaultdict
from contextlib import ExitStack
from itertools import islice

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from sqt.align import consensus
from sqt import SequenceReader

from .species import CDR3_REGEX

# minimum number of sequences needed for attempting to compute a consensus
MIN_CONSENSUS_SEQUENCES = 3


logger = logging.getLogger(__name__)

def add_arguments(parser):
	parser.add_argument('--groups-output', metavar='FILE', default=None,
		help='Write tab-separated table with groups to FILE')
	parser.add_argument('--plot-sizes', metavar='FILE', default=None,
		help='Plot group sizes to FILE (.png or .pdf)')
	parser.add_argument('--real-cdr3', action='store_true', default=False,
		help='Use the real CDR3 (detected with regex) instead of pseudo CDR3')
	parser.add_argument('--limit', default=None, type=int, metavar='N',
		help='Limit processing to the first N reads')
	parser.add_argument('--trim-g', action='store_true', default=False,
		help="Trim 'G' nucleotides at 5' end")
	parser.add_argument('--minimum-length', '-l', type=int, default=0,
		help='Minimum sequence length')
	parser.add_argument('--barcode-length', '-b', type=int, default=12,
		help="Length of 5' barcode (default: %(default)s")
	parser.add_argument('fastx', metavar='FASTA/FASTQ',
		help='FASTA or FASTQ file (can be gzip-compressed) with sequences')


class Graph:
	"""Graph that can find connected components"""
	def __init__(self, nodes):
		self._nodes = { node: [] for node in nodes }

	def add_edge(self, node1, node2):
		self._nodes[node1].append(node2)
		self._nodes[node2].append(node1)

	def connected_components(self):
		"""Return a list of connected components."""
		visited = set()
		components = []
		for node, neighbors in self._nodes.items():
			if node in visited:
				continue
			# Start a new component
			to_visit = [node]
			component = []
			while to_visit:
				n = to_visit.pop()
				if n in visited:
					continue
				visited.add(n)
				component.append(n)
				for neighbor in self._nodes[n]:
					if neighbor not in visited:
						to_visit.append(neighbor)
			components.append(component)
		return components


def hamming_neighbors(s):
	"""Return sequences that are at hamming distance 1 and return also s itself"""
	for i in range(len(s)):
		for c in 'ACGT':
			if s[i] != c:
				yield s[:i] + c + s[i+1:]
	yield s


def cluster_sequences(sequences):
	"""
	Single-linkage clustering. Two sequences are linked if
	- their (pseudo-) CDR3 sequences have a hamming distance of at most 1
	- and their lengths differs by at most 2.
	"""
	graph = Graph(sequences)
	cdr3_seqs = defaultdict(list)
	for s in sequences:
		cdr3_seqs[s.cdr3].append(s)
	for sequence in sequences:
		for neighbor in hamming_neighbors(sequence.cdr3):
			neighbor_seqs = cdr3_seqs.get(neighbor, [])
			for neighbor_seq in neighbor_seqs:
				if abs(len(sequence) - len(neighbor_seq)) <= 2:
					graph.add_edge(sequence, neighbor_seq)
	components = graph.connected_components()
	assert sum(len(component) for component in components) == len(sequences)
	return components


GROUPS_HEADER = ['barcode', 'cdr3', 'name', 'sequence']

def write_group(csvfile, barcode, sequences):
	for sequence in sequences:
		row = [barcode, sequence.cdr3, sequence.name.split(maxsplit=1)[0], sequence.sequence]
		csvfile.writerow(row)
	csvfile.writerow([])


def main(args):
	if args.barcode_length < 1:
		sys.exit("Barcode length must be positive")

	# Map barcodes to lists of sequences
	barcodes = defaultdict(list)
	n = 0
	too_short = 0
	cdr3s = set()
	regex_fail = 0
	with SequenceReader(args.fastx) as f:
		for record in islice(f, 0, args.limit):
			if len(record) < args.minimum_length:
				too_short += 1
				continue
			barcode = record.sequence[:args.barcode_length]
			unbarcoded = record[args.barcode_length:]
			if args.trim_g:
				# The RACE protocol leads to a run of non-template Gs in the beginning
				# of the sequence, after the barcode.
				unbarcoded.sequence = unbarcoded.sequence.lstrip('G')
				if unbarcoded.qualities:
					unbarcoded.qualities = unbarcoded.qualities[-len(unbarcoded.sequence):]

			if args.real_cdr3:
				match = CDR3_REGEX['VH'].search(unbarcoded.sequence)
				if match:
					cdr3 = match.group('cdr3')
				else:
					regex_fail += 1
					continue
			else:
				cdr3 = unbarcoded.sequence[-80:-60]
			unbarcoded.cdr3 = cdr3  # TODO slight abuse of Sequence objects
			barcodes[barcode].append(unbarcoded)
			cdr3s.add(cdr3)
			n += 1

	logger.info('%s sequences in input', n + too_short + regex_fail)
	logger.info('%s sequences long enough', n + regex_fail)
	if args.real_cdr3:
		logger.info('Using the real CDR3')
		logger.info('%s times (%.2f%%), the CDR3 regex matched', n, n / (n + regex_fail) * 100)
	else:
		logger.info('Using the pseudo CDR3')
	logger.info('%s unique CDR3s', len(cdr3s))
	logger.info('%s unique barcodes', len(barcodes))
	barcode_singletons = sum(1 for seqs in barcodes.values() if len(seqs) == 1)
	logger.info('%s barcodes used by only a single sequence (singletons)', barcode_singletons)

	with ExitStack() as stack:
		if args.groups_output:
			group_out = csv.writer(stack.enter_context(
				open(args.groups_output, 'w')), delimiter='\t', lineterminator='\n')
			group_out.writerow(GROUPS_HEADER)
		else:
			group_out = None
		too_few = 0
		n_clusters = 0
		n_singletons = 0
		n_consensus = 0
		n_ambiguous = 0
		sizes = []
		for barcode, sequences in barcodes.items():
			if len(sequences) != len(set(sequences)):
				logger.error('Duplicate sequence records detected')
				sys.exit(1)
			clusters = cluster_sequences(sequences)
			n_clusters += len(clusters)
			for cluster in clusters:
				sizes.append(len(cluster))
				if group_out:
					write_group(group_out, barcode, cluster)
				if len(cluster) == 1:
					n_singletons += 1
				if len(cluster) < MIN_CONSENSUS_SEQUENCES:
					too_few += 1
					sequence = cluster[0].sequence
					name = cluster[0].name
					cdr3 = cluster[0].cdr3
				else:
					cons = consensus({s.name: s.sequence for s in cluster}, threshold=0.501)
					if 'N' in cons:
						# Pick the first sequence as the output sequence
						sequence = cluster[0].sequence
						name = cluster[0].name
						cdr3 = cluster[0].cdr3
						n_ambiguous += 1
					else:
						sequence = cons
						n_consensus += 1
						cdr3 = Counter(cl.cdr3 for cl in cluster).most_common(1)[0][0]
						name = 'consensus{}'.format(n_consensus)
				print('>{} barcode={};cdr3={};size={};\n{}'.format(name, barcode,
					cdr3, len(cluster), sequence))

	logger.info('%d clusters (%d singletons)', n_clusters, n_singletons)
	logger.info('%d consensus sequences computed (from groups that had at least %d sequences)', n_consensus+n_ambiguous, MIN_CONSENSUS_SEQUENCES)
	logger.info('%d of those had no ambiguous bases', n_consensus)

	if args.plot_sizes:
		matplotlib.rcParams.update({'font.size': 14})
		fig = plt.figure()
		ax = fig.gca()
		v, _, _ = ax.hist(sizes, bins=100)
		ax.set_ylim(0, v[1:].max()*1.1)
		ax.set_xlabel('Group size')
		ax.set_ylabel('Read frequency')
		ax.set_title('Histogram of group sizes (>1)')
		ax.grid(axis='x')
		ax.tick_params(direction="outward", top=False, right=False)
		fig.set_tight_layout(True)
		fig.savefig(args.plot_sizes)
		logger.info('Plotted group sizes to %r', args.plot_sizes)

	if args.groups_output:
		logger.info('Groups written to %r', args.groups_output)
