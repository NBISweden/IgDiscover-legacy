"""
Group sequences by barcode and V/J assignment and print each group’s consensus
"""
from collections import Counter
import logging
from itertools import islice
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from sqt.ansicolor import bggreen
from sqt.align import multialign, consensus
from table import read_table_and_filter

logger = logging.getLogger(__name__)

def add_subcommand(subparsers):
	subparser = subparsers.add_parser('group', help=__doc__)
	subparser.set_defaults(func=group_command)
	subparser.add_argument('--groups-output', metavar='FILE', default=None,
		help='Write tab-separated table with groups to FILE')
	subparser.add_argument('--plot-sizes', metavar='FILE', default=None,
		help='Plot group sizes to FILE (.png or .pdf)')
	subparser.add_argument('--barcode-length', type=int, default=None)
	subparser.add_argument('--program', choices=('clustalo', 'muscle', 'muscle-fast', 'mafft'),
		default='muscle-fast',
		help='Program to use for computing the multiple alignment')
	subparser.add_argument('table', help='Table with parsed IgBLAST results.')
	return subparser


def write_groups(groups, path, limit=None):
	with open(path, 'w') as f:
		header = True
		for _, group in islice(groups, limit):
			group.to_csv(f, sep='\t', header=header, index=False)
			header = False
			print(file=f)


def colored_diff(seq, reference):
	"""
	Return a color-coded string of seq where the differences to refseq are highlighted.
	"""
	if pd.isnull(seq) or pd.isnull(reference):
		return seq
	colored = ''.join(bggreen(c1) if c1 != c2 else c1 for c1, c2 in zip(seq, reference))
	colored += bgreen(seq[len(reference):])
	return colored


def group_command(args):
	d = read_table_and_filter(args.table, log=True)
	barcode_length = args.barcode_length
	program = args.program
	assert barcode_length is not None
	lengths = Counter()
	d['barcode'] = [ s[:barcode_length] for s in d.sequence ]
	d['sequence'] = [ s[barcode_length:] for s in d.sequence ]

	# Group by barcode, V- and J gene assignment
	groups = d.groupby(('barcode', 'V_gene', 'J_gene', 'CDR3_nt'))
	sizes = groups['count'].sum()
	logger.info('Total no. of groups: %s', len(groups))
	for s in [1, 2, 3]:
		logger.info('Groups of size %s: %s', s, sum(1 for size in sizes if size == s))

	if args.groups_output:
		write_groups(groups, args.groups_output)
		logger.info('Groups written to %r', args.groups_output)

	# Plot group sizes if requested
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
		fig.tight_layout()
		fig.savefig(args.plot_sizes)

	# Compute consensus within each group
	n = 1
	for (barcode, v_gene, j_gene, cdr3_nt), group in groups:
		if len(group) == 1:
			cons = group['sequence'].iloc[0]
		else:
			sequences = dict()
			for i, sequence in enumerate(group['sequence']):
				sequences[str(i)] = sequence
			aligned = multialign(sequences, program=program)
			cons = consensus(aligned, threshold=0.65)
		size = group['count'].sum()
		print('>consensus{};barcode={};size={};\n{}'.format(n, barcode, size, cons))
		n += 1
		"""
		cons = consensus(aligned, keep_gaps=True, threshold=0.65)
		print(cons)
		for sequence in aligned.values():
			print(colored_diff(sequence, cons))
		print()
		"""
