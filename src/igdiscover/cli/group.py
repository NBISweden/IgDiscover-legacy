"""
Group sequences that share a unique molecular identifier (UMI, barcode)

Since the same UMI can sometimes be used by different sequences, the CDR3
sequence can further be used to distinguish sequences. You can choose between
using either a 'pseudo CDR3' sequence, which encompasses by default bases 80
to 61 counted from the 3' end. Or you can use the real CDR3 detected with a
regular expression.

If grouping by CDR3s is enabled, sequences with identical barcode and CDR3
must additionally have a similar length. If the length differs by more than
2 bp, they are put into different groups.

The barcode can be in the 5' end or the 3' end of the sequence.

Use --trim-g to remove initial runs of G at the 5' end (artifact from RACE protocol).
These are removed after the barcode is removed.

For all the found groups, one sequence is output to standard output (in FASTA
format). Which sequence that is depends on the group size:
- If the group consists of a single sequence, that sequence is output
- If the group consists of two sequences, one of them is picked randomly
- If the group has at least three sequences, a consensus is computed. The
  consensus is output if it contains no ambiguous bases. Otherwise, also here a
  random sequence is chosen.
"""
# NOTES
#
# - Different lengths of the initial G run cannot be used to distinguish sequences
#  since they can come from polymerase problems in homopolymers.
# - There are some indels in homopolymers (probably PCR problem)
# - There are also regular sequencing errors in the initial run of G nucleotides.
# - Some paired reads aren’t correctly merged into single reads. They end up being
#   too long.
# - When grouping by barcode and pseudo CDR3, sequence lengths vary within groups.
#   However, this affects only ~1% of sequences, so it is not necessary to compute
#   a multiple alignment. Just taking the consensus will drown the incorrect
#   sequences, at least if the group size is large.
# - It does not hurt to reduce the minimimum number of sequences per group for
#   taking a consensus to 2, but it also does not help much (results in 0.5% more
#   sequences). (The consensus is only successful if both sequences are identical.)
#   However, since this is also a simple way to deal with exact duplicates, we do
#   it anyway and can then skip the separate duplicate removal step (VSEARCH).

# TODO
# - Use pandas.DataFrame

import csv
import sys
import logging
from collections import Counter, defaultdict
from contextlib import ExitStack
from itertools import islice
import json
from typing import List, Iterable

import dnaio
from xopen import xopen

from ..align import consensus
from ..species import find_cdr3
from ..cluster import Graph
from ..utils import slice_arg

# minimum number of sequences needed for attempting to compute a consensus
MIN_CONSENSUS_SEQUENCES = 3


logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--real-cdr3', action='store_true', default=False,
        help='In addition to barcode, group sequences by real CDR3 (detected with regex).')
    group.add_argument('--pseudo-cdr3', nargs='?', default=None,
        type=slice_arg, const=slice(-80, -60), metavar='START:END',
        help='In addition to barcode, group sequences by pseudo CDR3. '
            'If START:END is omitted, use -80:-60.')
    arg('--groups-output', metavar='FILE', default=None,
        help='Write tab-separated table with groups to FILE')
    arg('--plot-sizes', metavar='FILE', default=None,
        help='Plot group sizes to FILE (.png or .pdf)')
    arg('--limit', default=None, type=int, metavar='N',
        help='Limit processing to the first N reads')
    arg('--trim-g', action='store_true', default=False,
        help="Trim 'G' nucleotides at 5' end")
    arg('--minimum-length', '-l', type=int, default=0,
        help='Minimum sequence length')
    arg('--barcode-length', '-b', type=int, default=12,
        help="Length of barcode. Positive for 5' barcode, negative for 3' barcode. Default: %(default)s")
    arg('--json', metavar="FILE", help="Write statistics to FILE")
    arg('fastx', metavar='FASTA/FASTQ',
        help='FASTA or FASTQ file (can be gzip-compressed) with sequences')


def hamming_neighbors(s: str) -> Iterable[str]:
    """Return sequences that are at hamming distance 1 and return also s itself"""
    for i in range(len(s)):
        for c in 'ACGT':
            if s[i] != c:
                yield s[:i] + c + s[i+1:]
    yield s


def cluster_by_cdr3(records, pseudo_cdr3: slice, real_cdr3: bool) -> List[List]:
    """
    Single-linkage clustering of sequences based on sequence length and CDR3 similarity.

    Two sequences are linked if
    - their (pseudo-) CDR3 sequences have a hamming distance of at most 1
    - and their full lengths (not the CDR3 length) differs by at most 2.

    Sequences for which no CDR3 could be detected are filtered out.
    """
    assert pseudo_cdr3 or real_cdr3

    records_cdr3s = []
    for record in records:
        if real_cdr3:
            match = find_cdr3(record.sequence, locus="IGH")
            if match:
                cdr3 = record.sequence[match[0]:match[1]]
            else:
                continue
        else:
            assert pseudo_cdr3
            cdr3 = record.sequence[pseudo_cdr3]
        records_cdr3s.append((record, cdr3))
    del records

    # logger.info('%s times (%.2f%%), the CDR3 regex matched', n, n / (n + regex_fail) * 100)

    if len(records_cdr3s) == 1:  # TODO check if this helps
        return [records_cdr3s]

    cdr3s = [rc[1] for rc in records_cdr3s]

    # Cluster by CDR3 sequences
    unique_cdr3s = set(cdr3s)
    sorted_cdr3s = sorted(unique_cdr3s)  # For reproducibility
    graph = Graph(sorted_cdr3s)
    for cdr3 in sorted_cdr3s:
        for neighbor in hamming_neighbors(cdr3):
            if neighbor in unique_cdr3s:
                graph.add_edge(cdr3, neighbor)
    cdr3_clusters = graph.connected_components()

    # Map each CDR3 sequence to a list of records with that CDR3
    cdr3_records = defaultdict(list)
    for record, cdr3 in records_cdr3s:
        cdr3_records[cdr3].append(record)

    # To compute the final clusters that also take length into account,
    # split the CDR3 clusters where necessary
    clusters = []
    for cdr3_cluster in cdr3_clusters:
        # Gather all records in this CDR3 cluster
        cluster_records = []
        for cdr3 in cdr3_cluster:
            for record in cdr3_records[cdr3]:
                cluster_records.append((record, cdr3))

        # For "clustering" by sequence length, sort by sequence length
        # and look at adjacent records: If their length differs too much,
        # the cluster needs to be split
        cluster_records.sort(key=lambda rc: len(rc[0].sequence))
        cluster = []
        prev_length = None
        for record, cdr3 in cluster_records:
            l = len(record.sequence)
            if prev_length is not None and l > prev_length + 2:
                # Start a new component
                clusters.append(cluster)
                cluster = []
            cluster.append((record, cdr3))
            prev_length = l
        if cluster:
            clusters.append(cluster)

    assert sum(len(component) for component in clusters) == len(records_cdr3s)
    assert all(clusters)  # Components must be non-empty
    return clusters


GROUPS_HEADER = ['barcode', 'cdr3', 'name', 'sequence']


def write_group(csvfile, barcode, sequences, write_cdr3: bool):
    for sequence, cdr3 in sequences:
        row = [barcode, sequence.name.split(maxsplit=1)[0], sequence.sequence]
        if write_cdr3:
            row[1:1] = [cdr3]
        csvfile.writerow(row)
    csvfile.writerow([])


def collect_barcode_groups(fastx, barcode_length, trim_g, limit, minimum_length):
    """
    fastx -- path to FASTA or FASTQ input

    """
    cdr3s = set()
    # Map barcodes to lists of sequences
    barcodes = defaultdict(list)
    n = 0
    too_short = 0
    with dnaio.open(fastx) as f:
        for record in islice(f, 0, limit):
            if len(record) < minimum_length:
                too_short += 1
                continue
            if barcode_length > 0:
                barcode = record.sequence[:barcode_length]
                unbarcoded = record[barcode_length:]
            else:
                barcode = record.sequence[barcode_length:]
                unbarcoded = record[:barcode_length]

            if trim_g:
                # The RACE protocol leads to a run of non-template Gs in the beginning
                # of the sequence, after the barcode.
                unbarcoded.sequence = unbarcoded.sequence.lstrip('G')
                if unbarcoded.qualities:
                    unbarcoded.qualities = unbarcoded.qualities[-len(unbarcoded.sequence):]

            barcodes[barcode].append(unbarcoded)
            n += 1

    logger.info('%s sequences in input', n + too_short)
    logger.info('%s sequences long enough', n)

    return barcodes


def plot_sizes(sizes, path):
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    import matplotlib
    import seaborn as sns
    sns.set()
    fig = Figure()
    matplotlib.rcParams.update({'font.size': 14})
    FigureCanvas(fig)
    ax = fig.add_subplot(111)
    v, _, _ = ax.hist(sizes, bins=100)
    ax.set_ylim(0, v[1:].max() * 1.1)
    ax.set_xlabel('Group size')
    ax.set_ylabel('Read frequency')
    ax.set_title('Histogram of group sizes (>1)')
    ax.grid(axis='x')
    ax.tick_params(direction='out', top=False, right=False)
    fig.set_tight_layout(True)
    fig.savefig(path)
    logger.info('Plotted group sizes to %r', path)


def main(args):
    if args.barcode_length == 0:
        sys.exit("Barcode length must be non-zero")

    barcodes = collect_barcode_groups(args.fastx, args.barcode_length, args.trim_g,
        args.limit, args.minimum_length)

    logger.info('%s unique barcodes', len(barcodes))
    barcode_singletons = sum(1 for seqs in barcodes.values() if len(seqs) == 1)
    logger.info('%s barcodes used by only a single sequence (singletons)', barcode_singletons)

    group_by_cdr3 = args.pseudo_cdr3 or args.real_cdr3

    with ExitStack() as stack:
        if args.groups_output:
            group_out = csv.writer(stack.enter_context(
                xopen(args.groups_output, 'w')), delimiter='\t', lineterminator='\n')
            group_out.writerow(GROUPS_HEADER)
        else:
            group_out = None
        too_few = 0
        n_clusters = 0
        n_singletons = 0
        n_consensus = 0
        n_ambiguous = 0
        sizes = []
        for barcode in sorted(barcodes):
            sequences = barcodes[barcode]
            if len(sequences) != len(set(sequence.name for sequence in sequences)):
                logger.error('Duplicate sequence records detected')
                sys.exit(1)
            if group_by_cdr3:
                clusters = cluster_by_cdr3(sequences, args.pseudo_cdr3, args.real_cdr3)
            else:
                # TODO it would be useful to do the length clustering that cluster_by_cdr3() does
                clusters = [[(record, None) for record in sequences]]
            n_clusters += len(clusters)
            for cluster in clusters:

                sizes.append(len(cluster))
                if group_out:
                    write_group(group_out, barcode, cluster, write_cdr3=group_by_cdr3)
                if len(cluster) == 1:
                    n_singletons += 1
                if len(cluster) < MIN_CONSENSUS_SEQUENCES:

                    sequence = cluster[0][0].sequence
                    name = cluster[0][0].name
                    cdr3 = cluster[0][1]
                    too_few += 1
                else:
                    cons = consensus({record.name: record.sequence for record, cdr3 in cluster}, threshold=0.501)
                    if 'N' in cons:
                        # Pick the first sequence as the output sequence
                        sequence = cluster[0][0].sequence
                        name = cluster[0][0].name
                        cdr3 = cluster[0][1]
                        n_ambiguous += 1
                    else:
                        sequence = cons
                        n_consensus += 1
                        name = 'consensus{}'.format(n_consensus)
                        cdr3 = Counter(cdr3 for record, cdr3 in cluster).most_common(1)[0][0]

                name = name.split(maxsplit=1)[0]
                if name.endswith(';'):
                    name = name[:-1]
                if group_by_cdr3:
                    print('>{};barcode={};cdr3={};size={};\n{}'.format(name, barcode,
                        cdr3, len(cluster), sequence))
                else:
                    print('>{};barcode={};size={};\n{}'.format(name, barcode,
                        len(cluster), sequence))

    logger.info('%d clusters (%d singletons)', n_clusters, n_singletons)
    logger.info('%d consensus sequences computed (from groups that had at least %d sequences)',
        n_consensus + n_ambiguous, MIN_CONSENSUS_SEQUENCES)
    logger.info('%d of those had no ambiguous bases', n_consensus)
    if args.groups_output:
        logger.info('Groups written to %r', args.groups_output)

    if args.json:
        sizes_counter = Counter(sizes)
        stats = {
            'unique_barcodes': len(barcodes),
            'barcode_singletons': barcode_singletons,
            'groups_written': n_clusters,
            'group_size_1': sizes_counter[1],
            'group_size_2': sizes_counter[2],
            'group_size_3plus': sum(v for k, v in sizes_counter.items() if k >= 3),
        }
        with open(args.json, 'w') as f:
            json.dump(stats, f, indent=2)
            print(file=f)

    if args.plot_sizes:
        plot_sizes(sizes, args.plot_sizes)
