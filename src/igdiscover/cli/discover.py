"""
Discover candidate new V genes within a single antibody library

Existing V sequences are grouped by their V gene assignment, and within each
group, consensus sequences are computed.
"""
import csv
import sys
import hashlib
import logging
import random
import itertools
import multiprocessing
from contextlib import ExitStack
from collections import namedtuple, Counter

import dnaio
import numpy as np
import pandas as pd
from tinyalign import edit_distance

from ..table import read_table, vdj_nt_column
from ..cluster import cluster_sequences, single_linkage
from ..utils import (unique_name, downsampled, SerialPool, Merger, has_stop,
    available_cpu_count, UniqueNamer)
from ..align import iterative_consensus, describe_nt_change

logger = logging.getLogger(__name__)


MINGROUPSIZE = 5
MINEXPRESSED = 10
MAXIMUM_SUBSAMPLE_SIZE = 1600


def add_arguments(parser):
    arg = parser.add_argument
    arg('--threads', '-j', type=int, default=min(4, available_cpu_count()),
        help='Number of threads. Default: no. of available CPUs, but at most 4')
    arg('--seed', type=int, default=None,
        help='Seed value for random numbers for reproducible runs.')
    arg('--consensus-threshold', '-t', metavar='PERCENT', type=float, default=60,
        help='Threshold for consensus computation. Default: %(default)s%%')
    arg('--gene', '-g', action='append', default=[],
        help='Compute consensus for this gene. Can be given multiple times. '
            'Default: Compute for all genes.')
    arg('--limit', type=int, default=None, metavar='N',
        help='Skip remaining genes as soon as at least N candidates were generated. Default: No limit')
    arg('--left', '-l', type=float, metavar='ERROR-RATE',
        help='For consensus, include only sequences that have at least this '
            'error rate (in percent). Default: %(default)s', default=0)
    arg('--right', '-r', type=float, metavar='ERROR-RATE',
        help='For consensus, include only sequences that have at most this '
            'error rate (in percent). Default: %(default)s', default=100)
    arg('--window-width', '-w', type=float, metavar='PERCENT',
        help='Compute consensus for all PERCENT-wide windows. Set to 0 to '
            'disable. Default: %(default)s', default=2)
    arg('--no-cluster', dest='cluster', action='store_false', default=True,
        help='Do not run linkage cluster analysis.')
    arg('--cluster-exact', metavar='N', type=int, default=0,
        help='Treat N exact occurrences of a sequence as a cluster. '
            'Default: Do not cluster exact occurrences')
    arg('--max-n-bases', type=int, default=0, metavar='MAXN',
        help='Remove rows that have more than MAXN "N" nucleotides. If >0, an '
            'N_bases column is added. Default: %(default)s')
    arg('--subsample', metavar='N', type=int, default=1000,
        help='When clustering, use N randomly chosen sequences. Default: %(default)s')
    arg('--ignore-J', action='store_true', default=False,
        help='Include also rows without J assignment or J%%SHM>0.')
    arg('--exact-copies', metavar='N', type=int, default=1,
        help='When subsampling, first pick rows whose V gene sequences'
            'have at least N exact copies in the input. Default: %(default)s')
    arg('--d-evalue', metavar='EVALUE', type=float, default=1E-4,
        help='For Ds_exact, require D matches with an E-value of '
            'at most EVALUE. Default: %(default)s')
    arg('--d-coverage', '--D-coverage', metavar='COVERAGE', type=float, default=70,
        help='For Ds_exact, require D matches with a minimum D '
            'coverage of COVERAGE (in percent). Default: %(default)s)')
    arg('--clonotype-diff', metavar='DIFFERENCES', type=int, default=6,
        help='When clustering CDR3s to computer the no. of clonotypes, allow DIFFERENCES '
            'between (nucleotide-)sequences. Default: %(default)s')

    arg('--table-output', '-o', metavar='DIRECTORY',
        help='Output tables for all analyzed genes to DIRECTORY. '
            'Files will be named <GENE>.tab.')
    arg('--database', metavar='FASTA', default=None,
        help='FASTA file with V genes. If provided, differences between consensus '
            'and database will be computed.')
    arg('--read-names', metavar='FILE',
        help='Write names of reads with exact matches used in discovering each candidate '
            'to FILE')
    arg('table', help='Table with parsed IgBLAST results')


Groupinfo = namedtuple('Groupinfo',
    'count unique_D unique_J unique_CDR3 shared_CDR3_ratio clonotypes read_names unique_barcodes')

SiblingInfo = namedtuple('SiblingInfo', 'sequence requested name group')


def safe_divide(x, y):
    return x / y if y != 0 else 0


class SiblingMerger(Merger):
    """
    Merge very similar consensus sequences into single entries. This could be
    seen as a type of clustering using very specific criteria. Two sequences
    are merged if one is the prefix of the other, allowing differences where
    one of the sequences has an 'N' base.
    """
    def merged(self, s, t):
        chars = []
        for c1, c2 in itertools.zip_longest(s.sequence, t.sequence):
            if c1 is None:
                c = c2
            elif c2 is None:
                c = c1
            elif c1 == 'N':
                c = c2
            elif c2 == 'N':
                c = c1
            elif c1 != c2:
                return None
            else:
                assert c1 == c2
                c = c1
            chars.append(c)
        seq = ''.join(chars)
        requested = s.requested or t.requested
        name = s.name + ';' + t.name
        # take union of groups
        group = pd.concat([s.group, t.group]).groupby(level=0).last()
        return SiblingInfo(seq, requested, name, group)


class Discoverer:
    """
    Discover candidates for novel V genes.
    """
    def __init__(
        self,
        database,
        windows,
        left,
        right,
        cluster,
        cluster_exact,
        table_output,
        consensus_threshold,
        downsample,
        clonotype_differences,
        cluster_subsample_size,
        max_n_bases,
        exact_copies,
        d_coverage,
        d_evalue,
        seed,
        cdr3_counts
    ):
        """
        :param cdr3_counts:
        """
        self.database = database
        self.windows = windows
        self.left = left
        self.right = right
        self.cluster = cluster
        self.cluster_exact = cluster_exact
        self.table_output = table_output
        self.consensus_threshold = consensus_threshold
        self.downsample = downsample
        self.clonotype_differences = clonotype_differences
        self.cluster_subsample_size = cluster_subsample_size
        self.max_n_bases = max_n_bases
        self.exact_copies = exact_copies
        self.d_coverage = d_coverage
        self.d_evalue = d_evalue
        self.seed = seed
        self.cdr3_counts = cdr3_counts

    def _sibling_sequence(self, gene, group):
        """
        For a given group, compute a consensus sequence over the V gene sequences
        in that group.

        If the found sibling is slightly longer or shorter than the version in
        the database, adjust it so it corresponds to the database version exactly.
        """
        sequence = iterative_consensus(list(group.V_nt), program='muscle-medium',
            threshold=self.consensus_threshold/100,
            maximum_subsample_size=self.downsample)
        if gene in self.database:
            database_sequence = self.database[gene]
            if sequence.startswith(database_sequence) or database_sequence.startswith(sequence):
                return database_sequence
        return sequence

    @staticmethod
    def _guess_chain(group):
        """
        Return a guess for the chain type of a given group
        """
        return Counter(group.locus).most_common()[0][0]

    @staticmethod
    def _guess_cdr3_start(group):
        """
        Return a guess for the CDR3 start within sequences in the given group
        """
        return Counter(group.V_CDR3_start).most_common()[0][0]

    def count_unique_d(self, table):
        g = table[(table.D_errors == 0) &
            (table.D_covered >= self.d_coverage) &
            (table.d_support <= self.d_evalue)]
        return len(set(s for s in g.d_call if s))

    @staticmethod
    def count_unique_barcodes(group):
        return len(set(s for s in group.barcode if s))

    def count_clonotypes(self, table):
        """
        Cluster sequences by edit distance and return the number of clusters.

        The sequences are first group by their J assignment and cluster
        numbers are computed within these groups separately, then summed up.
        """
        distance = self.clonotype_differences

        def linked(s, t):
            return edit_distance(s, t, distance) <= distance

        total = 0
        for j_gene, group in table.groupby('j_call'):
            sequences = list(set(s for s in group.cdr3 if s))
            components = single_linkage(sequences, linked)
            total += len(components)
        return total

    def _cluster_siblings(self, gene, table):
        """Find candidates by clustering sequences assigned to one gene"""

        if self.exact_copies > 1:
            # Preferentially pick those sequences (for subsampling) that have
            # multiple exact copies, then fill up with the others
            exact_group = table[table.copies >= self.exact_copies]
            indices = downsampled(list(exact_group.index), self.cluster_subsample_size)
            if len(indices) < self.cluster_subsample_size:
                not_exact_group = table[table.copies < self.exact_copies]
                indices.extend(downsampled(list(not_exact_group.index),
                    self.cluster_subsample_size - len(indices)))
        else:
            indices = downsampled(list(table.index), self.cluster_subsample_size)
        # Ignore CDR3 part of the V sequence for clustering
        sequences_no_cdr3 = list(table.V_no_CDR3.loc[indices])
        df, linkage, clusters = cluster_sequences(sequences_no_cdr3, MINGROUPSIZE)
        n_clusters = len(set(clusters))
        logger.info('%6d %s assignments generated %d cluster%s',
            len(table), gene, n_clusters, 's' if n_clusters != 1 else '')
        cluster_indices = [[] for _ in range(max(clusters) + 1)]
        for i, cluster_id in enumerate(clusters):
            cluster_indices[cluster_id].append(indices[i])

        cl = 0
        for ind in cluster_indices:
            group = table.loc[ind]
            if len(group) < MINGROUPSIZE:
                continue
            sibling = self._sibling_sequence(gene, group)
            cl += 1
            yield SiblingInfo(
                sequence=sibling, requested=False, name='cl{}'.format(cl), group=group)

    def _window_siblings(self, gene, table):
        """
        Find candidates by clustering sequences that have a similar number of differences
        to the reference sequence
        """
        for left, right in self.windows:
            left, right = float(left), float(right)
            group = table[(left <= table.V_SHM) & (table.V_SHM < right)]
            if len(group) < MINGROUPSIZE:
                continue
            sibling = self._sibling_sequence(gene, group)
            if left == int(left):
                left = int(left)
            if right == int(right):
                right = int(right)
            requested = (left, right) == (self.left, self.right)
            if (left, right) == (0, 100):
                name = 'all'
            else:
                name = '{}-{}'.format(left, right)
            yield SiblingInfo(sequence=sibling, requested=requested, name=name, group=group)

    def _cluster_exact_candidates(self, gene, table):
        index = 1
        for sequence, group in table.groupby('V_nt'):
            if len(group) >= self.cluster_exact:
                name = 'ex{}'.format(index)
                index += 1
                yield SiblingInfo(sequence, False, name, group)

    def _collect_siblings(self, gene, group):
        """
        gene -- gene name
        group -- pandas.DataFrame of sequences assigned to that gene

        Yield SiblingInfo objects
        """
        group = group.copy()
        # the original reference sequence for all the IgBLAST assignments in this group
        database_sequence = self.database.get(gene, None)
        database_sequence_found = False

        candidate_iterators = [self._window_siblings(gene, group)]
        if self.cluster:
            candidate_iterators.append(self._cluster_siblings(gene, group))
        if self.cluster_exact:
            candidate_iterators.append(self._cluster_exact_candidates(gene, group))
        for sibling in itertools.chain(*candidate_iterators):
            if sibling.sequence == database_sequence:
                database_sequence_found = True
            yield sibling

        if database_sequence:
            # If this is a database sequence and there are some exact occurrences of it,
            # add it to the list of candidates even if it has not been found as a
            # cluster.
            group_in_window = group[group.V_errors == 0]
            if len(group_in_window) >= MINEXPRESSED:
                if not database_sequence_found:
                    logger.info('Database sequence %r seems to be expressed, but is missing from '
                        'candidates. Re-adding it.', gene)
                yield SiblingInfo(database_sequence, False, 'db', group_in_window)

    def set_random_seed(self, name):
        """Set random seed depending on gene name and seed given to constructor"""
        h = hashlib.md5(name.encode()).digest()[:4]
        n = int.from_bytes(h, byteorder='big')
        random.seed(n + self.seed)

    def __call__(self, args):
        """
        Discover new V genes. args is a tuple (gene, group)
        gene -- name of the gene
        assignments -- a pandas DataFrame with the assignments to the gene
        """
        gene, assignments = args
        self.set_random_seed(gene)
        siblings = SiblingMerger()
        for sibling in self._collect_siblings(gene, assignments):
            siblings.add(sibling)

        candidates = []
        for sibling_info in siblings:
            sibling = sibling_info.sequence
            n_bases = sibling.count('N')
            if n_bases > self.max_n_bases:
                logger.debug('Sibling %s has too many N bases', sibling_info.name)
                continue

            # Sequence without the CDR3-covering part
            sibling_no_cdr3 = sibling[:self._guess_cdr3_start(assignments)]
            group_exact_v = assignments[assignments.V_no_CDR3 == sibling_no_cdr3]

            group_full_exact_v = assignments[assignments['VDJ_nt'].str.startswith(sibling, na=False)]
            groups = (
                ('window', sibling_info.group),
                ('exact', group_exact_v),
                ('full_exact', group_full_exact_v),
            )
            del sibling_no_cdr3

            # self.cdr3_counts are CDR3 counts of all CDR3s in the entire table.
            # We restrict this here to the counts for CDR3s belonging to
            # clusters other than the current one.
            other_cdr3_counts = self.cdr3_counts - Counter(s for s in sibling_info.group.cdr3 if s)
            info = dict()
            for key, group in groups:
                cdr3_counts = Counter(s for s in group.cdr3 if s)
                unique_cdr3 = len(cdr3_counts)
                shared_cdr3_ratio = safe_divide(len(other_cdr3_counts & cdr3_counts), unique_cdr3)
                unique_j = len(set(s for s in group.j_call if s))
                clonotypes = self.count_clonotypes(group)
                unique_d = self.count_unique_d(group)
                unique_barcodes = self.count_unique_barcodes(group)
                count = len(group.index)
                read_names = list(group.sequence_id)
                info[key] = Groupinfo(count=count, unique_D=unique_d, unique_J=unique_j,
                    unique_CDR3=unique_cdr3, shared_CDR3_ratio=shared_cdr3_ratio,
                    clonotypes=clonotypes, read_names=read_names,
                    unique_barcodes=unique_barcodes)
            if gene in self.database:
                database_diff = edit_distance(sibling, self.database[gene])
                database_changes = describe_nt_change(self.database[gene], sibling)
            else:
                database_diff = None
                database_changes = None

            # Build the Candidate
            sequence_id = gene if database_diff == 0 else unique_name(gene, sibling)
            chain = self._guess_chain(sibling_info.group)
            cdr3_start = self._guess_cdr3_start(sibling_info.group)
            ratio = safe_divide(info['exact'].count, info['exact'].unique_CDR3)

            # Apply some very light filtering on non-database sequences
            if database_diff > 0 and info['exact'].count < 2:
                continue

            candidate = Candidate(
                name=sequence_id,
                source=gene,
                chain=chain,
                cluster=sibling_info.name,
                cluster_size=info['window'].count,
                Js=info['window'].unique_J,
                CDR3s=info['window'].unique_CDR3,
                exact=info['exact'].count,
                full_exact=info['full_exact'].count,
                barcodes_exact=info['exact'].unique_barcodes,
                Ds_exact=info['exact'].unique_D,
                Js_exact=info['exact'].unique_J,
                CDR3s_exact=info['exact'].unique_CDR3,
                clonotypes=info['exact'].clonotypes,
                CDR3_exact_ratio=ratio,
                CDR3_shared_ratio=info['exact'].shared_CDR3_ratio,
                N_bases=n_bases,
                database_diff=database_diff,
                database_changes=database_changes,
                has_stop=has_stop(sibling),
                CDR3_start=cdr3_start,
                consensus=sibling,
                read_names=info['exact'].read_names,
            )
            candidates.append(candidate)
        return candidates


class Candidate(namedtuple('_Candidate', [
    'name',
    'source',
    'chain',
    'cluster',
    'cluster_size',
    'Js',
    'CDR3s',
    'exact',
    'full_exact',
    'barcodes_exact',
    'Ds_exact',
    'Js_exact',
    'CDR3s_exact',
    'clonotypes',
    'CDR3_exact_ratio',
    'CDR3_shared_ratio',
    'N_bases',
    'database_diff',
    'database_changes',
    'has_stop',
    'CDR3_start',
    'consensus',
    'read_names',
])):
    def formatted_dict(self):
        d = self._asdict()
        d['has_stop'] = int(d['has_stop'])
        for name in 'CDR3_exact_ratio', 'CDR3_shared_ratio':
            d[name] = '{:.2f}'.format(d[name])
        del d['read_names']
        return d


def count_prefixes(sequences):
    """
    Count how often each sequence occurs in the given list of
    sequences. If one sequence is the prefix of another one,
    they are considered to be 'identical'.

    Return a dictionary that maps sequence to count.

    >>> r = count_prefixes(['A', 'BB', 'CD', 'CDE', 'CDEF'])
    >>> r == {'A': 1, 'BB': 1, 'CD': 3, 'CDE': 3, 'CDEF': 3}
    True
    """
    sequences = sorted(sequences)
    sequences.append('GUARD')
    prev = 'X'
    start = 0
    count = dict()
    for i, s in enumerate(sequences):
        if not s.startswith(prev):
            # end of a run
            for j in range(start, i):
                count[sequences[j]] = i - start
            start = i
        prev = s
    return count


def main(args):
    if args.database:
        with dnaio.open(args.database) as sr:
            database = {record.name: record.sequence.upper() for record in sr}
    else:
        database = dict()

    if args.seed:
        seed = args.seed
    else:
        seed = random.randrange(10**6)
        logger.info('Use --seed=%d to reproduce this run', seed)

    table = read_table(args.table, usecols=('sequence_id', 'locus', 'v_call', 'd_call', 'j_call',
        'V_nt', 'cdr3', 'barcode', 'V_CDR3_start', 'V_SHM', 'J_SHM', 'D_covered', 'd_support',
        'V_errors', 'D_errors', 'J_errors', 'VDJ_nt'))
    table['V_no_CDR3'] = [s[:start] if start != 0 else s for s, start in
        zip(table.V_nt, table.V_CDR3_start)]

    logger.info('%s rows read', len(table))
    if not args.ignore_J:
        # Discard rows with any mutation within J at all
        table = table[table.J_SHM == 0][:]
        logger.info('%s rows remain after discarding J%%SHM > 0', len(table))

    if args.exact_copies > 1:
        multiplicities = count_prefixes(table.V_no_CDR3)
        table['copies'] = table.V_no_CDR3.map(multiplicities)
        logger.info('%s rows contain V sequences with at least %s copies',
            sum(table.copies >= args.exact_copies), args.exact_copies)

    columns = list(Candidate._fields)
    if not args.max_n_bases:
        columns.remove('N_bases')
    columns.remove('read_names')
    writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    genes = set(args.gene)
    if args.window_width:
        windows = [(start, start + args.window_width)
            for start in np.arange(0, 20, args.window_width)]
        logger.info('Using an error rate window of %.1f%% to %.1f%%', args.left, args.right)
        windows.append((args.left, args.right))
    else:
        windows = []

    groups = []
    for gene, group in table.groupby('v_call'):
        if genes and gene not in genes:
            continue
        if len(group) < MINGROUPSIZE:
            continue
        groups.append((gene, group))

    cdr3_counts = Counter(s for s in table.cdr3 if s)
    logger.info('%s unique CDR3s detected overall', len(cdr3_counts))
    discoverer = Discoverer(
        database,
        windows,
        args.left,
        args.right,
        args.cluster,
        args.cluster_exact,
        args.table_output,
        args.consensus_threshold,
        MAXIMUM_SUBSAMPLE_SIZE,
        clonotype_differences=args.clonotype_diff,
        cluster_subsample_size=args.subsample,
        max_n_bases=args.max_n_bases,
        exact_copies=args.exact_copies,
        d_coverage=args.d_coverage,
        d_evalue=args.d_evalue,
        seed=seed,
        cdr3_counts=cdr3_counts
    )

    Pool = SerialPool if args.threads == 1 else multiprocessing.Pool
    n_candidates = 0
    read_names_file = None
    with ExitStack() as stack:
        pool = stack.enter_context(Pool(args.threads))
        if args.read_names:
            read_names_file = stack.enter_context(open(args.read_names, 'w'))
        namer = UniqueNamer()
        for candidates in pool.imap(discoverer, groups, chunksize=1):
            for candidate in candidates:
                candidate = candidate._replace(name=namer(candidate.name))
                writer.writerow(candidate.formatted_dict())
                if read_names_file:
                    print(candidate.name, *candidate.read_names, sep='\t', file=read_names_file)
            n_candidates += len(candidates)
            if args.limit is not None and n_candidates >= args.limit:
                break
            sys.stdout.flush()
    logger.info('%s candidate sequences for %s gene(s) generated', n_candidates, len(groups))
