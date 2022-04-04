"""
Discover D and J genes

The most frequent D/J sequences are considered candidates.

For J genes, candidate sequences are merged if they overlap each other.

The result table is written to standard output. Use --fasta to also
generate FASTA output.
"""
import logging
import pandas as pd
from collections import defaultdict
from typing import List
import dnaio

from tinyalign import edit_distance
from ..utils import Merger, merge_overlapping, unique_name, is_same_gene, slice_arg, UniqueNamer
from ..table import read_table, read_table_chunks

logger = logging.getLogger(__name__)

MINIMUM_CANDIDATE_LENGTH = 5


def add_arguments(parser):
    arg = parser.add_argument
    arg('--database', metavar='FASTA',
        help='FASTA file with reference gene sequences')
    arg('--merge', default=None, action='store_true', help='Merge overlapping genes. '
        'Default: Enabled for D, disabled for J and V.')
    arg('--no-merge', dest='merge', action='store_false', help='Do not merge overlapping genes')
    arg('--gene', default='J', choices=('V', 'D', 'J'), type=str.upper,
        help='Which gene category to discover. Default: %(default)s')
    arg('--j-coverage', type=float, default=None, metavar='PERCENT',
        help='Require that the sequence covers at least PERCENT of the J gene. '
        'Default: 90 when --gene=J; 0 otherwise')
    arg('--allele-ratio', type=float, metavar='RATIO', default=0.2,
        help='Required allele ratio. Works only for genes named "NAME*ALLELE". Default: %(default)s')
    arg('--cross-mapping-ratio', type=float, metavar='RATIO', default=0.1,
        help='Ratio for detection of cross-mapping artifacts. Default: %(default)s')
    arg('--min-count', metavar='N', type=int, default=None,
        help='Omit candidates with fewer than N exact occurrences in the input table. '
            'Default: 1 for J; 10 for D; 100 for V')
    arg('--no-perfect-matches', dest='perfect_matches', default=True, action='store_false',
        help='Do not filter out sequences for which the V assignment (or J for --gene=V) '
            'has at least one error')

    # --gene=D options
    arg('--d-core-length', metavar='L', type=int, default=6,
        help='Use only D core regions that have at least length L (only '
            'applies when --gene=D). Default: %(default)s')
    arg('--d-core', type=slice_arg, default=slice(2, -2),
        help='D core region location (only applies when --gene=D). '
            'Default: %(default)s')

    arg('--fasta', help='Write discovered sequences to FASTA file')
    arg('table', help='Table with parsed and filtered IgBLAST results')


class Candidate:
    __slots__ = ('name', 'sequence', 'exact_occ', 'max_count', 'other_genes', 'db_name',
        'db_distance', 'cdr3s', 'missing')

    def __init__(self, name, sequence, exact_occ=0, max_count=0, cdr3s=None, other_genes=None,
            db_name=None, db_distance=None):
        self.name = name
        self.sequence = sequence
        self.exact_occ = exact_occ
        self.max_count = max_count  # an upper bound for the exact_occ
        self.cdr3s = cdr3s if cdr3s is not None else set()
        self.other_genes = other_genes if other_genes is not None else set()
        self.db_name = db_name
        self.db_distance = db_distance
        self.missing = ''

    @property
    def unique_CDR3(self):
        return len(self.cdr3s)

    def __repr__(self):
        return 'Candidate({sequence!r}, exact_occ={exact_occ}, max_count={max_count}, ...)'.format(
            sequence=self.sequence,
            exact_occ=self.exact_occ,
            max_count=self.max_count,
        )


class OverlappingSequenceMerger(Merger):
    """
    Merge sequences that overlap
    """
    def merged(self, s, t):
        """
        Merge two sequences if they overlap. If they should not be merged,
        None is returned.
        """
        m = merge_overlapping(s.sequence, t.sequence)
        if m is not None:
            return Candidate(s.name, m, max_count=t.max_count + s.max_count)

        return None


class AlleleRatioMerger(Merger):
    """
    Discard sequences with too low allele ratio
    """
    def __init__(self, allele_ratio, cross_mapping_ratio):
        super().__init__()
        self._allele_ratio = allele_ratio
        self._cross_mapping_ratio = cross_mapping_ratio

    def merged(self, s, t):
        """
        Merge two sequences if they overlap. If they should not be merged,
        None is returned.
        """
        # TODO copy-and-pasted from germlinefilter
        #
        # Check allele ratio. Somewhat similar to cross-mapping, but
        # this uses sequence names to decide whether two genes can be
        # alleles of each other and the ratio is between the CDR3s_exact
        # values
        if self._allele_ratio and is_same_gene(s.name, t.name):
            for u, v in [(s, t), (t, s)]:
                if v.unique_CDR3 == 0:
                    continue
                ratio = u.unique_CDR3 / v.unique_CDR3
                if ratio < self._allele_ratio:
                    # logger.info('Allele ratio %.4f too low for %r compared to %r',
                    #     ratio, u.name, v.name)
                    return v

        if self._cross_mapping_ratio:
            # When checking for cross mapping, ignore overhanging bases in the 5' end.
            # Example:
            # ---ACTACGACTA...
            # XXX|||||X||||
            # ATTACTACTACTA...
            if len(t.sequence) < len(s.sequence):
                t, s = s, t  # s is now the shorter sequence
            t_seq = t.sequence[len(t.sequence) - len(s.sequence):]
            s_seq = s.sequence
            dist = edit_distance(s_seq, t_seq, 1)
            if dist > 1:
                return None
            total_occ = (s.exact_occ + t.exact_occ)
            if total_occ == 0:
                return None
            for u, v in [(s, t), (t, s)]:
                ratio = u.exact_occ / total_occ
                if ratio < self._cross_mapping_ratio:
                    # u is probably a cross-mapping artifact of the higher-expressed v
                    logger.info('%r is a cross-mapping artifact of %r (ratio %.4f)',
                        u.name, v.name, ratio)
                    return v

        return None


def filter_by_allele_ratio(table, allele_ratio):
    arm = AlleleRatioMerger(allele_ratio, cross_mapping_ratio=None)
    renamed_counts = table.reset_index().rename(columns={'gene': 'name'})
    arm.extend(renamed_counts.itertuples(index=False))
    counts = pd.DataFrame(list(arm), columns=renamed_counts.columns) \
        .rename(columns={'name': 'gene'}) \
        .set_index('gene')
    return counts


def count_occurrences(candidates, table_path, search_columns, other_gene, other_errors, merge,
        perfect_matches):
    """
    Count how often each candidate sequence occurs in the input table.
    The columns named in search_columns are concatenated and searched.
    This circumvents inaccurate IgBLAST alignment boundaries.
    Rows where the column named by other_errors is not zero are ignored.

    The input table is read in chunks to reduce memory usage.

    candidates -- list of candidates
    table_path -- path to input table
    search_columns -- which columns to search. The contained strings are
        concatenated and then searched.
    merge -- If True, stop searching for other candidates in a single row
        after one candidate has been found.

    The following attributes of the candidates are updated:

    - exact_occ
    - other_genes
    - cdr3s

    Return the updated list of candidates.
    """
    candidates_map = {c.sequence: c for c in candidates}

    search_order = [c.sequence for c in candidates]
    cols = [other_gene, 'V_errors', 'J_errors', 'cdr3'] + search_columns

    search_cache = defaultdict(list)  # map haystack sequence to list of candidates that occur in it
    for chunk in read_table_chunks(table_path, usecols=cols, chunksize=10000):
        if perfect_matches:
            chunk = chunk[chunk[other_errors] == 0].copy()
        # concatenate search columns
        if len(chunk) == 0:  # TODO that this is needed is possibly a pandas bug
            continue
        chunk['haystack'] = chunk.loc[:, search_columns].astype(str).sum(axis=1)
        chunk['haystack'] = chunk['haystack'].str.replace('(', '', regex=False).replace(')', '', regex=False)

        for row in chunk.itertuples():
            if row.haystack not in search_cache:
                for needle in search_order:
                    if needle in row.haystack:
                        search_cache[row.haystack].append(candidates_map[needle])
            for candidate in search_cache[row.haystack]:
                candidate.exact_occ += 1  # TODO += row.count?
                candidate.other_genes.add(getattr(row, other_gene))
                candidate.cdr3s.add(row.cdr3)
                if merge:
                    # When overlapping candidates have been merged,
                    # there will be no other pattern that is a
                    # substring of the current search pattern.
                    break
    return candidates_map.values()


def discard_substring_occurrences(candidates):
    """
    Filter a candidate list by discarding candidates whose sequences are
    substrings of another candidate’s sequence
    """
    # Shorter sequences first
    candidates = sorted(candidates, key=lambda c: len(c.sequence))
    for i, short in enumerate(candidates):
        for long in candidates[i+1:]:
            if short.sequence in long.sequence:
                break
        else:
            # no substring occurrence - keep this candidate
            yield short


def sequence_candidates(table, column, minimum_length, core=slice(None, None), min_occ=3):
    """
    Generate candidates by clustering all sequences in a column
    (V_nt, D_nt or J_nt). At least min_occ occurrences are required.

    core -- a slice object. If given, the strings in the column are
           sliced before being clustered.
    """
    for sequence, occ in table[column].str[core].value_counts().items():
        if len(sequence) >= minimum_length and occ >= min_occ:
            yield Candidate(None, sequence, max_count=occ)


def count_unique_cdr3(table):
    return len(set(s for s in table.cdr3 if s))


def count_unique_gene(table, gene_type):
    return len(set(s for s in table[gene_type.lower() + '_call'] if s))


def compute_expressions(table, gene_type):
    """Compute expression counts of known genes"""
    assert gene_type in {'V', 'D', 'J'}
    columns = ('gene', 'count', 'unique_CDR3')
    for gt in 'V', 'D', 'J':
        if gene_type != gt:
            columns += ('unique_' + gt,)
    gene_column = gene_type.lower() + '_call'
    rows = []
    for gene, group in table.groupby(gene_column):
        if gene == '':
            continue
        unique_cdr3 = count_unique_cdr3(group)
        row = dict(gene=gene, count=len(group), unique_CDR3=unique_cdr3)
        for gt in 'V', 'D', 'J':
            if gene_type != gt:
                row['unique_' + gt] = count_unique_gene(group, gene_type=gt)
        rows.append(row)
    counts = pd.DataFrame(rows, columns=columns).set_index('gene')
    return counts


def make_whitelist(table, database, gene_type: str, allele_ratio=None) -> List[str]:
    """
    Return a list of sequences that represent expressed alleles
    """
    assert gene_type in {'V', 'D', 'J'}
    # Compute expression counts in the same way the 'igdiscover count' command would do it
    counts = compute_expressions(table, gene_type)
    if allele_ratio:
        counts = filter_by_allele_ratio(counts, allele_ratio)
    names = list(counts.index)

    # Construct whitelist from all expressed alleles
    database = {r.name: r.sequence for r in database}
    sequences = [database[name] for name in names]

    return sequences


def print_table(candidates, other_gene, missing):
    columns = ['name', 'exact_occ', other_gene + 's', 'CDR3s', 'database', 'database_diff', 'sequence']
    if missing:
        columns.append('missing')
    print(*columns, sep='\t')
    for candidate in candidates:
        columns = [
            candidate.name,
            candidate.exact_occ,
            len(set(candidate.other_genes)),
            candidate.unique_CDR3,
            candidate.db_name if candidate.db_name is not None else '',
            candidate.db_distance if candidate.db_distance is not None else -1,
            candidate.sequence
        ]
        if missing:
            columns.append(candidate.missing)
        print(*columns, sep='\t')


def main(args):
    if args.database:
        with dnaio.open(args.database) as fr:
            database = list(fr)
        logger.info('Read %d sequences from %r', len(database), args.database)
    else:
        database = None
    column = {'V': 'V_nt', 'J': 'J_nt', 'D': 'D_region'}[args.gene]
    other = 'V' if args.gene in ('D', 'J') else 'J'
    other_gene = other.lower() + '_call'
    other_errors = other + '_errors'
    table = read_table(args.table,
        usecols=['count', 'v_call', 'd_call', 'j_call', 'V_errors', 'J_errors', 'J_covered', column, 'cdr3'])
    logger.info('Table with %s rows read', len(table))

    if args.j_coverage is None and args.gene == 'J':
        args.j_coverage = 90
    if args.j_coverage:
        table = table[table['J_covered'] >= args.j_coverage]
        logger.info('Keeping %s rows that have J_covered >= %s', len(table), args.j_coverage)
    if args.perfect_matches:
        table = table[table[other_errors] == 0]
        logger.info('Keeping %s rows that have no %s mismatches', len(table), other)

    if args.merge is None:
        args.merge = args.gene == 'D'
    if args.min_count is None:
        args.min_count = {'J': 1, 'D': 10, 'V': 100}[args.gene]  # TODO J is fine, but are D and V?

    if args.gene == 'D':
        candidates = sequence_candidates(
            table, column, minimum_length=args.d_core_length, core=args.d_core)
    elif args.gene == 'J':
        candidates = sequence_candidates(
            table, column, minimum_length=MINIMUM_CANDIDATE_LENGTH)
    else:
        candidates = sequence_candidates(
            table, column, minimum_length=MINIMUM_CANDIDATE_LENGTH)

    candidates = list(candidates)
    logger.info('Collected %s unique %s sequences', len(candidates), args.gene)

    # Add whitelisted sequences
    if database:
        whitelist = make_whitelist(table, database, args.gene, args.allele_ratio)
        missing_whitelisted = set(whitelist) - set(c.sequence for c in candidates)
        for sequence in missing_whitelisted:
            candidates.append(Candidate(None, sequence))
        logger.info('Added %d whitelisted sequence%s',
            len(missing_whitelisted), 's' if len(missing_whitelisted) != 1 else '')

    candidates = list(discard_substring_occurrences(candidates))
    logger.info('Removing candidate sequences that occur within others results in %s candidates',
        len(candidates))
    candidates = [candidate for candidate in candidates if 'N' not in candidate.sequence]
    logger.info('Removing candidates containing "N" results in %s candidates',
        len(candidates))

    if args.merge:
        logger.info('Merging overlapping sequences ...')
        # Merge candidate sequences that overlap. If one candidate is longer than
        # another, this is typically a sign that IgBLAST has not extended the
        # alignment long enough.
        merger = OverlappingSequenceMerger()
        for candidate in candidates:
            merger.add(candidate)
        logger.info('After merging overlapping %s sequences, %s remain', args.gene, len(merger))
        candidates = list(merger)
        logger.info('%d candidates', len(candidates))
    del table

    # Assign names etc.
    if database:
        for candidate in candidates:
            distances = [(edit_distance(db.sequence, candidate.sequence), db) for db in database]
            candidate.db_distance, closest = min(distances, key=lambda x: x[0])
            candidate.db_name = closest.name

            if candidate.db_distance == 0:
                candidate.name = closest.name
            else:
                # Exact db sequence not found, is there one that contains
                # this candidate as a substring?
                for db_record in database:
                    index = db_record.sequence.find(candidate.sequence)
                    if index == -1:
                        continue
                    if args.gene == 'D':
                        start = db_record.sequence.find(candidate.sequence)
                        prefix = db_record.sequence[:start]
                        suffix = db_record.sequence[start + len(candidate.sequence):]
                        candidate.missing = '{}...{}'.format(prefix, suffix)
                    else:
                        # Replace this record with the full-length version
                        candidate.sequence = db_record.sequence
                        candidate.db_distance = 0
                    candidate.name = db_record.name
                    break
                else:
                    candidate.name = unique_name(closest.name, candidate.sequence)
    else:
        for candidate in candidates:
            candidate.name = unique_name(args.gene, candidate.sequence)

    logger.info('Counting occurrences ...')
    if args.gene == 'D':
        search_columns = ['np1', 'D_region', 'np2']
    elif args.gene == 'J':
        search_columns = ['np2', 'J_nt']
    else:
        search_columns = ['sequence']
    candidates = count_occurrences(candidates, args.table, search_columns, other_gene, other_errors,
        args.merge, args.perfect_matches)

    # Filter by allele ratio
    if args.allele_ratio or args.cross_mapping_ratio:
        arm = AlleleRatioMerger(args.allele_ratio, args.cross_mapping_ratio)
        arm.extend(candidates)
        candidates = list(arm)
        logger.info('After filtering by allele ratio and/or cross-mapping ratio, %d candidates remain',
            len(candidates))

    candidates = sorted(candidates, key=lambda c: c.name)
    candidates = [c for c in candidates if c.exact_occ >= args.min_count or c.db_distance == 0]
    namer = UniqueNamer()
    for candidate in candidates:
        candidate.name = namer(candidate.name)
    print_table(candidates, other_gene, missing=args.gene == 'D')

    if args.fasta:
        with open(args.fasta, 'w') as f:
            for candidate in sorted(candidates, key=lambda r: r.name):
                print('>{}\n{}'.format(candidate.name, candidate.sequence), file=f)

    logger.info('Wrote %d genes', len(candidates))
