"""
Query a table of assigned sequences by clonotype

Two sequences have the same clonotype if
- their V and J assignments are the same
- the length of their CDR3 is identical
- the difference between their CDR3s (in terms of mismatches)
  is not higher than a given threshold (by default 1)

Clonotypes for the query sequences are determined and sequences
in the input table that have this clonotype are reported.

The table is written to standard output.
"""
import logging
from collections import defaultdict
from contextlib import ExitStack

import pandas as pd
from xopen import xopen

from ..table import read_table
from ..utils import slice_arg
from .clonotypes import is_similar_with_junction, CLONOTYPE_COLUMNS, augment_group

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--minimum-count', '-c', metavar='N', default=1, type=int,
        help='Discard all rows with count less than N. Default: %(default)s')
    arg('--cdr3-core', default=None,
        type=slice_arg, metavar='START:END',
        help='START:END defines the non-junction region of CDR3 '
            'sequences. Use negative numbers for END to count '
            'from the end. Regions before and after are considered to '
            'be junction sequence, and for two CDR3s to be considered '
            'similar, at least one of the junctions must be identical. '
            'Default: no junction region.')
    arg('--mismatches', default=1, type=float,
        help='No. of allowed mismatches between CDR3 sequences. '
            'Can also be a fraction between 0 and 1 (such as 0.15), '
            'interpreted relative to the length of the CDR3 (minus the front non-core). '
            'Default: %(default)s')
    arg('--aa', default=False, action='store_true',
        help='Count CDR3 mismatches on amino-acid level. Default: Compare nucleotides.')
    arg('--summary', metavar='FILE',
        help='Write summary table to FILE')
    arg('reftable', help='Reference table with parsed and filtered '
        'IgBLAST results (filtered.tab)')
    arg('querytable', help='Query table with IgBLAST results (assigned.tab or filtered.tab)')


def collect(querytable, reftable, mismatches, cdr3_core_slice, cdr3_column):
    """
    Find all queries from the querytable in the reftable.

    Yield tuples (query_rows, similar_rows) where the query_rows is a list
    with all the rows that have the same result. similar_rows is a DataFrame
    whose rows are the ones matching the query.
    """

    # The vjlentype is a "clonotype without CDR3 sequence" (only V, J, CDR3 length)
    # Determine set of vjlentypes to query
    query_vjlentypes = defaultdict(list)
    for row in querytable.itertuples():
        vjlentype = (row.v_call, row.j_call, len(row.cdr3))
        query_vjlentypes[vjlentype].append(row)

    groupby = ['v_call', 'j_call', 'CDR3_length']
    for vjlentype, vjlen_group in reftable.groupby(groupby):
        # (v_gene, j_gene, cdr3_length) = vjlentype
        if vjlentype not in query_vjlentypes:
            continue

        # Collect results for this vjlentype. The results dict
        # maps row indices (into the vjlen_group) to each query_row,
        # allowing us to group identical results together.
        results = defaultdict(list)
        for query_row in query_vjlentypes.pop(vjlentype):
            cdr3 = getattr(query_row, cdr3_column)
            # Save indices of the rows that are similar to this query
            indices = tuple(index for index, r in enumerate(vjlen_group.itertuples())
                if is_similar_with_junction(cdr3, getattr(r, cdr3_column), mismatches, cdr3_core_slice))
            results[indices].append(query_row)

        # Yield results, grouping queries that lead to the same result
        for indices, query_rows in results.items():
            if not indices:
                for query_row in query_rows:
                    yield ([query_row], reftable.head(0))
                continue

            similar_group = vjlen_group.iloc[list(indices), :].copy()
            yield (query_rows, similar_group)

    # Yield result tuples for all the queries that have not been found
    for queries in query_vjlentypes.values():
        for query_row in queries:
            yield ([query_row], reftable.head(0))


def main(args):
    usecols = CLONOTYPE_COLUMNS
    # TODO backwards compatibility
    if ('FR1_aa_mut' not in pd.read_table(args.querytable, nrows=0).columns or
            'FR1_aa_mut' not in pd.read_table(args.reftable, nrows=0).columns):
        usecols = [col for col in usecols if not col.endswith('_aa_mut')]
    querytable = read_table(args.querytable, usecols=usecols)
    querytable = querytable[usecols]  # reorder columns
    # Filter empty rows (happens sometimes)
    querytable = querytable[querytable.v_call != '']
    logger.info('Read query table with %s rows', len(querytable))
    reftable = read_table(args.reftable, usecols=usecols)
    reftable = reftable[usecols]
    logger.info('Read reference table with %s rows', len(reftable))
    if args.minimum_count > 1:
        reftable = reftable[reftable['count'] >= args.minimum_count]
        logger.info('After filtering out rows with count < %s, %s rows remain', args.minimum_count,
            len(reftable))
    for tab in querytable, reftable:
        tab.insert(5, 'CDR3_length', tab['cdr3'].apply(len))

    if len(querytable) > len(reftable):
        logger.warning('The reference table is smaller than the '
            'query table! Did you swap query and reference?')

    cdr3_column = 'cdr3_aa' if args.aa else 'cdr3'
    summary_columns = ['FR1_SHM', 'CDR1_SHM', 'FR2_SHM', 'CDR2_SHM', 'FR3_SHM', 'V_SHM', 'J_SHM',
        'V_aa_mut', 'J_aa_mut']
    summary_columns.extend(col for col in usecols if col.endswith('_aa_mut'))
    with ExitStack() as stack:
        if args.summary:
            summary_file = stack.enter_context(xopen(args.summary, 'w'))
            print('name', 'size', *('avg_' + s for s in summary_columns), sep='\t', file=summary_file)
        else:
            summary_file = None

        # Header
        print(*reftable.columns, sep='\t')

        for query_rows, result_table in collect(querytable, reftable, args.mismatches,
                args.cdr3_core, cdr3_column):
            assert len(query_rows) >= 1
            if summary_file:
                for query_row in query_rows:
                    print(query_row.sequence_id, len(result_table), sep='\t', end='', file=summary_file)
                    for col in summary_columns:
                        mean = result_table[col].mean() if len(result_table) > 0 else 0
                        print('\t{:.2f}'.format(mean), end='', file=summary_file)
                    print(file=summary_file)

            for query_row in query_rows:
                print('# Query: {}'.format(query_row.sequence_id), '', *(query_row[3:]), sep='\t')
            if len(result_table) > 0:
                print(result_table.to_csv(sep='\t', header=False, index=False))
            else:
                print()
