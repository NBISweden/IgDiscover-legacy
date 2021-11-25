"""
Filter table with parsed IgBLAST results

Discard the following rows in the table:
- no J assigned
- stop codon found
- V gene coverage less than 90%
- J gene coverage less than 60%
- V gene E-value greater than 1E-3

The filtered table is printed to standard output.
"""
import logging
import json
import pandas as pd

from ..table import read_table_chunks

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--v-coverage', type=float, default=90, metavar='PERCENT',
        help='Require that the sequence covers at least PERCENT of the V gene. '
        'Default: %(default)s')
    arg('--v-evalue', type=float, default=1E-3, metavar='EVALUE',
        help='Require that the E-value for the V gene match is at most EVALUE. '
        'Default: %(default)s')
    arg('--j-coverage', type=float, default=60, metavar='PERCENT',
        help='Require that the sequence covers at least PERCENT of the J gene. '
        'Default: %(default)s')
    arg('--json', metavar='FILE', help='Write statistics to FILE')
    arg('table', help='Table with filtered IgBLAST results.')


class FilteringStatistics:
    __slots__ = ('total', 'has_vj_assignment', 'has_no_stop', 'good_v_evalue', 'good_v_coverage',
    'good_j_coverage', 'has_cdr3')

    def __init__(self):
        self.total = 0
        self.has_vj_assignment = 0
        self.has_no_stop = 0
        self.good_v_evalue = 0
        self.good_v_coverage = 0
        self.good_j_coverage = 0
        self.has_cdr3 = 0

    def __iadd__(self, other):
        for att in self.__slots__:
            v = getattr(self, att)
            setattr(self, att, v + getattr(other, att))
        return self

    def asdict(self):
        d = dict()
        for att in self.__slots__:
            d[att] = getattr(self, att)
        return d


def filtered_table(table,
        v_gene_coverage,  # at least
        j_gene_coverage,  # at least
        v_gene_evalue,  # at most
    ):
    """
    Discard the following rows in the table:
    - no J assigned
    - stop codon found
    - V gene coverage less than v_gene_coverage
    - J gene coverage less than j_gene_coverage
    - V gene E-value greater than v_gene_evalue

    Return the filtered table.
    """
    stats = FilteringStatistics()
    stats.total = len(table)
    # Both V and J must be assigned
    # (Note v_call and j_call columns use empty strings instead of NA)
    filtered = table[(table['v_call'] != '') & (table['j_call'] != '')][:]
    stats.has_vj_assignment = len(filtered)
    filtered['v_call'] = pd.Categorical(filtered['v_call'])

    # Filter out sequences that have a stop codon
    filtered = filtered[filtered.stop_codon == 'F']
    stats.has_no_stop = len(filtered)

    # Filter out sequences with a too low V gene hit E-value
    filtered = filtered[filtered.v_support <= v_gene_evalue]
    stats.good_v_evalue = len(filtered)

    # Filter out sequences with too low V gene coverage
    filtered = filtered[filtered.V_covered >= v_gene_coverage]
    stats.good_v_coverage = len(filtered)

    # Filter out sequences with too low J gene coverage
    filtered = filtered[filtered.J_covered >= j_gene_coverage]
    stats.good_j_coverage = len(filtered)

    stats.has_cdr3 = sum(filtered['cdr3'] != '')
    return filtered, stats


def main(args):
    first = True
    written = 0
    stats = FilteringStatistics()
    for chunk in read_table_chunks(args.table, chunksize=10000):
        filtered, chunk_stats = filtered_table(chunk, v_gene_coverage=args.v_coverage,
            j_gene_coverage=args.j_coverage, v_gene_evalue=args.v_evalue)
        stats += chunk_stats
        print(filtered.to_csv(sep='\t', index=False, header=first), end='')
        first = False
        written += len(filtered)

    logger.info('%s rows in input table', stats.total)
    logger.info('%s rows have both V and J assignment', stats.has_vj_assignment)
    logger.info('%s of those do not have a stop codon', stats.has_no_stop)
    logger.info('%s of those have an E-value of at most %s', stats.good_v_evalue, args.v_evalue)
    logger.info('%s of those cover the V gene by at least %s%%', stats.good_v_coverage, args.v_coverage)
    logger.info('%s of those cover the J gene by at least %s%%', stats.good_j_coverage, args.j_coverage)
    logger.info('%d rows written', written)
    logger.info('%s rows have a recognized CDR3 (these are not filtered)', stats.has_cdr3)

    if args.json:
        with open(args.json, 'w') as f:
            json.dump(stats.asdict(), f, indent=2)
            print(file=f)
