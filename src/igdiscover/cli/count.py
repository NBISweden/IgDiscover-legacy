"""
Compute expression counts

This command takes a table of filtered IgBLAST results (filtered.tab.gz),
filters it by various criteria and then counts how often specific
genes are named.

By default, all genes with non-zero expression are listed, sorted by
name. Use --database to change this.
"""
import sys
import logging
from typing import Optional

import dnaio
from ..table import read_table
from ..utils import natural_sort_key
from .discoverjd import filter_by_allele_ratio, compute_expressions

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--gene', choices=('V', 'D', 'J'), default='V',
        help='Which gene type: Choose V, D or J. Default: Default: %(default)s')
    arg('--database', metavar='FASTA',
        help='Compute expressions for the sequences that are named in the FASTA '
        'file. Only the sequence names in the file are used! This is the only '
        'way to also include genes with an expression of zero.')
    arg('--plot', metavar='FILE',
        help='Plot expressions to FILE (PDF or PNG)')

    group = parser.add_argument_group('Input table filtering')
    group.add_argument('--d-evalue', type=float, default=None,
        help='Maximal allowed E-value for D gene match. Default: 1E-4 '
        'if --gene=D, no restriction otherwise.')
    group.add_argument('--d-coverage', '--D-coverage', type=float, default=None,
        help='Minimum D coverage (in percent). Default: 70 '
        'if --gene=D, no restriction otherwise.')
    group.add_argument('--d-errors', type=int, default=None,
        help='Maximum allowed D errors. Default: No limit.')

    group = parser.add_argument_group('Expression counts table filtering')
    group.add_argument('--allele-ratio', type=float, metavar='RATIO', default=None,
        help='Required allele ratio. Works only for genes named "NAME*ALLELE". '
        'Default: Do not check allele ratio.')

    arg('table', help='Table with parsed and filtered IgBLAST results')


def plot_counts(counts, gene_type):
    """Plot expression counts. Return a Figure object"""
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    import matplotlib
    import seaborn as sns
    import numpy as np
    sns.set_style('white')
    sns.set_style({'patch.force_edgecolor': False})

    fig = Figure(figsize=((50 + len(counts) * 5) / 25.4, 210/25.4))
    matplotlib.rcParams.update({'font.size': 14})
    FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.set_title('{} gene usage'.format(gene_type))
    ax.set_xlabel('{} gene'.format(gene_type))
    ax.set_ylabel('Count')
    ax.set_xticks(np.arange(len(counts)))
    ax.set_xticklabels(counts.index, rotation='vertical')
    ax.grid(axis='y')
    ax.set_xlim(-1, len(counts))
    ax.bar(np.arange(len(counts)), counts['count'])
    fig.set_tight_layout(True)
    return fig


def main(args):
    run_count(**vars(args))


def run_count(
    database: Optional[str],
    table: str,
    plot: Optional[str],
    gene: str,
    d_evalue: Optional[float],
    d_coverage: Optional[float],
    d_errors: Optional[int],
    allele_ratio: Optional[float],
):
    if database and allele_ratio:
        logger.error('--database and --allele-ratio cannot be used at the same time.')
        sys.exit(1)

    if database:
        with dnaio.open(database) as fr:
            gene_names = [record.name for record in fr]
        gene_names.sort(key=natural_sort_key)
    else:
        gene_names = None

    usecols = ['v_call', 'd_call', 'j_call', 'V_errors', 'D_errors', 'J_errors', 'D_covered',
        'd_support', 'cdr3']
    table = read_table(table, usecols=usecols)
    logger.info('Table with %s rows read', len(table))

    # Set default filters depending on gene
    if gene == 'd':
        if d_evalue is None:
            d_evalue = 1E-4
        if d_coverage is None:
            d_coverage = 70
    if d_evalue is not None:
        table = table[table.D_evalue <= d_evalue]
        logger.info('%s rows remain after requiring D E-value <= %s', len(table), d_evalue)

    if d_coverage is not None:
        table = table[table.D_covered >= d_coverage]
        logger.info('%s rows remain after requiring D coverage >= %s', len(table), d_coverage)

    if d_errors is not None:
        table = table[table.D_errors <= d_errors]
        logger.info('%s rows remain after requiring D errors <= %s', len(table), d_errors)

    logger.info('Computing expression counts for %s genes', gene)
    counts = compute_expressions(table, gene)

    # Make sure that always all gene names are listed even if no sequences
    # were assigned.
    if gene_names:
        counts = counts.reindex(gene_names, fill_value=0)

    logger.info('Computed expressions for %d genes', len(counts))
    if allele_ratio is not None:
        counts = filter_by_allele_ratio(counts, allele_ratio)
        logger.info(
            'After filtering by allele ratio, %d genes remain',
            len(counts))

    # logger.info('%d sequences were not assigned a %s gene', unassigned, args.gene)
    counts.to_csv(sys.stdout, sep='\t')
    logger.info('Wrote expression count table')

    if plot:
        plot_counts(counts, gene.upper()).savefig(plot)
        logger.info("Wrote %s", plot)
