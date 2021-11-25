"""
Plot histograms of differences to reference V gene

For each gene, a histogram is plotted that shows how often a sequence was
assigned to that gene at a certain percentage difference.
"""
import sys
import logging
import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure

from ..table import read_table

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--minimum-group-size', '-m', metavar='N', default=None, type=int,
        help="Plot only genes with at least N assigned sequences. "
        "Default: 0.1%% of assigned sequences or 100, whichever is smaller.")
    arg('--max-j-shm', metavar='VALUE', type=float, default=None,
        help='Use only rows with J%%SHM >= VALUE')
    arg('--multi', metavar='PDF', default=None,
        help='Plot individual error frequency histograms (for each V gene) to this PDF file')
    arg('--boxplot', metavar='PDF', default=None,
        help='Plot a single page with box(en)plots of V SHM for multiple genes')
    arg('table', metavar='FILTERED.TAB.GZ', help='Table with parsed IgBLAST results')


def plot_difference_histogram(group, gene_name, bins=np.arange(20.1)):
    """
    Plot a histogram of percentage differences for a specific gene.
    """
    exact_matches = group[group.V_SHM == 0]
    cdr3s_exact = len(set(s for s in exact_matches.cdr3 if s))
    js_exact = len(set(exact_matches.j_call))

    fig = Figure(figsize=(100/25.4, 60/25.4))
    ax = fig.gca()
    ax.set_xlabel('Percentage difference')
    ax.set_ylabel('Frequency')
    fig.suptitle('Gene ' + gene_name, y=1.08, fontsize=16)
    ax.set_title('{:,} sequences assigned'.format(len(group)))

    ax.text(0.25, 0.95,
        '{:,} ({:.1%}) exact matches\n  {} unique CDR3\n  {} unique J'.format(
            len(exact_matches), len(exact_matches) / len(group),
            cdr3s_exact, js_exact),
        transform=ax.transAxes, fontsize=10,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.5),
        horizontalalignment='left', verticalalignment='top')
    _ = ax.hist(list(group.V_SHM), bins=bins)
    return fig


def main(args):
    table = read_table(args.table, usecols=['v_call', 'j_call', 'V_SHM', 'J_SHM', 'cdr3'])
    if not args.multi and not args.boxplot:
        print('Don’t know what to do', file=sys.stderr)
        sys.exit(2)

    # Discard rows with any mutation within J at all
    logger.info('%s rows read', len(table))
    if args.max_j_shm is not None:
        # Discard rows with too many J mutations
        table = table[table.J_SHM <= args.max_j_shm][:]
        logger.info('%s rows remain after keeping only those with J%%SHM <= %s',
            len(table), args.max_j_shm)

    if args.minimum_group_size is None:
        total = len(table)
        minimum_group_size = min(total // 1000, 100)
        logger.info('Skipping genes with less than %s assignments', minimum_group_size)
    else:
        minimum_group_size = args.minimum_group_size

    # Genes with high enough assignment count
    all_genes = table['v_call'].unique()
    genes = sorted(table['v_call'].value_counts().loc[lambda x: x >= minimum_group_size].index)
    gene_set = set(genes)

    logger.info('%s out of %s genes have enough assignments', len(genes), len(all_genes))
    if args.multi:
        with PdfPages(args.multi) as pages:
            for gene, group in table.groupby('v_call'):
                if gene not in gene_set:
                    continue
                fig = plot_difference_histogram(group, gene)
                FigureCanvasPdf(fig).print_figure(pages, bbox_inches='tight')

        logger.info('Wrote %r', args.multi)

    if args.boxplot:
        aspect = 1 + len(genes) / 32
        g = sns.catplot(x='v_call', y='V_SHM', kind='boxen', order=genes, data=table,
            height=2 * 2.54, aspect=aspect, color='g')
        # g.set(ylim=(-.1, None))
        g.set(ylabel='% V SHM (nt)')
        g.set(xlabel='V gene')
        g.set_xticklabels(rotation=90)
        g.savefig(args.boxplot)
        logger.info('Wrote %r', args.boxplot)
