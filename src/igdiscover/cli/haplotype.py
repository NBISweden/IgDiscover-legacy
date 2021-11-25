"""
Determine haplotypes based on co-occurrences of alleles
"""
import sys
import logging
from typing import List, Tuple, Iterator
from itertools import product
from argparse import ArgumentParser
import pandas as pd
import dnaio
from ..table import read_table

logger = logging.getLogger(__name__)


# The second-most expressed allele of a gene must be expressed at at least this
# fraction of the highest-expressed allele in order for the gene to be considered
# heterozygous.
HETEROZYGOUS_THRESHOLD = 0.1

#
EXPRESSED_RATIO = 0.1


def add_arguments(parser: ArgumentParser):
    arg = parser.add_argument
    arg('--v-gene', help='V gene to use for haplotyping J. Default: Auto-detected')
    arg('--d-evalue', type=float, default=1E-4,
        help='Maximal allowed E-value for D gene match. Default: %(default)s')
    arg('--d-coverage', '--D-coverage', type=float, default=65,
        help='Minimum D coverage (in percent). Default: %(default)s%%)')
    arg('--restrict', metavar='FASTA',
        help='Restrict analysis to the genes named in the FASTA file. '
            'Only the sequence names are used!')
    arg('--order', metavar='FASTA', default=None,
        help='Sort the output according to the order of the records in '
            'the given FASTA file.')
    arg('--plot', metavar='FILE', default=None,
        help='Write a haplotype plot to FILE')
    arg('--structure-plot', metavar='FILE', default=None,
        help='Write a haplotype structure plot (counts binarized 0 and 1) to FILE')
    arg('table', help='Table with parsed and filtered IgBLAST results')


def expression_counts(table: pd.DataFrame, gene_type: str) -> Iterator[pd.DataFrame]:
    """
    Yield DataFrames for each gene with gene and allele as the row index and columns 'name' and
    'count'. For example, when 'name' is VH1-1*01, gene would be 'VH1-1' and allele
    would be '01'.
    """
    counts = table.groupby(gene_type + '_gene').size()
    names, _, alleles = zip(*[s.partition('*') for s in counts.index])
    expressions = pd.DataFrame(
        {'gene': names, 'allele': alleles, 'count': counts, 'name': counts.index},
        columns=['gene', 'allele', 'name', 'count']).set_index(['gene', 'allele'])
    del alleles

    # Example expressions at this point:
    #
    #                             name  count
    # gene       allele
    # IGHV1-18   01        IGHV1-18*01    166
    #            03        IGHV1-18*03      1
    # IGHV1-2    02         IGHV1-2*02     85
    #            04         IGHV1-2*04     16
    # IGHV1-24   01        IGHV1-24*01      5

    logger.info('Heterozygous %s genes:', gene_type)
    for _, alleles in expressions.groupby(level='gene'):
        # Remove alleles that have too low expression relative to the highest-expressed allele
        max_expression = alleles['count'].max()
        alleles = alleles[alleles['count'] >= HETEROZYGOUS_THRESHOLD * max_expression]
        if len(alleles) >= 2:
            logger.info('%s with alleles %s -- Counts: %s',
                alleles.index[0][0],
                ', '.join(alleles['name']),
                ', '.join(str(x) for x in alleles['count'])
            )
        yield alleles


class HeterozygousGene:
    def __init__(self, name: str, alleles: List[str]):
        """
        name -- name of this gene, such as 'VH4-4'
        alleles -- list of its alleles, such as ['VH4-4*01', 'VH4-4*02']
        """
        self.name = name
        self.alleles = alleles


def compute_coexpressions(table: pd.DataFrame, gene_type1: str, gene_type2: str):
    assert gene_type1 != gene_type2
    coexpressions = table.groupby(
        (gene_type1 + '_gene', gene_type2 + '_gene')).size().to_frame()
    coexpressions.columns = ['count']
    return coexpressions


def cooccurrences(coexpressions, het_alleles: Tuple[str, str], target_groups):
    """
    het_alleles -- a pair of alleles of a heterozygous gene,
    such as ('IGHJ6*02', 'IGHJ6*03').
    """
    assert len(het_alleles) == 2

    haplotype = []
    for target_alleles in target_groups:
        is_expressed_list = []
        names = []
        counts = []
        for target_allele, _ in target_alleles.itertuples(index=False):
            ex = []
            for het_allele in het_alleles:
                try:
                    e = coexpressions.loc[(het_allele, target_allele), 'count']
                except KeyError:
                    e = 0
                ex.append(e)
            ex_total = sum(ex) + 1  # +1 avoids division by zero
            ratios = [x / ex_total for x in ex]
            is_expressed = [ratio >= EXPRESSED_RATIO for ratio in ratios]
            if is_expressed != [False, False]:
                is_expressed_list.append(is_expressed)
                names.append(target_allele)
                counts.append(ex)
        if len(is_expressed_list) == 1:
            is_expressed = is_expressed_list[0]
            if is_expressed == [True, False]:
                haplotype.append((names[0], '', 'deletion', counts[0]))
            elif is_expressed == [False, True]:
                haplotype.append(('', names[0], 'deletion', counts[0]))
            elif is_expressed == [True, True]:
                haplotype.append((names[0], names[0], 'homozygous', counts[0]))
            else:
                assert False
        elif is_expressed_list == [[True, False], [False, True]]:
            haplotype.append((names[0], names[1], 'heterozygous', (counts[0][0], counts[1][1])))
        elif is_expressed_list == [[False, True], [True, False]]:
            haplotype.append((names[1], names[0], 'heterozygous', (counts[0][1], counts[1][0])))
        else:
            type_ = 'unknown'
            # Somewhat arbitrary criteria for a "duplication":
            # 1) one heterozygous allele, 2) at least three alleles in total
            n_true = sum(x.count(True) for x in is_expressed_list)
            if ([True, False] in is_expressed_list or [False, True] in is_expressed_list) and n_true > 2:
                type_ = 'duplication'
            for is_expressed, name, count in zip(is_expressed_list, names, counts):
                haplotype.append((
                    name if is_expressed[0] else '',
                    name if is_expressed[1] else '',
                    type_,
                    count,
                ))
    return haplotype


class HaplotypePair:
    """Haplotype pair for a single gene type (V/D/J)"""

    def __init__(self, haplotype, gene_type, het1, het2):
        self.haplotype = haplotype
        self.gene_type = gene_type
        self.het1 = het1
        self.het2 = het2

    def sort(self, order: List[str]) -> None:
        """
        Sort the haplotype
        order -- list a names of genes in the desired order
        """
        gene_order = {name: i for i, name in enumerate(order)}

        def keyfunc(hap):
            name = hap[0] if hap[0] else hap[1]
            gene, _, allele = name.partition('*')
            try:
                allele = int(allele)
            except ValueError:
                allele = 999
            try:
                index = gene_order[gene]
            except KeyError:
                logger.warning('Gene %s not found in gene order file, placing it at the end',
                    gene)
                index = 1000000
            return index * 1000 + allele

        self.haplotype = sorted(self.haplotype, key=keyfunc)

    def switch(self):
        """Swap the two haplotypes"""
        self.het2, self.het1 = self.het1, self.het2
        haplotype = []
        for name1, name2, type_, counts in self.haplotype:
            assert len(counts) == 2
            counts = counts[1], counts[0]
            haplotype.append((name2, name1, type_, counts))
        self.haplotype = haplotype

    def to_tsv(self, header: bool=True) -> str:
        lines = []
        if header:
            lines.append('\t'.join(['haplotype1', 'haplotype2', 'type', 'count1', 'count2']))
        lines.append(
            '# {} haplotype from {} and {}'.format(self.gene_type, self.het1, self.het2))
        for h1, h2, typ, count in self.haplotype:
            lines.append('\t'.join([h1, h2, typ, str(count[0]), str(count[1])]))
        return '\n'.join(lines) + '\n'


def plot_haplotypes(blocks: List[HaplotypePair], show_unknown: bool=False, binarize: bool=False):
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch

    colormap = dict(homozygous='cornflowerblue', heterozygous='lightgreen', deletion='gold',
        duplication='crimson', unknown='gray')
    if not show_unknown:
        del colormap['unknown']
    labels = [[], []]
    heights = [[], []]
    names = []
    colors = []
    positions = []
    pos = 0
    for block in blocks:
        for hap1, hap2, type_, (count1, count2) in block.haplotype:
            if not show_unknown and type_ == 'unknown':
                continue
            if binarize:
                count1 = 1 if hap1 else 0
                count2 = 1 if hap2 else 0
            label = hap1 if hap1 else hap2
            hap1 = ('*' + hap1.partition('*')[2]) if hap1 else ''
            hap2 = ('*' + hap2.partition('*')[2]) if hap2 else ''
            heights[0].append(count1)
            heights[1].append(count2)
            labels[0].append(hap1)
            labels[1].append(hap2)
            names.append(label.partition('*')[0])
            colors.append(colormap[type_])
            positions.append(pos)
            pos += 1
        pos += 2

    n = len(names)
    assert len(labels[0]) == len(labels[1]) == len(colors) == len(heights[0]) == len(heights[1]) == n

    fig = Figure(figsize=(8, 24))
    FigureCanvas(fig)
    axes = fig.subplots(ncols=2)
    for i in 0, 1:
        axes[i].barh(y=positions, width=heights[i], color=colors)
        axes[i].set_yticklabels(labels[i])
        axes[i].set_ylim((-0.5, max(positions) + 0.5))

    # Add center axis
    ax_center = axes[0].twinx()
    ax_center.set_yticks(axes[0].get_yticks())
    ax_center.set_ylim(axes[0].get_ylim())
    ax_center.set_yticklabels(names)

    # Synchronize x limits on both axes (has no effect if binarize is True)
    max_x = max(axes[i].get_xlim()[1] for i in range(2))
    for ax in axes:
        ax.set_xlim(right=max_x)
    axes[0].invert_xaxis()

    for ax in axes[0], axes[1], ax_center:
        ax.invert_yaxis()
        for spine in ax.spines.values():
            spine.set_visible(False)
        if binarize:
            ax.set_xticks([])
        else:
            ax.set_axisbelow(True)
            ax.grid(True, axis='x')
        ax.set_yticks(positions)
        ax.tick_params(left=False, right=False, labelsize=12)
    axes[1].tick_params(labelleft=False, labelright=True)
    if not binarize:
        axes[0].spines['right'].set_visible(True)
        axes[1].spines['left'].set_visible(True)

    # Legend
    legend_patches = [Patch(color=col, label=label) for label, col in colormap.items()]
    fig.legend(handles=legend_patches, loc='upper center', ncol=len(legend_patches))
    fig.tight_layout()
    fig.subplots_adjust(top=1 - 2 / len(names))

    return fig


def read_and_filter(path: str, d_evalue: float, d_coverage: float):
    usecols = ['v_call', 'd_call', 'j_call', 'V_errors', 'D_errors', 'J_errors', 'D_covered',
        'd_support']
    # Support reading a table without D_errors
    try:
        table = read_table(path, usecols=usecols)
    except ValueError:
        usecols.remove('D_errors')
        table = read_table(path, usecols=usecols)
    logger.info('Table with %s rows read', len(table))

    table = table[table.V_errors == 0]
    logger.info('%s rows remain after requiring V errors = 0', len(table))
    table = table[table.J_errors == 0]
    logger.info('%s rows remain after requiring J errors = 0', len(table))
    table = table[table.d_support <= d_evalue]
    logger.info('%s rows remain after requiring D E-value <= %s', len(table), d_evalue)
    table = table[table.D_covered >= d_coverage]
    logger.info('%s rows remain after requiring D coverage >= %s', len(table), d_coverage)
    if 'D_errors' in table.columns:
        table = table[table.D_errors == 0]
        logger.info('%s rows remain after requiring D errors = 0', len(table))
    return table


def main(args):
    if args.order is not None:
        with dnaio.open(args.order) as sr:
            gene_order = [r.name for r in sr]
    else:
        gene_order = None

    table = read_and_filter(args.table, args.d_evalue, args.d_coverage)

    if args.restrict is not None:
        with dnaio.open(args.restrict) as sr:
            restrict_names = set(r.name for r in sr)
        table = table[table['v_call'].map(lambda name: name in restrict_names)]
        logger.info('After restricting to V genes named in %r, %d rows remain', args.restrict,
            len(table))
        if len(table) == 0:
            logger.error('No rows remain, cannot continue')
            sys.exit(1)

    expressions = dict()
    het_expressions = dict()  # these are also sorted, most highly expressed first
    for gene_type in 'VDJ':
        ex = list(expression_counts(table, gene_type))
        expressions[gene_type] = ex
        het_ex = [e for e in ex if len(e) == 2]
        if het_ex:
            # Pick most highly expressed
            het_expressions[gene_type] = sorted(het_ex, key=lambda e: e['count'].sum(), reverse=True)[:5]
        else:
            # Force at least something to be plotted
            het_expressions[gene_type] = [None]

    if args.v_gene:
        het_ex = [e for e in expressions['V'] if len(e) == 2]
        for ex in het_ex:
            if (args.v_gene in ex.index and not ex.loc[args.v_gene].empty) or (args.v_gene in ex['name'].values):
                het_expressions['V'] = [ex]
                break
        else:
            logger.error('The gene or allele %s was not found in the list of heterozygous V genes. '
                'It cannot be used with the --v-gene option.', args.v_gene)
            sys.exit(1)

    block_lists = []

    # We want to avoid using a gene classifed as 'duplicate' for haplotyping, but the
    # classification is only known after we have done the haplotyping, so we try it
    # until we found a combination that works
    products = list(product(het_expressions['J'], het_expressions['V']))
    for attempt, (het_j, het_v) in enumerate(products):
        best_het_genes = {
            'V': het_v,
            'D': het_expressions['D'][0] if het_expressions['D'] else None,
            'J': het_j,
        }
        for gene_type in 'VDJ':
            bhg = best_het_genes[gene_type]
            text = bhg.index[0][0] if bhg is not None else 'none found'
            logger.info('Heterozygous %s gene to use for haplotyping: %s',
                gene_type, text)

        # Create HaplotypePair objects ('blocks') for each gene type
        blocks = []

        for target_gene_type, het_gene in (
            ('J', 'V'),
            ('D', 'J'),
            ('V', 'J'),
        ):
            het_alleles = best_het_genes[het_gene]
            if het_alleles is None:
                continue
            coexpressions = compute_coexpressions(table, het_gene, target_gene_type)
            target_groups = expressions[target_gene_type]
            het1, het2 = het_alleles['name']
            haplotype = cooccurrences(coexpressions, (het1, het2), target_groups)
            block = HaplotypePair(haplotype, target_gene_type, het1, het2)
            if gene_order:
                block.sort(gene_order)
            blocks.append(block)

        if het_j is None or het_v is None:
            break
        het_used = set(sum([list(h['name']) for h in best_het_genes.values() if h is not None], []))
        het_is_duplicate = False
        for block in blocks:
            if block.gene_type == 'D':
                continue

            # This nested loop is put in a separate generator so we can easily 'break'
            # out of both loops at the same time.
            def nameiter():
                for name1, name2, type_, _ in block.haplotype:
                    for name in name1, name2:
                        yield name, type_

            for name, type_ in nameiter():
                if name in het_used and type_ != 'heterozygous':
                    het_is_duplicate = True
                    if not args.v_gene:
                        logger.warning('%s not classified as "heterozygous" during haplotyping, '
                            'attempting to use different alleles', name)
                    break
            if het_is_duplicate:
                break
        block_lists.append(blocks)
        if not het_is_duplicate:
            break
    else:
        if not args.v_gene:
            logger.warning('No other alleles remain, using first found solution')
        blocks = block_lists[0]

    # Get the phasing right across blocks (i.e., swap J haplotypes if necessary)
    assert len(blocks) in (0, 1, 3)
    if len(blocks) == 3:
        j_hap, d_hap, v_hap = blocks
        assert j_hap.gene_type == 'J'
        assert d_hap.gene_type == 'D'
        assert v_hap.gene_type == 'V'
        assert d_hap.het1 == v_hap.het1
        assert d_hap.het2 == v_hap.het2
        for name1, name2, _, _ in j_hap.haplotype:
            if (name1, name2) == (v_hap.het2, v_hap.het1):
                j_hap.switch()
                break

    # Print the table
    header = True
    for block in blocks:
        print(block.to_tsv(header=header))
        header = False

    # Create plots if requested
    if args.plot:
        fig = plot_haplotypes(blocks, show_unknown=True)
        fig.savefig(args.plot)
    if args.structure_plot:
        fig = plot_haplotypes(blocks, binarize=True)
        fig.savefig(args.structure_plot)
