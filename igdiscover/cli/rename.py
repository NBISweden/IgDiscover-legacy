"""
Rename and reorder records in a FASTA file

Sequences can be renamed according to the sequences in a template file.
The sequences in the target file will get the name that they have in the
template file. Sequences are considered to be equivalent if one is a prefix of
the other.

Sequences can also be ordered by name, either alphabetically or by
the order given in a template file. For comparison, only the 'gene part'
of the name is used (before the '*'). For example, for 'VH4-4*01', the name 'VH4-4'
is looked up in the template. Alphabetic order is the default. Use --no-sort
to disable sorting entirely.
"""
import sys
import logging
import dnaio

from ..utils import natural_sort_key

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--no-sort', dest='sort', action='store_false', default=True,
        help='Do not sort sequences by name')
    arg('--not-found', metavar='TEXT', default=' (not found)',
        help='Append this text to the record name when the sequence was not found '
        'in the template. Default: %(default)r')
    arg('--rename-from', metavar='TEMPLATE',
        help='FASTA template file with correctly named sequences. If a sequence in '
            'the target file is identical to one in the template, it is assigned the '
            'name of the sequence in the template.')
    arg('--order-by', metavar='TEMPLATE',
        help='FASTA template that contains genes in the desired order. '
            'If a name contains a "*" (asterisk), only the preceding part '
            'is used. Thus, VH4-4*01 and VH4-4*02 are equivalent.')
    arg('target', help='FASTA file with to-be renamed sequences')


class PrefixDict:
    """
    A dict that maps strings to values, but where a prefix of a key is enough
    to retrieve the value.
    """
    def __init__(self, items):
        self._items = []
        for k, v in items:
            self.add(k, v)

    def add(self, k, v):
        try:
            self[k]
            raise ValueError('Key {!r} already exists'.format(k))
        except KeyError:
            self._items.append((k, v))

    def __getitem__(self, key):
        found = None
        for seq, value in self._items:
            if seq.startswith(key) or key.startswith(seq):
                if found is not None:
                    # TODO don't use keyerror here
                    raise KeyError('Key {!r} is ambiguous'.format(key))
                found = value
        if found is None:
            raise KeyError(key)
        return found

    def get(self, key, default=None):
        try:
            v = self[key]
        except KeyError:
            return default
        else:
            return v

    def __len__(self):
        return len(self._items)


class GeneMissing(Exception):
    pass


def gene_name(record):
    record_name, _, record_comment = record.name.partition(' ')
    return record_name.partition('*')[0]


def sorted_by_gene(records, gene_order):
    d = {name: i for i, name in enumerate(gene_order)}

    def keyfunc(record):
        gene = gene_name(record)
        try:
            return d[gene]
        except KeyError:
            raise GeneMissing(gene)

    return sorted(records, key=keyfunc)


def main(args):
    if args.rename_from:
        with dnaio.open(args.rename_from) as fr:
            template = PrefixDict([])
            for record in fr:
                try:
                    template.add(record.sequence.upper(), record.name)
                except ValueError:
                    logger.error('Sequences in entry %r and %r are duplicate',
                        record.name, template[record.sequence.upper()])
        logger.info('Read %d entries from template', len(template))
    else:
        template = None

    if args.order_by:
        with dnaio.open(args.order_by) as fr:
            gene_order = [gene_name(r) for r in fr]
    else:
        gene_order = None

    with dnaio.open(args.target) as fr:
        sequences = list(fr)

    # Rename
    renamed = 0
    if template is not None:
        for record in sequences:
            name = template.get(record.sequence.upper())
            if name is None:
                name = record.name + args.not_found
            else:
                renamed += 1
            # Replace record’s name, leaving comment intact
            record_name, _, record_comment = record.name.partition(' ')
            if record_comment:
                record.name = name + ' ' + record_comment
            else:
                record.name = name

    # Reorder
    if gene_order:
        try:
            sequences = sorted_by_gene(sequences, gene_order)
        except GeneMissing as e:
            logger.error('Gene "%s" not found in the --order-by template file', e)
            sys.exit(1)
    elif args.sort:
        sequences = sorted(sequences, key=lambda s: natural_sort_key(s.name))
    for record in sequences:
        print('>{}\n{}'.format(record.name, record.sequence))
    logger.info('Wrote %s FASTA records (%d sequences found in template)', len(sequences), renamed)
