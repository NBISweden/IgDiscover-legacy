"""
Dereplicate sequences and remove barcodes

The difference to the 'group' subcommand is that that one
also computes consensus sequences from groups with identical
barcode/CDR3. This one does not.

The barcode can be in the 5' end or the 3' end of the sequence.

Use --trim-g to remove initial runs of G at the 5' end (artifact from RACE protocol).
These are removed after the barcode is removed.
"""

import logging
from collections import defaultdict
from itertools import islice
import json
import dnaio

logger = logging.getLogger(__name__)


def add_arguments(parser):
    arg = parser.add_argument
    arg('--limit', default=None, type=int, metavar='N',
        help='Limit processing to the first N reads')
    arg('--trim-g', action='store_true', default=False,
        help="Trim 'G' nucleotides at 5' end")
    arg('--minimum-length', '-l', type=int, default=0,
        help='Minimum sequence length')
    arg('--barcode-length', '-b', type=int, default=0,
        help="Length of barcode. Positive for 5' barcode, negative for 3' barcode. "
             "Default: %(default)s")
    arg('--json', metavar="FILE", help="Write statistics to FILE")
    arg('fastx', metavar='FASTA/FASTQ',
        help='FASTA or FASTQ file (can be gzip-compressed) with sequences')


def main(args):
    barcode_length = args.barcode_length
    too_short = 0
    n = 0
    sequences = defaultdict(list)  # maps sequences to a list of Sequence objects containing them
    with dnaio.open(args.fastx) as f:
        for record in islice(f, 0, args.limit):
            n += 1
            if len(record) < args.minimum_length:
                too_short += 1
                continue
            sequences[record.sequence].append(record)

    n_written = 0
    for records in sequences.values():
        # If there are multiple records with the same sequence, pick the first
        record = records[0]

        if barcode_length >= 0:
            barcode = record.sequence[:barcode_length]
            unbarcoded = record[barcode_length:]
        else:
            barcode = record.sequence[barcode_length:]
            unbarcoded = record[:barcode_length]

        if args.trim_g:
            # The RACE protocol leads to a run of non-template Gs in the beginning
            # of the sequence, after the barcode.
            unbarcoded.sequence = unbarcoded.sequence.lstrip('G')
            if unbarcoded.qualities:
                unbarcoded.qualities = unbarcoded.qualities[-len(unbarcoded.sequence):]

        name = record.name.split(maxsplit=1)[0]
        if name.endswith(';'):
            name = name[:-1]

        if barcode_length:
            print('>{};barcode={};size={};\n{}'.format(name, barcode, len(records),
                unbarcoded.sequence))
        else:
            print('>{};size={};\n{}'.format(name, len(records),
                unbarcoded.sequence))

        n_written += 1

    logger.info('%s sequences processed', n)
    logger.info('%s sequences long enough', n - too_short)
    logger.info('%s dereplicated sequences written', n_written)

    if args.json:
        stats = {
            'groups_written': n_written,
        }
        with open(args.json, 'w') as f:
            json.dump(stats, f, indent=2)
            print(file=f)
