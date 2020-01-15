"""
Compare two FASTA files based on sequences

The order of records in the two files does not matter.

Exit code:
    2 if duplicate sequences or duplicate record names were found
    1 if there are any lost or gained records or sequence differences
    0 if the records are identical, but allowing for different record names

"""
import sys
import logging
import numpy as np
from scipy.optimize import linear_sum_assignment
import dnaio
from tinyalign import hamming_distance

logger = logging.getLogger(__name__)

do_not_show_cpustats = 1


def add_arguments(parser):
    arg = parser.add_argument
    arg('--color', default='auto', choices=('auto', 'never', 'always'),
        help='Whether to colorize output')
    arg('a', help='FASTA file with expected sequences')
    arg('b', help='FASTA file with actual sequences')


RED = "\x1b[0;31m"
GREEN = "\x1b[0;32m"
RESET = "\x1b[0m"


def red(s):
    return RED + s + RESET


def green(s):
    return GREEN + s + RESET


def check_duplicate_names(records):
    names = set()
    for record in records:
        if record.name in names:
            yield record.name
        names.add(record.name)


def check_exact_duplicate_sequences(records):
    sequences = dict()
    for record in records:
        if record.sequence in sequences:
            yield record.name, sequences[record.sequence]
        else:
            sequences[record.sequence] = record.name


def compare(a, b):
    """Return cost of comparing a to b"""

    l = min(len(a.sequence), len(b.sequence))
    length_diff = max(len(a.sequence), len(b.sequence)) - l
    dist_prefixes = hamming_distance(a.sequence[:l], b.sequence[:l])
    dist_suffixes = hamming_distance(a.sequence[-l:], b.sequence[-l:])

    return 5 * min(dist_prefixes, dist_suffixes) + length_diff


def pair_up_identical(a_records, b_records):
    identical = []
    b_map = {record.sequence: record for record in b_records}
    a_rest = []
    for a in a_records:
        if a.sequence in b_map:
            identical.append((a, b_map[a.sequence]))
            del b_map[a.sequence]
        else:
            a_rest.append(a)

    return identical, a_rest, list(b_map.values())


def pair_up(a_records, b_records, max_cost=20):
    # Pair up identical sequences first
    identical, a_records, b_records = pair_up_identical(a_records[:], b_records[:])

    # Compare all vs all and fill in a score matrix
    m = len(a_records)
    n = len(b_records)
    cost = np.zeros((m, n), dtype=int)
    for i, a in enumerate(a_records):
        for j, b in enumerate(b_records):
            cost[i, j] = compare(a, b)

    # Solve minimum weighted bipartite matching
    assignment = linear_sum_assignment(cost)
    similar = []
    a_similar = set()
    b_similar = set()
    for i, j in zip(*assignment):
        if cost[i, j] <= max_cost:
            similar.append((a_records[i], b_records[j]))
            a_similar.add(i)
            b_similar.add(j)

    a_only = [a for i, a in enumerate(a_records) if i not in a_similar]
    b_only = [b for j, b in enumerate(b_records) if j not in b_similar]

    return a_only, b_only, identical, similar


def format_indel(a, b, colored: bool):
    if len(a) > len(b):
        assert len(b) == 0
        s = '{-' + a + '}'
        return red(s) if colored else s
    elif len(b) > len(a):
        assert len(a) == 0
        s = '{+' + b + '}'
        return green(s) if colored else s
    else:
        return ''


def print_similar(a, b, colored: bool):
    l = min(len(a.sequence), len(b.sequence))
    dist_prefixes = hamming_distance(a.sequence[:l], b.sequence[:l])
    dist_suffixes = hamming_distance(a.sequence[-l:], b.sequence[-l:])
    if dist_prefixes <= dist_suffixes:
        a_prefix = ''
        b_prefix = ''
        a_common = a.sequence[:l]
        b_common = b.sequence[:l]
        a_suffix = a.sequence[l:]
        b_suffix = b.sequence[l:]
    else:
        a_prefix = a.sequence[:-l]
        b_prefix = b.sequence[:-l]
        a_common = a.sequence[-l:]
        b_common = b.sequence[-l:]
        a_suffix = ''
        b_suffix = ''

    s = format_indel(a_prefix, b_prefix, colored)
    edits = []
    for i, (ac, bc) in enumerate(zip(a_common, b_common)):
        if ac != bc:
            if colored:
                s = '{' + red(ac) + ' → ' + green(bc) + '}'
            else:
                s = '{' + ac + ' → ' + bc + '}'
            edits.append(s)
        else:
            edits.append(ac)
    s += ''.join(edits)

    s += format_indel(a_suffix, b_suffix, colored)

    print('~', a.name, '--', b.name)
    print(s)
    print()


def main(args):
    if args.color == 'auto':
        colored = sys.stdout.isatty()
    elif args.color == 'never':
        colored = False
    else:
        assert args.color == 'always'
        colored = True

    with dnaio.open(args.a) as f:
        a_records = list(f)
    with dnaio.open(args.b) as f:
        b_records = list(f)

    has_duplicate_names = False
    for records, path in ((a_records, args.a), (b_records, args.b)):
        dups = list(check_duplicate_names(records))
        if dups:
            has_duplicate_names = True
            print('Duplicate record names found in', path)
            for name in dups:
                print('-', name)

    has_duplicate_sequences = False
    for record, path in ((a_records, args.a), (b_records, args.b)):
        dups = list(check_exact_duplicate_sequences(records))
        if dups:
            has_duplicate_sequences = True
            print('Duplicate sequences found in', path)
            for name, name_orig in dups:
                print('-', name, 'is identical to earlier record', name_orig)

    only_a, only_b, identical, similar = pair_up(a_records, b_records)
    different_name = [(a, b) for a, b in identical if a.name != b.name]

    # Summary
    print('{} vs {} records. {} lost, {} gained, {} identical, {} different name, {} similar'.format(
        len(a_records), len(b_records), len(only_a), len(only_b),
        len(identical) - len(different_name), len(different_name),
        len(similar)))

    # Report what has changed
    if only_a:
        print()
        print('## Only in A')
        for record in only_a:
            print('-', record.name)
    if only_b:
        print()
        print('## Only in B')
        for record in only_b:
            print('+', record.name)
    if different_name:
        print()
        print('## Different name (sequence identical)')
        for a, b in different_name:
            print('=', a.name, '--', b.name)
    if similar:
        print()
        print('## Similar')
        for a, b in similar:
            print_similar(a, b, colored)

    if has_duplicate_names or has_duplicate_sequences:
        sys.exit(2)
    if only_a or only_b or similar:
        sys.exit(1)
    # different name is fine for success
    sys.exit(0)
