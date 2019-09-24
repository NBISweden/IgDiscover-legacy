import subprocess
from collections import Counter, OrderedDict
from itertools import groupby
from io import BytesIO

import dnaio

from .utils import available_cpu_count, downsampled


def multialign(sequences, program="mafft", threads=available_cpu_count()):
    """
    Wrapper for multiple sequence alignment tools. Currently supported are
    * ClustalO (http://www.clustal.org/omega/),
    * MAFFT (http://mafft.cbrc.jp/alignment/software/)
    * MUSCLE (http://www.drive5.com/muscle/)

    A package using the libclustalo library directly also exists:
    https://github.com/benchling/clustalo-python/blob/master/clustalo.c
    It is not for Python 3.

    sequences -- a dictionary mapping names to sequences.
        Use an OrderedDict if order matters.

    program -- must be "clustalo", "mafft", "muscle", "muscle-medium",
        "muscle-fast". The latter calls MUSCLE with parameters that make it run
        faster (but less accurate). "muscle-medium" is in between muscle and
        muscle-fast.

    threads -- number of threads to use for those programs that support it.
        By default, set to the number of processors.
    """
    if program == "mafft":
        args = ["mafft", "--quiet", "--thread", str(threads), "-"]
    elif program == "clustalo":
        args = ["clustalo", "--threads=" + str(threads), "--infile=-"]
    elif program == "muscle":
        args = ["muscle", "-quiet", "-in", "-", "-out", "-"]
    elif program == "muscle-fast":
        args = ["muscle", "-quiet", "-maxiters", "1", "-diags", "-in", "-", "-out", "-"]
    elif program == "muscle-medium":
        args = ["muscle", "-quiet", "-maxiters", "2", "-diags", "-in", "-", "-out", "-"]
    else:
        raise ValueError("program {!r} not supported".format(program))

    fasta_data = "".join(
        ">{}\n{}\n".format(name, seq) for name, seq in sequences.items()
    ).encode()
    result = subprocess.check_output(args, input=fasta_data)
    with dnaio.open(BytesIO(result)) as f:
        aligned = list(f)

    return {record.name: record.sequence.upper() for record in aligned}


def consensus(aligned, threshold=0.7, ambiguous='N', keep_gaps=False):
    """
    Compute a consensus from multialign() output, allowing degraded sequences
    in the 3' end.

    aligned -- a dict mapping names to sequences or a list of sequences
    keep_gaps -- whether the returned sequence contains gaps (-)
    """
    n = len(aligned)
    result = []
    if hasattr(aligned, 'values'):
        sequences = aligned.values()
    else:
        sequences = aligned

    active = int(len(aligned) * 0.05)
    for i, chars in enumerate(reversed(list(zip(*sequences)))):
        counter = Counter(chars)
        active = max(n - counter['-'], active)
        assert counter['-'] >= n - active
        counter['-'] -= n - active
        char, freq = counter.most_common(1)[0]
        if i >= 10:  # TODO hard-coded
            active = n
        if freq / active >= threshold:
            if keep_gaps or char != '-':
                result.append(char)
        else:
            result.append(ambiguous)
    return ''.join(result[::-1])


def iterative_consensus(sequences, program='muscle-medium', threshold=0.6,
        subsample_size=200, maximum_subsample_size=1600):
    """
    Compute a consensus sequence of the given sequences, but do not use all
    sequences if there are many: First, try to compute the consensus from a
    small subsample. If there are 'N' bases, increase the subsample size and
    repeat until either there are no more 'N' bases, all available sequences
    have been used or maximum_subsample_size is reached.
    """
    while True:
        sample = downsampled(sequences, subsample_size)
        aligned = multialign(OrderedDict(enumerate(sample)), program=program)

        cons = consensus(aligned, threshold=threshold).strip('N')
        if 'N' not in cons:
            # This consensus is good enough
            break
        if len(sequences) <= subsample_size:
            # We have already used all the sequences that are available
            break
        subsample_size *= 2
        if subsample_size > maximum_subsample_size:
            break
    return cons


def describe_nt_change(s: str, t: str):
    """
    Describe changes between two nucleotide sequences

    >>> describe_nt_change('AAA', 'AGA')
    '2A>G'
    >>> describe_nt_change('AAGG', 'AATTGG')
    '2_3insTT'
    >>> describe_nt_change('AATTGG', 'AAGG')
    '3_4delTT'
    >>> describe_nt_change('AATTGGCG', 'AAGGTG')
    '3_4delTT; 7C>T'
    >>> describe_nt_change('AAGCTAA', 'AACTGAA')
    '3G>C; 4C>T; 5T>G'
    """
    alignment = align_global(s, t, insertion=-4, deletion=-4)
    row1 = alignment.ref_row.replace('\0', '-')
    row2 = alignment.query_row.replace('\0', '-')
    changes = []

    def grouper(c):
        c1, c2 = c
        if c1 == c2:
            return 'MATCH'
        elif c1 == '-':
            return 'INS'
        elif c2 == '-':
            return 'DEL'
        else:
            return 'SUBST'

    index = 1
    for event, group in groupby(zip(row1, row2), grouper):
        if event == 'MATCH':
            index += len(list(group))
        elif event == 'SUBST':
            # ungroup
            for c1, c2 in group:
                change = '{}{}>{}'.format(index, c1, c2)
                changes.append(change)
                index += 1
        elif event == 'INS':
            inserted = ''.join(c[1] for c in group)
            change = '{}_{}ins{}'.format(index-1, index, inserted)
            changes.append(change)
        elif event == 'DEL':
            deleted = ''.join(c[0] for c in group)
            change = '{}_{}del{}'.format(index, index + len(deleted) - 1, deleted)
            changes.append(change)
            index += len(deleted)
    return '; '.join(changes)


class Alignment:
    def __init__(self, ref_row, query_row, ref_start, ref_stop, query_start, query_stop, score, errors):
        assert len(ref_row) == len(query_row)
        self.ref_row = ref_row
        self.ref_start = ref_start
        self.ref_stop = ref_stop
        self.query_row = query_row
        self.query_start = query_start
        self.query_stop = query_stop
        self.score = score
        self.errors = errors


def align_global_core(ref: str, query: str, match=1, mismatch=-2, insertion=-2, deletion=-2):
    """
    Compute an optimal global alignment between strings *ref* and *query*.
    An alignment is optimal if it has maximal score.

    Return an Alignment object.
    """
    m = len(ref)
    n = len(query)
    # DP Matrix:
    #            query (j)
    #          ----------> n
    #         |
    # ref (i) |
    #         |
    #         V
    #         m

    # the DP matrix is stored column-major
    scores = [[0] * (n + 1) for _ in range(m+1)]
    backtrace = [[0] * (n + 1) for _ in range(m+1)]

    start_in_ref = False
    start_in_query = False
    stop_in_ref = False
    stop_in_query = False

    UP = 0
    LEFT = 1
    DIAG = 2
    GAPCHAR = '\0'

    # initialize first column
    for i in range(m + 1):
        scores[i][0] = 0 if start_in_ref else i * deletion
        backtrace[i][0] = UP

    # initialize first row
    for j in range(n + 1):
        scores[0][j] = 0 if start_in_query else j * insertion
        backtrace[0][j] = LEFT

    # fill the entire DP matrix
    # outer loop goes over columns
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            bt = DIAG
            score = scores[i-1][j-1] + (match if (ref[i-1] == query[j-1]) else mismatch)
            tmp = scores[i-1][j] + insertion
            if tmp > score:
                bt = UP
                score = tmp
            tmp = scores[i][j-1] + deletion
            if tmp > score:
                bt = LEFT
                score = tmp
            scores[i][j] = score
            backtrace[i][j] = bt

    # initialize best score and its position to the bottomright cell
    best_i = m
    best_j = n
    best_score = scores[m][n]

    if stop_in_query:
        # search also in last row
        for j in range(n + 1):
            if scores[m][j] >= best_score:
                best_score = scores[m][j]
                best_i = m
                best_j = j

    if stop_in_ref:
        # search also in last column
        for i in range(m + 1):
            if scores[i][n] >= best_score:
                best_i = i
                best_j = n
                best_score = scores[i][n]

    # trace back
    p1 = []
    p2 = []
    i = m
    j = n

    # first, walk from the lower right corner to the
    # position where we found the maximum score

    errors = 0

    # if gaps are currently errors, this is 1, otherwise it's 0
    gaps_are_errors = 0 if stop_in_query else 1
    if i == best_i:  # we are in the last row
        while j > best_j:
            p1.append(GAPCHAR)
            j -= 1
            p2.append(query[j])
            errors += gaps_are_errors
    else:  # we are in the last column
        gaps_are_errors = 0 if stop_in_ref else 1
        while i > best_i:
            i -= 1
            p1.append(ref[i])
            p2.append(GAPCHAR)
            errors += gaps_are_errors

    assert i == best_i and j == best_j

    # the actual backtracing
    # The alignments are constructed in reverse
    # and this is undone afterwards.
    while i > 0 and j > 0:
        direction = backtrace[i][j]
        if direction == DIAG:
            i -= 1
            j -= 1
            if ref[i] != query[j]:
                errors += 1
            p1.append(ref[i])
            p2.append(query[j])
        elif direction == LEFT:
            errors += 1
            p1.append(GAPCHAR)
            j -= 1
            p2.append(query[j])
        elif direction == UP:
            i -= 1
            p1.append(ref[i])
            p2.append(GAPCHAR)
            errors += 1
        else:
            assert False, "DP table corrupt"

    start1 = i if start_in_ref else 0
    start2 = j if start_in_query else 0

    errors += (i - start1) + (j - start2)

    while j > 0:
        p1.append(GAPCHAR)
        j -= 1
        p2.append(query[j])
    while i > 0:
        i -= 1
        p1.append(ref[i])
        p2.append(GAPCHAR)
    assert i == 0 and j == 0

    p1 = "".join(p1[::-1])
    p2 = "".join(p2[::-1])

    return Alignment(p1, p2, start1, best_i, start2, best_j, best_score, errors)


def align_global(ref: str, query: str, match=1, mismatch=-2, insertion=-2, deletion=-2):
    m = len(ref)
    n = len(query)
    start = 0
    while start < m and start < n and ref[start] == query[start]:
        start += 1

    stop_ref = m
    stop_query = n
    while stop_ref > start and stop_query > start and ref[stop_ref-1] == query[stop_query-1]:
        stop_ref -= 1
        stop_query -= 1

    assert 0 <= start <= stop_ref <= m
    assert 0 <= start <= stop_query <= n

    alignment = align_global_core(
        ref[start:stop_ref], query[start:stop_query], match, mismatch, insertion, deletion)

    alignment.ref_row = ref[:start] + alignment.ref_row + ref[stop_ref:]
    assert alignment.ref_start == 0
    assert alignment.ref_stop == stop_ref - start
    alignment.ref_stop = m
    alignment.query_row = query[:start] + alignment.query_row + query[stop_query:]
    assert alignment.query_start == 0
    assert alignment.query_stop == stop_query - start
    alignment.query_stop = n
    alignment.score += start + m - stop_ref

    return alignment
