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
    >>> describe_nt_change('TTGGTCAAGCCTGGAGGGTCCCTGAGACTCTCCTGTGCAGCCTCT', 'GGCTCCCAGATACTCTCCTCTGCAGCCTCT')
    '1_15delTTGGTCAAGCCTGGA; 18G>C; 23T>A; 26G>T; 35G>C'
    """
    alignment = align_affine(s, t)
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
    for event, group in groupby(zip(alignment.ref_row, alignment.query_row), grouper):
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


def align_affine(
    ref: str, query: str, gap_open=-6, gap_extend=-1, mismatch=-3, match=1
):
    """
    Return an optimal global alignment of strings ref and query, using affine gap penalties.

    The default scores are those that BWA uses.
    """
    #          query (n, j)
    #        --------------->
    #       |
    #   ref | (m, i)
    #       |
    #       V
    #
    m = len(ref)
    n = len(query)

    if not (gap_open <= 0 and gap_extend <= 0 and mismatch <= 0 and match >= 0):
        raise ValueError(
            "gap_open/gap_extend/mismatch scores must be <= 0; "
            "match score must be >= 0"
        )

    # Initialize three dynamic programming tables:
    #
    # - M(i, j) is the score of an optimal alignment between ref[:i] and query[:j]
    #   that does not end with an insertion nor a deletion
    # - H(i, j) is the score of an optimal alignment between ref[:i] and query[:j]
    #   that aligns query[j-1] to a gap
    # - V(i, j) is the score of an optimal alignment between ref[:i] and query[:j]
    #   that aligns ref[i-1] to a gap

    def zero_matrix():
        return [[0] * (n + 1) for _ in range(m + 1)]

    M = zero_matrix()
    H = zero_matrix()
    V = zero_matrix()
    gapchar = "-"

    inf = float("inf")
    H[0][0] = V[0][0] = -inf

    # initialize first column
    for i in range(1, m + 1):
        M[i][0] = -inf
        H[i][0] = -inf
        V[i][0] = gap_open + (i - 1) * gap_extend

    # initialize first row
    for j in range(1, n + 1):
        M[0][j] = -inf
        H[0][j] = gap_open + (j - 1) * gap_extend
        V[0][j] = -inf

    for i in range(1, m + 1):
        refchar = ref[i - 1]
        for j in range(1, n + 1):
            diag_score = match if refchar == query[j - 1] else mismatch
            M[i][j] = diag_score + max(
                M[i - 1][j - 1], V[i - 1][j - 1], H[i - 1][j - 1]
            )
            # The terms where V is updated from H and H updated from V can be
            # removed if one thinks that insertions should not follow indels
            # directly.
            V[i][j] = max(
                H[i - 1][j] + gap_open, V[i - 1][j] + gap_extend, M[i - 1][j] + gap_open
            )
            H[i][j] = max(
                H[i][j - 1] + gap_extend, V[i][j - 1] + gap_open, M[i][j - 1] + gap_open
            )

    ref_row = []
    query_row = []
    i = m
    j = n
    errors = 0

    TRACE_M, TRACE_V, TRACE_H = "MVH"
    optimal_score = max(M[m][n], V[m][n], H[m][n])
    if optimal_score == M[m][j]:
        state = TRACE_M
    elif optimal_score == V[m][j]:
        state = TRACE_V
    else:
        state = TRACE_H

    # Compute traceback
    while i > 0 or j > 0:
        assert i >= 0 and j >= 0
        if state == TRACE_M:
            m_score = M[i - 1][j - 1]
            v_score = V[i - 1][j - 1]
            h_score = H[i - 1][j - 1]
            c1 = ref[i - 1]
            c2 = query[j - 1]
            ref_row.append(c1)
            query_row.append(c2)
            if c1 != c2:
                errors += 1
            i -= 1
            j -= 1
        elif state == TRACE_V:
            m_score = M[i - 1][j] + gap_open
            v_score = V[i - 1][j] + gap_extend
            h_score = H[i - 1][j] + gap_open
            errors += 1
            ref_row.append(ref[i - 1])
            query_row.append(gapchar)
            i -= 1
        else:
            assert state == TRACE_H
            m_score = M[i][j - 1] + gap_open
            v_score = V[i][j - 1] + gap_open
            h_score = H[i][j - 1] + gap_extend
            errors += 1
            ref_row.append(gapchar)
            query_row.append(query[j - 1])
            j -= 1

        # Update state
        if h_score > m_score and h_score > v_score:
            state = TRACE_H
        elif v_score > m_score:
            state = TRACE_V
        else:
            state = TRACE_M

    ref_row, query_row = "".join(ref_row[::-1]), "".join(query_row[::-1])
    return Alignment(ref_row, query_row, 0, m, 0, n, optimal_score, errors)
