import subprocess
from collections import Counter, OrderedDict
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
