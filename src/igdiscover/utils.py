"""
Some utility functions that work on sequences and lists of sequences.
"""
import os
import sys
import re
import random
import hashlib
import resource
from collections import defaultdict
from typing import List

import dnaio
import numpy as np
from tinyalign import edit_distance
from cutadapt.align import Aligner

from .dna import GENETIC_CODE


def downsampled(population, size):
    """
    Return a random subsample of the population.

    Uses reservoir sampling. See https://en.wikipedia.org/wiki/Reservoir_sampling
    """
    sample = population[:size]
    for index in range(size, len(population)):
        r = random.randint(0, index)
        if r < size:
            sample[r] = population[index]
    return sample


def distances(sequences, band=0.2):
    """
    Compute all pairwise edit distances and return a square matrix.

    Entry [i,j] in the matrix is the edit distance between sequences[i]
    and sequences[j].
    """
    # Pre-compute distances between unique sequences
    unique_sequences = list(set(sequences))
    unique_distances = dict()  # maps (seq1, seq2) tuples to edit distance
    maxdiff = max((int(len(s) * band) for s in sequences), default=0)  # TODO double-check this
    for i, s in enumerate(unique_sequences):
        for j, t in enumerate(unique_sequences):
            if i < j:
                dist = min(maxdiff+1, edit_distance(s, t, maxdiff=maxdiff))
                unique_distances[(t, s)] = dist
                unique_distances[(s, t)] = dist

    # Fill the result matrix
    m = np.zeros((len(sequences), len(sequences)), dtype=float)
    for i, s in enumerate(sequences):
        for j, t in enumerate(sequences):
            if i < j:
                d = 0 if s == t else unique_distances[(s, t)]
                m[j, i] = m[i, j] = d
    return m


def sequence_hash(s, digits=4):
    """
    For a string, return a 'fingerprint' that looks like 'S1234' (the character
    'S' is fixed). The idea is that this allows one to quickly see whether two
    sequences are not identical.
    """
    h = int(hashlib.md5(s.encode()).hexdigest()[-4:], base=16)
    return 'S' + str(h % 10**digits).rjust(digits, '0')


def unique_name(name, sequence):
    """Create a unique name based on the current name and the sequence

    The returned name looks like name_S1234. If the current name contains
    already a S_.... suffix, it is removed before the new suffix is appended.

    name -- current name
    """
    return '{}_{}'.format(name.rsplit('_S', 1)[0], sequence_hash(sequence))


class SerialPool:
    """
    An alternative to multiprocessing.Pool that runs things in serial for
    easier debugging
    """
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def imap(self, func, iterable, chunksize):
        for i in iterable:
            yield func(i)


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    """
    Use this function as sorting key to sort in 'natural' order.

    >>> names = ['file10', 'file1.5', 'file1']
    >>> sorted(names, key=natural_sort_key)
    ['file1', 'file1.5', 'file10']

    Source: http://stackoverflow.com/a/16090640/715090
    """
    return [int(text) if text.isdigit() else text.lower()
        for text in re.split(_nsre, s)]


class UniqueNamer:
    """
    Assign unique names by appending letters to already seen names.
    """
    def __init__(self):
        self._names = set()

    def __call__(self, name):
        ext = 'A'
        new_name = name
        while new_name in self._names:
            if ext == '[':
                raise ValueError('Too many duplicate names')
            new_name = name + ext
            ext = chr(ord(ext) + 1)
        self._names.add(new_name)
        return new_name


class Merger:
    """
    Merge similar items into one. To specify what "similar" means, implement
    the merged() method in a subclass.
    """
    def __init__(self):
        self._items = []

    def add(self, item):
        # This method could possibly be made simpler if the graph structure
        # was made explicit.

        items = []
        for existing_item in self._items:
            m = self.merged(existing_item, item)
            if m is None:
                items.append(existing_item)
            else:
                item = m
        items.append(item)
        self._items = items

    def extend(self, iterable):
        for i in iterable:
            self.add(i)

    def __iter__(self):
        if self._items and hasattr(self._items, 'name'):
            yield from sorted(self._items, key=lambda x: x.name)
        else:
            yield from self._items

    def __len__(self):
        return len(self._items)

    def merged(self, existing_item, item):
        """
        If existing_item and item can be returned, this method must return
        a new item that represents a merged version of both. If they cannot
        be merged, it must return None.
        """
        raise NotImplementedError("not implemented")


def relative_symlink(src, dst, force=False):
    """
    Create a symbolic link in any directory.

    force -- if True, then overwrite an existing file/symlink
    """
    if force:
        try:
            os.remove(dst)
        except FileNotFoundError:
            pass
    target = os.path.relpath(os.path.abspath(src), start=os.path.dirname(dst))
    os.symlink(target, dst)


def nt_to_aa(s, _get=GENETIC_CODE.get):
    """Translate a nucleotide sequence to an amino acid sequence"""
    return ''.join([_get(s[i:i+3], '*') for i in range(0, len(s), 3)])


def has_stop(sequence):
    """
    Return a boolean indicating whether the sequence has an internal stop codon.
    An incomplete codon at the end is allowed.

    >>> has_stop('GGG')
    False
    >>> has_stop('TAA')
    True
    >>> has_stop('GGGAC')
    False
    """
    s = sequence[:len(sequence) // 3 * 3]
    return '*' in nt_to_aa(s)


def plural_s(n):
    return 's' if n != 1 else ''


class FastaValidationError(Exception):
    pass


def validate_fasta(path):
    """
    Ensure that the FASTA file is suitable for use with makeblastdb.
    Raise a FastaValidationError if any of the following are true:

    - a record is empty
    - a record name occurs more than once
    - a sequence occurs more than once
    """
    with dnaio.open(path) as sr:
        records = list(sr)

    names = set()
    sequences = dict()
    for r in records:
        if len(r.sequence) == 0:
            raise FastaValidationError("Record {!r} is empty".format(r.name))
        if r.name in names:
            raise FastaValidationError("Record name {!r} occurs more than once".format(r.name))
        s = r.sequence.upper()
        if s in sequences:
            raise FastaValidationError("Records {!r} and {!r} contain the same sequence".format(
                r.name, sequences[s]))
        sequences[s] = r.name
        names.add(r.name)


def find_overlap(s, t, min_overlap=1):
    """
    Detect if s and t overlap.

    Returns:

    None if no overlap was detected.
    0 if s is a prefix of t or t is a prefix of s.
    Positive int gives index where t starts within s.
    Negative int gives -index where s starts within t.

    >>> find_overlap('ABCDE', 'CDE')
    2
    >>> find_overlap('CDE', 'ABCDEFG')
    -2
    >>> find_overlap('ABC', 'X') is None
    True
    """
    aligner = Aligner(s, max_error_rate=0, min_overlap=min_overlap)
    result = aligner.locate(t)
    if result is None:
        return None
    s_start, _, t_start, _, _, _ = result
    return s_start - t_start


def merge_overlapping(s, t):
    """
    Return merged sequences or None if they do not overlap
    The minimum overlap is 50% of the length of the shorter sequence.
    """
    i = find_overlap(s, t, min_overlap=max(1, min(len(s), len(t)) // 2))
    if i is None:
        return None

    if i >= 0:
        # positive: index of t in s
        if i + len(t) < len(s):
            # t is in s
            return s
        return s[:i] + t
    if -i + len(s) < len(t):
        # s is in t
        return t
    return t[:-i] + s


def get_cpu_time():
    """Return CPU time used by process and children"""
    if sys.platform != 'linux':
        return None
    rs = resource.getrusage(resource.RUSAGE_SELF)
    rc = resource.getrusage(resource.RUSAGE_CHILDREN)
    return rs.ru_utime + rs.ru_stime + rc.ru_utime + rc.ru_stime


def slice_arg(s):
    """
    Parse a string that describes a slice with start and end.

    >>> slice_arg('2:-3')
    slice(2, -3, None)

    >> slice_arg(':-3')
    slice(None, -3, None)

    >> slice_arg('2:')
    slice(2, None, None)
    """

    start, end = s.split(':')
    start = None if start == '' else int(start)
    end = None if end == '' else int(end)
    return slice(start, end)


def is_same_gene(name1: str, name2: str):
    """
    Compare gene names to find out whether they are alleles of each other.
    Both names must have a '*' in them
    """
    return '*' in name1 and '*' in name2 and name1.split('*')[0] == name2.split('*')[0]


class ChimeraFinder:
    def __init__(self, sequences: List[str], min_length: int=10):
        self._sequences = sequences
        self._min_length = min_length
        self._build_index()

    def _build_index(self):
        min_length = self._min_length
        # Create two dictionaries that map all prefixes and suffixes to indexes of all
        # sequences they occur in
        prefixes = defaultdict(list)
        suffixes = defaultdict(list)

        for i, sequence in enumerate(self._sequences):
            if len(sequence) < min_length:
                continue
            for stop in range(min_length, len(sequence) + 1):
                prefix = sequence[:stop]
                prefixes[prefix].append(i)
            for start in range(0, len(sequence) + 1 - min_length):
                suffix = sequence[start:]
                suffixes[suffix].append(i)
        self._prefixes = prefixes
        self._suffixes = suffixes

    def find_exact(self, query: str):
        """
        Find out whether the query string can be explained as a concatenation of
        a prefix of one of the strings plus a suffix of one of the strings in
        the given list of sequences. Both the prefix and the suffix must have a
        length of at least min_length.

        If the answer is yes, a tuple (prefix_length, prefix_indices, suffix_indices)
        is returned where prefix_length is the length of the query prefix,
        prefix_indices is a list of int indices into the sequences of all possible
        prefix sequences that match, and suffix_indices is the same for the suffix.

        The prefix_length returned is the first that yields a result. More are
        possible.

        If the answer is no, None is returned.
        """
        min_length = self._min_length
        for split_index in range(min_length, len(query) + 1 - min_length):
            prefix = query[:split_index]
            suffix = query[split_index:]
            if prefix in self._prefixes and suffix in self._suffixes:
                return (split_index, self._prefixes[prefix], self._suffixes[suffix])

        return None


def available_cpu_count():
    """
    Number of available virtual or physical CPUs on this system
    """
    if sys.platform != 'linux':
        try:
            import multiprocessing
            return multiprocessing.cpu_count()
        except (ImportError, NotImplementedError):
            return 1
    return len(os.sched_getaffinity(0))
