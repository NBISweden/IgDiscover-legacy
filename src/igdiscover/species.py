"""
Species-specific code, such as lists of motifs and regular expressions.

Some refactoring is needed to make this module actually usable for many
species. Right now, it works for - at least - human, rhesus monkey and mouse.
"""
import re
from .utils import nt_to_aa


# Regular expressions for CDR3 detection
#
# The idea comes from D’Angelo et al.: The antibody mining toolbox.
# http://dx.doi.org/10.4161/mabs.27105
# The heavy-chain regex was taken directly from there, but the difference
# is that we express everything in terms of amino acids, not nucleotides.
# This simplifies the expressions and makes them more readable.
#
_CDR3_REGEX = {
    # Heavy chain
    'IGH': re.compile("""
        [FY] [FHVWY] C
        (?P<cdr3>
            [ADEGIKMNRSTV] .{3,31}
        )
        W[GAV]
        """, re.VERBOSE),

    # Light chain, kappa
    'IGK': re.compile("""
        [FSVY] [CFHNVY] [CDFGLSW]
        (?P<cdr3>
            .{4,15}
        )
        [FLV][GRV]
        """, re.VERBOSE),

    # Light chain, lambda
    'IGL': re.compile("""
        # the negative lookahead assertion ensures that the rightmost start is found
        [CDY](?![CDY][CFHSY][CFGW])[CFHSY][CFGW]
        (?P<cdr3>
            .{4,15}
        )
        [FS]G
        """, re.VERBOSE)
}

_CDR3_VH_ALTERNATIVE_REGEX = re.compile("""
        C
        (?P<cdr3> . [RK] .{3,30})
        [WF]G.G
""", re.VERBOSE)


def find_cdr3(sequence, locus):
    """
    Find the CDR3 in the given sequence, assuming it comes from the given locus
    (chain type). If the locus is not one of 'IGH', 'IGK', 'IGL', return None.

    Return a tuple (start, stop) if found, None otherwise.
    """
    try:
        regex = _CDR3_REGEX[locus]
    except KeyError:
        return None
    matches = []
    for offset in 0, 1, 2:
        aa = nt_to_aa(sequence[offset:])
        match = regex.search(aa)
        if not match and locus == "IGH":
            match = _CDR3_VH_ALTERNATIVE_REGEX.search(aa)
        if match:
            start, stop = match.span('cdr3')
            matches.append((start * 3 + offset, stop * 3 + offset))
    return min(matches, default=None)


# The following code is used for detecting CDR3 start sites within V
# reference sequences and CDR3 end sites within J reference sequences.


# Matches the start of the CDR3 within the end of a VH sequence
_CDR3_START_VH_REGEX = re.compile("""
    [FY] [FHVWY] C
    (?P<cdr3_start>
        [ADEGIKMNRSTV*] | $
    )
    """, re.VERBOSE)


_CDR3_START_VH_ALTERNATIVE_REGEX = re.compile("""
    C
    (?P<cdr3_start> . [RK])
    """, re.VERBOSE)


_CDR3_START_REGEXES = {
    "IGK": re.compile('[FSVY][CFHNVY][CDFGLSW]'),
    "IGL": re.compile('[CDY](?![CDY][CFHSY][CFGW])[CFHSY][CFGW]'),
    "TRG": re.compile('[YFH]C'),  # TODO test whether this also works for alpha and beta
    "TRD": re.compile('[YFH]C'),
}


def _cdr3_start_heavy(aa):
    head, tail = aa[:-15], aa[-15:]
    match = _CDR3_START_VH_REGEX.search(tail)
    if not match:
        match = _CDR3_START_VH_ALTERNATIVE_REGEX.search(tail)
    if not match:
        return None
    return len(head) + match.start('cdr3_start')


def cdr3_start(nt, locus):
    """
    Find CDR3 start location within a V gene (Ig or TCR)

    nt -- nucleotide sequence of the gene
    locus -- one of "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"
    """
    aa = nt_to_aa(nt)
    if locus == "IGH":
        start = _cdr3_start_heavy(aa)
        if start is None:
            return None
        return 3 * start
    if locus in ("IGK", "IGL", "TRG", "TRD"):
        head, tail = aa[:-15], aa[-15:]
        match = _CDR3_START_REGEXES[locus].search(tail)
        if match:
            return 3 * (len(head) + match.end())
        else:
            return None
    elif locus in ("TRA", "TRB"):
        head, tail = aa[:-8], aa[-8:]
        pos = tail.find('C')
        if pos == -1:
            return None
        else:
            return 3 * (len(head) + pos + 1)


# Matches after the end of the CDR3 within a J sequence
_CDR3_END_REGEXES = {
    "IGH": re.compile("W[GAV]"),
    "IGK": re.compile("FG"),
    "IGL": re.compile("FG"),
    "TRA": re.compile("FG"),
    "TRB": re.compile("FG"),
    "TRG": re.compile("FG"),
    "TRD": re.compile("FG"),
}


def cdr3_end(nt, locus):
    """
    Find the position of the CDR3 end within a J sequence

    nt -- nucleotide sequence of the J gene
    locus -- one of "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"
    """
    regex = _CDR3_END_REGEXES[locus]
    for frame in 0, 1, 2:
        aa = nt_to_aa(nt[frame:])
        match = regex.search(aa)
        if match:
            return match.start() * 3 + frame
    return None


# When searching for the CDR3, start this many bases to the left of the end of
# the V match.
CDR3_SEARCH_START = 30
