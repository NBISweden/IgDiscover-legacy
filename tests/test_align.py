from igdiscover.align import consensus, multialign, align_affine


def test_consensus():
    assert consensus((
        "TATTACTGTGCGAG---",
        "TATTACTGTGCGAGAGA",
        "TATTACTGTGCGAGAGA",
        "TATTACTGTGCGAGAG-",
        "TATTACTGTGCGAGAG-",
        "TATTACTGTGCGAG---",
        "TATTACTGTGCGAGA--",
    )) == \
        "TATTACTGTGCGAGAGA"


def test_multialign():
    result = multialign({
        "seq1": "TATTACTGTGCGAG",
        "seq2": "TCTTACGTGCTAG",
    }, program="muscle")
    assert result == {
        "seq1": "TATTACTGTGCGAG",
        "seq2": "TCTTAC-GTGCTAG",
    }


# The following examples come from Flouri et al. (2015), http://dx.doi.org/10.1101/031500
# Gap scores were adjusted to take into account that a gap of
# length k costs gap_open + gap_extend * (k-1) in our implementation.


def test_flouri1():
    ref = 'AAAGGG'
    query = 'TTAAAAGGGGTT'
    a = align_affine(ref, query, match=0, mismatch=-1, gap_open=-6, gap_extend=-1)
    assert -15 == a.score
    assert 'AAA------GGG' == a.ref_row
    assert 'TTAAAAGGGGTT' == a.query_row


def test_flouri2():
    ref = 'CGCCTTAC'
    query = 'AAATTTGC'
    a = align_affine(ref, query, match=10, mismatch=-40, gap_open=-41, gap_extend=-1)
    assert -72 == a.score
    assert 'CGCCTTA------C' == a.ref_row
    assert '------AAATTTGC' == a.query_row


def test_flouri3():
    ref = 'TAAATTTGC'
    query = 'TCGCCTTAC'
    a = align_affine(ref, query, match=10, mismatch=-30, gap_open=-41, gap_extend=-1)
    assert -62 == a.score
    assert '------TAAATTTGC' == a.ref_row
    assert 'TCGCCTTA------C' == a.query_row


def test_flouri4():
    ref = 'AAATTTGC'
    query = 'CGCCTTAC'
    a = align_affine(ref, query, match=10, mismatch=-30, gap_open=-41, gap_extend=-1)
    assert -72 == a.score
    assert '------AAATTTGC' == a.ref_row
    assert 'CGCCTTA------C' == a.query_row


def test_flouri5():
    ref = 'AGAT'
    query = 'CTCT'
    a = align_affine(ref, query, match=10, mismatch=-30, gap_open=-29, gap_extend=-1)
    assert '---AGAT' == a.ref_row
    assert 'CTC---T' == a.query_row
