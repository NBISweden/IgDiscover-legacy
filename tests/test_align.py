from igdiscover.align import consensus, multialign


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
