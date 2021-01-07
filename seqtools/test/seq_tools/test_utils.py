from seqtools.seq_tools import utils
import pytest


@pytest.mark.parametrize(("source", "expected"), [\
    (list(range(0, 9)), [(0, 1, 2), (3, 4, 5), (6, 7, 8)]),  #
    (list(range(0, 8)), [(0, 1, 2), (3, 4, 5), (6, 7)]),  #
                                                  ])
def test_split_per_n(source, expected):
    for s, e in zip(utils.split_per_n(source, n=3), expected):
        assert tuple(s) == tuple(e)