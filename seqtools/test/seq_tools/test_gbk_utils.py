from seqtools.seq_tools import gbk_utils
from pathlib import Path
from Bio.Seq import Seq

import pytest


def load(gbkname) -> Path:
    """test_get_creature_name.

    Args:
        gbk:
        expected:
    """
    return Path(Path(__file__).resolve().parent / "in" / gbkname)


@pytest.mark.parametrize(("gbk", "expected"), [
    ("NC_039811.gbk", "45628"),
    ("NC_050359.gbk", "1130134"),
    ("LT671463.gbk", "269814"),
    ("NC_012920.gbk", "9606")
])
def test_get_taxonID(gbk, expected):
    assert gbk_utils.get_taxonID(load(gbk)) == expected


@pytest.mark.parametrize(("gbk", "expected"), [
    ("NC_039811.gbk", "Gecarcoidea natalis isolate GN2 mitochondrion, partial genome"),
    ("NC_050359.gbk", "Haementeria acuecueyetzin mitochondrion, partial genome"),
    ("LT671463.gbk", "Kudoa iwatai mitochondrion, chromosome 2, complete sequence"),
    ("NC_012920.gbk", "Homo sapiens mitochondrion, complete genome")
])
def test_get_definition(gbk, expected):
    assert gbk_utils.get_definition(load(gbk)) == expected


@pytest.mark.parametrize(("gbk", "expected"), [
    ("NC_039811.gbk", "Gecarcoidea natalis"),
    ("NC_050359.gbk", "Haementeria acuecueyetzin"),
    ("LT671463.gbk", "Kudoa iwatai"),
    ("NC_012920.gbk", "Homo sapiens")
])
def test_get_creature_name(gbk, expected):
    assert gbk_utils.get_creature_name(load(gbk)) == expected


@pytest.mark.parametrize(("gbk", "expected"), [
    ("NC_012920.gbk", True),
    ("no_seq_creature.gbk", False)
])
def test_has_seq(gbk, expected):
    assert gbk_utils.has_seq(load(gbk)) == expected


@pytest.mark.parametrize(("name", "expected"), [
    ("Homo sapiens", False),
    ("foo x bar", True),
    ("includexinname", False),
    ("include xinname", False)
])
def test_is_mongrel(name, expected):
    assert gbk_utils.is_mongrel(name) == expected


@pytest.mark.parametrize(("gbk", "expected"), [
    ("NC_039811.gbk", False),
    ("NC_050359.gbk", False),
    ("LT671463.gbk", False),
    ("NC_012920.gbk", True),
])
def test_is_complete_genome(gbk, expected):
    assert gbk_utils.is_complete_genome(
        gbk_utils.get_definition(load(gbk))) == expected


@pytest.mark.parametrize(("source", "expected"), [
    ("ATGCATGC", "ATGCATGC"),
    ("AAAANNNN", "AAAA"),
    ("aTaT", "aTaT"),
    ("aaaannnn", "aaaa")
])
def test_to_only_atgc(source, expected):
    assert gbk_utils.to_only_actg(source) == Seq(expected)
