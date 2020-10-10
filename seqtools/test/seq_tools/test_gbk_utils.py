from seqtools.seq_tools.gbk_utils import get_taxonID
from pathlib import Path

import pytest


@pytest.fixture()
def human_mitchondria_gbk_path():
    return Path(Path(__file__).resolve().parent / "in" / "NC_012920.gbk")


def test_get_taxonID(human_mitchondria_gbk_path):
    assert get_taxonID(human_mitchondria_gbk_path) == "9606"
