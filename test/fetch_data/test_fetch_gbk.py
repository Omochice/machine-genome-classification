from src.fetch_data import fetch_gbk
from pathlib import Path
import pytest

project_dir = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(("source", "expected"), [
    ("NC_031807.1", "NC_031807.1"),
    ("NC_031807.1/KU146531", "NC_031807.1"),
    ("MT:NC_010195.2/JF275060", "NC_010195.2")
])
def test_extract_usable_accession(source, expected):
    cliant = fetch_gbk.OrganelleFetchCliant(email="example@example.com")
    """複数のaccから使える方を選べているか"""
    assert cliant._extract_usable_accession(source) == expected


def test_pp():
    pass
