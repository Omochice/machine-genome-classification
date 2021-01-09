from src.data_preparation import fetch_gbk
from pathlib import Path
import pytest
import yaml


@pytest.fixture(scope="session")
def load_setting() -> dict:
    project_dir = Path(__file__).resolve().parents[2]
    with open(project_dir/"setting.yml") as f:
        config = yaml.safe_load(f)
    return config


@pytest.mark.parametrize(("source", "expected"), [
    ("NC_031807.1", "NC_031807.1"),
    ("NC_031807.1/KU146531", "NC_031807.1"),
    ("MT:NC_010195.2/JF275060", "NC_010195.2")
])
def test_extract_usable_accession(source, expected):
    """複数のaccから使える方を選べているか"""
    assert fetch_gbk.extract_usable_accession(source) == expected


def test_fetch(load_setting):
    accessions = ("NC_012920", "NC_031807")
    cliant = fetch_gbk.OrganelleFetchCliant(load_setting["email"])
    assert len([r for r in cliant.fetch(accessions)]) == len(accessions)
