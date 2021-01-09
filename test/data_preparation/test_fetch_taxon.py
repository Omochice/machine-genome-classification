import pytest

from src.data_preparation import fetch_taxon

from pathlib import Path
import yaml


@pytest.fixture(scope="session")
def load_setting():
    project_dir = Path(__file__).resolve().parents[2]
    with open(project_dir/"setting.yml") as f:
        config = yaml.safe_load(f)
    return config


@pytest.fixture(scope="session")
def test_in() -> Path:
    return Path(__file__).resolve().parents[1] / "in"


def test_fetch(load_setting, test_in):
    cliant = fetch_taxon.TaxonFetchCliant(load_setting["email"])
    gbks = [file for file in test_in.glob("gbk/*.gbk")]
    d = {}
    for r in cliant.fetch(gbks):
        d.update(r)
    assert len(d) == len(gbks)
