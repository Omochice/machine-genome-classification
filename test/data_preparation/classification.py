import pytest

from pathlib import Path
from src.data_preparation import classification
import json

project_dir = Path(__file__).resolve().parents[2]

with open(project_dir / "test" / "in" / "json" / "classification.json") as f:
    classification_test_dict = json.load(f)


@pytest.mark.parametrize(("data"), classification_test_dict)
def test_get_class(data):
    assert classification.get_class(
        data["priority"], data["data"], "class1") == data["expected"]
