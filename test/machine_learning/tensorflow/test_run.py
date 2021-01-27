from src.machine_learning.tensorflow import run
import numpy as np
import pytest


@pytest.mark.parametrize(
    ("source", "expected"),
    [(("a", "a", "a", "b", "b", "c"),
      ([[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 1, 0], [0, 1, 0], [1, 0, 0]])),
     (("a", "a", "b", "b"), ([[1, 0], [1, 0], [0, 1], [0, 1]]))])
def test_get_sorted_class(source, expected):
    arr = np.array(source, dtype=object)
    expected = np.array(expected)
    assert (run.to_categorical(arr) == expected).all()
