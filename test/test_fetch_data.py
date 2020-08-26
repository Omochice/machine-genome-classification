import unittest
from src.fetch_data.fetch import validate_replicons, parge_csv
from pathlib import Path


class TestFetchData(unittest.TestCase):
    def test_fetch(self):
        """csvのパーサのtest
        """
        def helper(expected: list, actual: Path):
            self.assertEqual(expected,
                             validate_replicons(parge_csv(actual, "Replicons")))

        "0 -> 単純な読み取り, 1 -> 複数のもの, 2 -> MT:の除去"

        source = Path(__file__).resolve().parent / "in" / "0.csv"
        helper(["NC_038132.1", "NC_020356.1", "AC_0123.1"], source)


if __name__ == '__main__':
    unittest.main()