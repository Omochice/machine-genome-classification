from itertools import groupby
from typing import Iterable, Iterator


def split_per_n(it: Iterable, n: int) -> Iterator:
    for _, item in groupby(enumerate(it), lambda x: x[0] // n):
        yield (x[1] for x in item)