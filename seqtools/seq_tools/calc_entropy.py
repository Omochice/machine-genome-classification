from collections import Counter
from pathlib import Path
from typing import List

import numpy as np
from Bio import SeqIO

from . import gbk_utils


def compute_triplets(gbk_files: List[Path]) -> Counter:
    """使用するgbkファイルの配列から重みを算出する

    Args:
        gbk_files (List[Path]): 使用するgbkファイルのpathの集合体

    Returns:
        Counter: 重み
    """
    triplets = Counter()

    for gbk in gbk_files:
        individudal_triplets = Counter()
        for record in SeqIO.parse(gbk, "genbank"):
            individudal_triplets += Counter(gbk_utils.window_serach(record.seq))
        triplets += individudal_triplets
    return triplets


def calc_self_entropy(triplets: Counter) -> dict:
    """自己情報量に基づいて3塩基のdictのvalueをエントロピーに書き換える

    Args:
        triplets (Counter): 3塩基の出現個数がまとまったCounter

    Returns:
        dict: 3塩基に対するエントロピーがまとまったdict
    """
    bases = {"A", "T", "G", "C"}
    entropies = {}
    for triplet, content in triplets.items():
        first_2_base = triplet[:2]
        denominator = sum([triplets[first_2_base + base] for base in bases])
        entropies[triplet] = -1 * np.log2(content / denominator)
    return entropies


def calc_self_entropies(gbk_files: List[Path]) -> dict:
    triplets = compute_triplets(gbk_files)
    return calc_self_entropy(triplets)