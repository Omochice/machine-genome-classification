import re
from pathlib import Path
from typing import Iterator

from Bio import Seq, SeqIO


def get_taxonID(path: Path) -> str:
    """対象のgbkファイルのrecordからtaxonIDを抽出する

    Args:
        record (SeqRecord.SeqRecord): 対象のrecord

    Returns:
        str: db_xrefに記載されたTaxonID
    """
    for record in SeqIO.parse(path, "genbank"):
        for feature in record.features:
            if feature.type == "source":
                db_xref = feature.qualifiers["db_xref"][0]
                taxonID = db_xref.split(":")[1]
                return taxonID
    raise NotFoundTaxonIDError(f"Not Found taxonID in {path}")


class NotFoundTaxonIDError(Exception):
    pass


def window_serach(sequence: Seq.Seq,
                  each: int = 3,
                  overhang: str = None,
                  atgc_only: bool = True) -> Iterator[str]:
    """入力された配列をwindow searchする

    Args:
        sequence (Seq.Seq): 入力配列
        each (int, optional): 何文字ごとに切り出すか. Defaults to 3.
        overhang (str, optional): <each>文字に満たない始端,終端を必要とするか. {None, "before", "after", "both"}のうちのどれか. Defaults to None.
        atgc_only (bool, optional): ATGCのみに整形するかどうか. Defaults to True.

    Yields:
        Iterator[str]: 切り出した<each>文字を返すIterator
    """

    if atgc_only:
        formatted = "".join(re.findall(r"[ATGC]", str(sequence.upper())))
    else:
        formatted = str(sequence.upper())

    if overhang in {"before", "both"}:
        for i in range(1, each):
            yield formatted[:i]

    for i in range(0, len(formatted) - each + 1):
        yield formatted[i:i + each]

    if overhang in {"after", "both"}:
        for i in range(-1 * each + 1, 0, 1):
            yield formatted[i:]


if __name__ == "__main__":
    file = Path(
        "/home/mochi/workspace/master_thesis/test/testdata/output/NC_005958.gbk")

    for record in SeqIO.parse(file, "genbank"):
        # print(get_taxonID(record))

        # print(type(record.seq))
        for window in window_serach(record.seq):
            print(window)
            # pass
