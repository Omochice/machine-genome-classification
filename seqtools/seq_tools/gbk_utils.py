import re
from pathlib import Path
from typing import Iterator, AnyStr, Tuple, Dict
from os import PathLike
from collections import Counter

from Bio import Seq, SeqIO, SeqRecord


def get_taxonID(path: PathLike) -> str:
    """対象のgbkファイルのrecordからtaxonIDを抽出する

    Args:
        record (PathLike): 対象のファイルのpath

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


def get_definition(path: PathLike) -> str:
    """get definition.

    Args:
        path (Path): 対象のファイルのpath

    Returns:
        str: 生物の学名
    """
    for record in SeqIO.parse(path, "genbank"):
        return record.description


def get_creature_name(path: PathLike) -> str:
    """get creature name.

    Args:
        path (PathLike): 対象のファイルのpath

    Returns:
        str:
    """
    for record in SeqIO.parse(path, "genbank"):
        return record.annotations["organism"]


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


def get_atgc_rate(seq: Seq) -> Dict[str, int]:
    pretty_seq = re.sub("[^ATGC]", "", seq.upper())
    n_seq = len(pretty_seq)
    counter = Counter(pretty_seq)
    rate = {k: v / n_seq for k, v in counter.items()}
    return rate


def to_only_actg(seq: AnyStr) -> Seq.Seq:
    """Change source str like object to {ATGCatgc} only format.

    Args:
        seq (AnyStr): seq

    Returns:
        Seq.Seq:
    """

    return Seq.Seq(re.sub("[^ATGCatgc]", "", str(seq)))


def has_seq(gbk: PathLike) -> bool:
    """与えられたgbkファイルが有効な配列長を持つかどうかを返す.

    Args:
        gbk (PathLike): gbkファイルへのpath

    Returns:
        bool:
    """
    return any([len(to_only_actg(rec.seq)) for rec in SeqIO.parse(gbk, "genbank")])


def is_mongrel(name: str) -> bool:
    """'~ x ~'で書かれる雑種かどうかを返す.

    Args:
        name (str): 生物種

    Returns:
        bool:
    """
    return " x " in name


def is_complete_genome(definition: str) -> bool:
    """完全なミトコンドリアゲノムかどうかを返す.

    Args:
        definition (str): definition

    Returns:
        bool:
    """
    return "mitochondrion, complete genome" in definition


contig_pattern = re.compile(r"join(\(complement)?\((\w*(\.\d)?):(\d+)\.\.(\d+)\)\)?")


def parse_contig(contig: str) -> dict:
    """return contig doscription

    Args:
        contig [str]: A contig string.

    Returns:
        dict: A dict of contig informations.
    """
    match = contig_pattern.match(contig)
    is_complement = bool(match.groups(1))
    accession = match.groups(2)[1].split(".")[0]
    start = int(match.group(4)) - 1
    end = int(match.group(5))

    return {
        "accession": accession,
        "is_complement": is_complement,
        "start": start,
        "end": end
    }


if __name__ == "__main__":
    file = Path(
        "/home/mochi/workspace/master_thesis/test/testdata/output/NC_005958.gbk")

    for record in SeqIO.parse(file, "genbank"):
        # print(get_taxonID(record))

        # print(type(record.seq))
        for window in window_serach(record.seq):
            print(window)
            # pass
