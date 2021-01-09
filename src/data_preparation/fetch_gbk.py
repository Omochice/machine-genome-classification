from argparse import ArgumentParser, Namespace
from dataclasses import InitVar, dataclass
from os import PathLike
from pathlib import Path
from typing import Iterable, Iterator

import yaml
from collections import deque
from Bio import Entrez, SeqIO, SeqRecord
from seqtools.seq_tools import utils


@dataclass
class OrganelleFetchCliant:
    email: InitVar[str]
    n_once: int = 1000

    def __post_init__(self, email: str) -> None:
        Entrez.email = email

    def fetch(self, accessions: Iterable[str]) -> Iterator[SeqRecord.SeqRecord]:
        """Fetch genbank records.

        Args:
            accessions (Iterable[str]): An iterable object that have accession numbers.

        Yields:
            Iterator[SeqRecord.SeqRecord]: An SeqRecord iterator.
        """
        for use_accs in utils.split_per_n(accessions, self.n_once):
            use_accs = tuple(use_accs)
            queue = deque()
            queue.append(use_accs)
            # ダウンロードができなければクエリを半分にする
            while queue:
                query = queue.popleft()
                try:
                    for r in self._efetch(query):
                        yield r
                except Exception:
                    if n := len(query) > 1:
                        queue.append(query[:n // 2])
                        queue.append(query[n // 2:])
                    else:
                        for acc in query:
                            print(f"Cannot fetch {acc}.")
                            yield SeqRecord.SeqRecord(seq="", name=acc)

    def _efetch(self, accs: Iterable[str]) -> Iterator[SeqRecord.SeqRecord]:
        with Entrez.efetch(db="nucleotide", id=accs, rettype="gb",
                           retmode="text") as handle:
            for rst_record in SeqIO.parse(handle, "genbank"):
                yield rst_record

    # def fetch(self, source: PathLike, dst: PathLike) -> Iterator[SeqRecord.SeqRecord]:
    #     """csvファイルを元にGenBankからgbkファイルを取得する.

    #     Args:
    #         source (PathLike): sourceとなるcsvへのpath
    #         dst (PathLike): 取得したファイルの保存先

    #     Returns:
    #         None:
    #     """
    #     # preprocess for download
    #     dst = Path(dst)
    #     dst.mkdir(parents=True, exist_ok=True)

    #     # format replicons for Entrez.efetch
    #     records = pd.read_csv(source)[self.column_title].map(
    #         lambda x: self._extract_usable_accession(x))

    #     # download (1000 records at once in default)
    #     for use_gbk in utils.split_per_n(records, n=self.n_once):
    #         queries = [r for r in use_gbk if not self._is_exist_gbk(r, dst)]
    #         if not queries:
    #             continue
    #         else:
    #             yield self.fetch_from_acc(queries)

    # def fetch_from_acc(self, accs: Iterable[str]) -> Iterator[SeqRecord.SeqRecord]:
    #     with Entrez.efetch(db="nucleotide", id=accs, rettype="gb",
    #                        retmode="text") as handle:
    #         for rst_record in SeqIO.parse(handle, "genbank"):
    #             yield rst_record


def extract_usable_accession(concated_acc: str) -> str:
    """repliconsから必要なものだけ抜き出す
    ex) "NC_012920" -> "NC_012920"
        "MT:NC_012920" -> "NC_012920"
        "NC_012920/NN_00001" -> "NC_012920"
        "MT:NC_012920/NN_00001" -> "NC_012920"
    Args:
        concated_acc (str): 複数のaccession番号が結合したもの
    Returns:
        str:
    """
    return concated_acc.split("/")[0].split(":")[-1]


def is_exist_gbk(acc: str, search_root: PathLike) -> bool:
    """Check gbk file exists already.

    Args:
        acc (str): accession number.
        search_root (PathLike): the root dir tu search.

    Returns:
        bool:
    """
    search_root = Path(search_root).resolve()
    acc = acc.split(".")[0]
    return bool(list(search_root.glob(f"**/{acc}.gbk")))


def parser() -> Namespace:
    """setting for argparser

    Returns:
        Namespace: args namespace.
    """
    usage = f"Usage: python {__file__} [-s source_csv] [-d destination]"
    argparser = ArgumentParser(usage=usage)
    argparser.add_argument("-s", "--souse_csv", help="Source csv file path.")
    argparser.add_argument("-d",
                           "--destination",
                           help="Destination dir path of downloaded file.")

    args = argparser.parse_args()
    return args


def main(source: Path, dst: Path, email):
    cliant = OrganelleFetchCliant(email)
    cliant.fetch(source, dst)


if __name__ == "__main__":
    args = parser()
    project_dir = Path(__file__).resolve().parents[2]

    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    source = Path(args.source_csv or config["source_csv"]).resolve()
    dst = Path(args.destination or Path(config["destination"]) / "gbk").resolve()

    main(source, dst, config["email"])
