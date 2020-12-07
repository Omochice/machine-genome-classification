from argparse import ArgumentParser, Namespace
from dataclasses import dataclass
import pandas as pd
from pathlib import Path
from os import PathLike
import yaml
from Bio import Entrez, SeqIO


@dataclass
class OrganelleFetchCliant:
    email: str
    n_once: int = 1000
    column_title: str = "Replicons"

    def fetch(self, source: PathLike, dst: PathLike) -> None:
        """csvファイルを元にGenBankからgbkファイルを取得する.

        Args:
            source (PathLike): sourceとなるcsvへのpath
            dst (PathLike): 取得したファイルの保存先

        Returns:
            None:
        """
        # preprocess for download
        dst = Path(dst) / "gbk"
        dst.mkdir(parents=True, exist_ok=True)

        # format replicons for Entrez.efetch
        records = pd.read_csv(source)[self.column_title].map(
            lambda x: self._extract_usable_accession(x))

        # download (1000 records at once in default)
        for i in range(0, len(records), self.n_once):
            queries = records[i: i+self.n_once]
            with Entrez.efetch(db="nucleotide",
                               id=queries,
                               rettype="gb",
                               retmode="text") as handle:
                for rst_record in SeqIO.parse(handle, "genbank"):
                    acc = rst_record.name
                    save_dst = dst / f"{acc}.gbk"
                    if save_dst.exists():
                        continue
                    else:
                        with open(save_dst, "w")as f:
                            SeqIO.write(rst_record, f, "genbank")

    def _extract_usable_accession(self, concated_acc: str) -> str:
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


if __name__ == "__main__":
    args = parser()
    project_dir = Path(__file__).resolve().parents[2]

    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    source = Path(args.input_csv or config["source_csv"]).resolve()
    dst = Path(args.destination or config["data_dst"]).resolve()

    cliant = OrganelleFetchCliant().fetch(source, dst, email=config["email"])
