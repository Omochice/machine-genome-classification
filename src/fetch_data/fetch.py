import os
from argparse import ArgumentParser, Namespace
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO


class OrganelleFetchCliant:
    """オルガネラデータベース(https://www.ncbi.nlm.nih.gov/genome/browse#!/organelles/)で
    取得したcsvをもとにgenbankからデータをfetchするクラス
    """
    def __init__(self, limit: int = 1000, column_title: str = "Replicons") -> None:
        """1度のアクセスで何件を取得するか初期設定

        Args:
            limit (int): 1度のアクセスで何件を取得するか. Defaults to 1000.
            columns_title: Accession番号が書いてある列名. Degaults to "Replicons".
        """
        self.limit = limit
        self.use_column_title = "Replicons"
        Entrez.email = os.environ["email"]

    def fetch(self, source: Path, dst: Path) -> None:
        """csvファイルをもとにGenbankからデータを取得する

        Args:
            source (Path): データが記載されたcsvへのパス
            dst (Path): 取得したデータの保存先ディレクトリ
        """

        dst.mkdir(parents=True, exist_ok=True)

        # csvの解析
        records = validate_replicons(parge_csv(source, self.use_column_title))

        # データの取得
        n_records = len(records)
        for i in range(0, n_records, self.limit):
            queries = records[i:i + 1000]
            with Entrez.efetch(db="nucleotide",
                               id=queries,
                               rettype="gb",
                               retmode="text") as handle:
                for r in SeqIO.parse(handle, "genbank"):
                    acc_num = r.name
                    save_path = dst / f"{acc_num}.gbk"
                    if save_path.exists():    # すでに取得しているのならば
                        continue
                    with open(save_path, "w") as f:
                        SeqIO.write(r, f, "genbank")


def parge_csv(source: Path, use_column_title: str) -> pd.Series:
    """csvの使う列を返す

    Args:
        source (Path): csvへのパス
        use_column_title (str): 使う列の列名

    Returns:
        pd.Series: 列データ
    """
    df = pd.read_csv(source)
    return df[use_column_title]


def validate_replicons(replicons: pd.Series) -> list:
    """Replicon列を整形する
    2つのアクセッション番号が書いてれば前者を優先、MT: はアクセッション番号と関係ないので削除

    Args:
        replicons (pd.Series): Replicon列データ

    Returns:
        list: 整形後の列データ
    """
    pretty = [complex_txt.split("/")[0].split(":")[-1] for complex_txt in replicons]
    return pretty


# ここからパーサの設定
def parser() -> Namespace:
    """setting for argparser

    Returns:
        Namespace: args namespace.
    """
    usage = f"Usage: python {__file__} <csv> <destination>"
    argparser = ArgumentParser(usage=usage)
    argparser.add_argument("csv", help="Source csv file path.")
    argparser.add_argument("destination",
                           help="Destination dir path of downloaded file.")

    args = argparser.parse_args()
    return args


if __name__ == "__main__":
    args = parser()
    print(type(args))

    source = Path(args.csv).resolve()
    dst = Path(args.destination).resolve()

    cliant = OrganelleFetchCliant().fetch(source, dst)
