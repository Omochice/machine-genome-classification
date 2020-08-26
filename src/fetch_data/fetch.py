from argparse import ArgumentParser, Namespace
from pathlib import Path

import pandas as pd
import yaml
from Bio import Entrez, SeqIO


class OrganelleFetchCliant:
    """オルガネラデータベース(https://www.ncbi.nlm.nih.gov/genome/browse#!/organelles/)で
    取得したcsvをもとにgenbankからデータをfetchするクラス
    """
    def __init__(self, n_once: int = 1000, column_title: str = "Replicons") -> None:
        """1度のアクセスで何件を取得するか初期設定

        Args:
            n_once (int): 1度のアクセスで何件を取得するか. Defaults to 1000.
            columns_title: Accession番号が書いてある列名. Degaults to "Replicons".
        """
        self.n_once = n_once
        self.use_column_title = column_title
        with open(Path(__file__).resolve().parents[2] / "setting.yml") as f:
            Entrez.email = yaml.safe_load(f)["email"]

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
        for i in range(0, n_records, self.n_once):
            queries = records[i:i + self.n_once]
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
    usage = f"Usage: python {__file__} [-i input_csv] [-d destination]"
    argparser = ArgumentParser(usage=usage)
    argparser.add_argument("-i", "--input_csv", help="Source csv file path.")
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

    source = Path(args.csv or config["source_csv"]).resolve()
    dst = Path(args.destination or config["gbk_destination"]).resolve()

    cliant = OrganelleFetchCliant().fetch(source, dst)
