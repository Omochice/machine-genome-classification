import json
import shutil
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import List, Optional

import yaml


def concate_json(jsons: List[Path]) -> dict:
    """複数に別れてる分類情報が記載されたjsonを結合する

    Args:
        jsons (List[Path]): 対象のjosnのパスのlist

    Returns:
        dict: 結合したdict
    """
    d = {}
    for j in jsons:
        with open(j) as f:
            d.update(json.load(f))
    return d


def get_taxon(d: dict, focus_rank: str, priority: list) -> Optional[str]:
    """与えられた生物のtaxon dictを調べ、優先順位が上でfocus_rankが含まれているものを返す

    Args:
        d (dict): taxon dict
        focus_rank (str): 調べる対象の分類階層. ex. class, phylo
        priority (list): データベースの優先順位 

    Returns:
        Optional[str]: 見つかった分類階層名かNone
    """
    def serch_rank(k: list) -> Optional[str]:
        """taxon dictを再帰的に調べる

        Args:
            k (list): 調べる残りのデータベース

        Returns:
            Optional[str]: 見つかればstr, そうでなければNone
        """
        if not k:    #どのデータベースもrankを持っていない
            return None
        else:
            binomial_name = d[k[0]].get(focus_rank, None)
            if binomial_name is None:
                return serch_rank(k[1:])
            else:
                return binomial_name

    return serch_rank(priority)


def move_gbk(path: Path, dst_base: Path, invalids: dict, taxon_dict: dict,
             config: dict) -> None:
    """1つのgbkに対して使えるかどうか(雑種ではないetc)を調べ対応するディレクトリに移動する

    Args:
        path (Path): 対象となるgbk
        dst_base (Path): 送り先のベース 
        invalids (dict): 使わない雑種などをまとめたdict
        taxon_dict (dict): 分類情報が乗ったdict
        config (dict): プロジェクトのconfig
    """
    stem = path.stem
    if stem in invalids.keys():
        dst_base = dst_base / "invalid" / invalids[stem]
    else:
        try:
            taxon = get_taxon(taxon_dict[stem]["taxon"], config["focus_rank"],
                              config["priority"])
        except KeyError:
            dst_base = dst_base / "invalid" / "no_taxon"
        else:
            if taxon is not None:
                dst_base = dst_base / "valid" / taxon
            else:
                dst_base = dst_base / "invalid" / "no_class"
    dst_base.mkdir(parents=True, exist_ok=True)
    shutil.move(str(path), str(dst_base / path.name))


def classification(target_dir: Path) -> None:
    """targer_dirに入っているファイルを分類する

    Args:
        target_dir (Path): 対象ファイルが入っているdirectory
    """
    project_dir = Path(__file__).resolve().parents[2]
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)
    taxon_path = Path(config["taxoninfo_destination"])

    with open(config["invalid_creatures"]) as f:
        invalids = json.load(f)
    taxon_dict = concate_json(list(taxon_path.glob("result_*.json")))
    targets = [t for t in target_dir.glob("*.gbk") if t.is_file]
    for target in targets:
        move_gbk(target, target_dir, invalids, taxon_dict, config)


def extract_small_class(target_dir: Path, limit: int, invalid_path: Path) -> None:
    """生物種が少ないクラスは使うのに適していないので除外する

    Args:
        target_dir (Path): 対象のディレクトリ
        limit (int): これ以下の生物しか含まないクラスは除外する
    """
    targets = [c for c in target_dir.glob("*") if c.is_dir]
    for c in targets:
        if len([creature for creature in c.glob("*") if creature.is_file]) <= limit:
            shutil.move(str(c), invalid_path / "small_class")


def parser() -> Namespace:
    """setting for argparser

    Returns:
        Namespace: args namespace.
    """
    usage = f"Usage: python {__file__} <target_dir>"
    argparser = ArgumentParser(usage=usage)
    argparser.add_argument("target_dir", help="Target directory for classification.")

    args = argparser.parse_args()
    return args


import yaml
if __name__ == "__main__":
    args = parser()
    project_dir = Path(__file__).resolve().parents[2]
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)
    target_dir = Path(args.target_dir).resolve()
    classification(target_dir)
    extract_small_class(
        Path(config["gbk_destination"]) / "valid", config["use_limit"],
        Path(config["gbk_destination"]) / "invalid")
