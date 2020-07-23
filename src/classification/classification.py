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


def classification(target_dir: Path) -> None:
    """targer_dirに入っているファイルを分類する

    Args:
        target_dir (Path): 対象ファイルが入っているdirectory
    """
    project_dir = Path(__file__).resolve().parents[2]
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)
    taxon_path = Path(config["taxoninfo_destination"])

    d = concate_json(list(taxon_path.glob("result_*.json")))
    targets = [t for t in target_dir.glob("*") if t.is_file]
    valid = target_dir / "valid"
    invalid = target_dir / "invalid"
    valid.mkdir(exist_ok=True)
    invalid.mkdir(exist_ok=True)
    no_classes = []
    for target in targets:
        stem = target.stem
        try:
            taxon = d[stem]["taxon"]
        except KeyError:
            dst = invalid / target.name
        else:
            use_class = get_taxon(taxon, config["focus_rank"], config["priority"])
            if use_class is None:    # どのデータベースにもclassが記載されていない
                no_classes.append(stem)
                dst = invalid / target.name
            else:
                class_root = valid / use_class
                class_root.mkdir(exist_ok=True)
                dst = class_root / target.name
        shutil.move(str(target), str(dst))    # dstはpathlikeみたいだけど一応strをかけた


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


if __name__ == "__main__":
    args = parser()
    target_dir = Path(args.target_dir).resolve()
    classification(target_dir)
