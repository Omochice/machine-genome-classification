import json
import shutil
from pathlib import Path
from typing import List, Optional, Iterable
from os import PathLike

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


def get_class(priority: list, taxon_info: dict, use_rank: str) -> Optional[str]:
    """dictからpriorityに基づいて最初の一致を取り出す.

    Args:
        priority (list): priority of db
        taxon_info (dict): dict of db to taxon dict
        use_rank (str): use rank

    Returns:
        Optional[str]:
    """
    for db_name in priority:
        if use_rank in taxon_info[db_name].keys():
            return taxon_info[db_name][use_rank]
    return None


def move_gbk(source: PathLike, class_name: Optional[str]) -> None:
    """move gbk file.

    Args:
        source (PathLike): source gbk file path
        class_name (Optional(str)): class name

    Returns:
        None:
    """
    gbk_dir = Path(source).parents[1]
    valid_dir = gbk_dir / "valid"
    invalid_dir = gbk_dir / "invalid"

    valid_dir.mkdir(parents=True, exist_ok=True)
    invalid_dir.mkdir(parents=True, exist_ok=True)

    if class_name is None:
        # classがなかったら
        dst = invalid_dir / "noclass" / source.name
        dst.parent.mkdir(exist_ok=True)
        shutil.move(source, dst)
    else:
        dst = valid_dir / class_name / source.name
        dst.parent.mkdir(exist_ok=True)
        shutil.move(source, dst)


def get_taxon(accession: str, priority: Iterable[str], taxon_dict: dict) -> str:
    """the wrapper to get_class.

    Args:
        accession (str): accession
        priority (Iterable[str]): priority
        taxon_dict (dict): taxon_dict

    Returns:
        str:
    """

    return get_class(priority, taxon_dict[accession]["taxon"], use_rank="class")


def classification(config: dict) -> None:
    """分類情報を元にディレクトリに分ける.

    Args:
        config (dict): config dict

    Returns:
        None:
    """
    taxon_path = Path(config["destination"]) / "taxon"
    valid_gbk_dir = Path(config["destination"]) / "gbk" / "valid"

    no_class_features = []

    taxon_dict = concate_json(list(taxon_path.glob("result_*.json")))

    for accession, value in taxon_dict.items():
        focused_gbkfile = valid_gbk_dir / f"{accession}.gbk"

        focused_taxon_name = get_class(config["priority"], value["taxon"])
        if focused_taxon_name is None:
            no_class_features.append(accession)
        move_gbk(focused_gbkfile, focused_taxon_name)
    return no_class_features


def main(config: dict):
    classification(config)


if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parents[2]
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)
    main(config)
