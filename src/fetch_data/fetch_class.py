import json
import re
from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Iterator, List, Tuple

import requests
import yaml
from Bio import Entrez, SeqIO

from seqtools.seq_tools import gbk_utils


class TaxonFetchCliant:
    def __init__(self, n_once: int = 100) -> None:
        self.n_once = n_once
        with open(Path(__file__).resolve().parents[2] / "setting.yml") as f:
            Entrez.email = yaml.safe_load(f)["email"]

    def fetch_taxon_infos(self, pathes: List[Path], dst: Path,
                          invalids_path: str) -> None:
        """インターネットから生物の分類情報を取得する

        Args:
            pathes (List[Path]): 取得対象の生物の情報が記載されたgenbank形式のpathの集合体
            dst (Path): 取得結果を書き出すディレクトリ

        """
        dst.mkdir(exist_ok=True, parents=True)
        # 集合体から雑種とかを取り除く
        valid_creatures, reasons = self.extract_invalid_creature(pathes)
        with open(invalids_path, "w") as f:
            json.dump(reasons, f, indent=4)

        # ncbi, globalnameresolverをそれぞれ検索する
        for i in range(0, len(valid_creatures), self.n_once):
            results = {}
            use_pathes = valid_creatures[i:i + self.n_once]
            taxon_ids = list(map(gbk_utils.get_taxonID, use_pathes))
            names = list(map(self._get_creature_name, use_pathes))
            for path_name, binomial_name, res_ncbi, res_gnr in zip(
                    use_pathes, names, fetch_taxon_from_NCBI(taxon_ids),
                    fetch_taxon_from_GNR(names)):
                taxon = {**{"NCBI Taxonomy": res_ncbi}, **res_gnr}
                stem = path_name.stem
                results[stem] = {
                    "binomial_name": binomial_name,
                    "accession": stem,
                    "taxon": taxon
                }
            # n_once件取得したらファイルに書き出す
            with open(dst / f"result_{i}.json", "w") as f:
                json.dump(results, f, indent=4)

    def extract_invalid_creature(self,
                                 creatures: List[Path]) -> Tuple[List[Path], dict]:
        """入力されたgbkのパスの集合体から使用できるものを抽出する

        Args:
            creatures (List[Path]): gbkのpathの集合体

        Returns:
            Tuple[list[Path], dict]: 使用できるpathのlist, 使用できないデータの理由とpathのdict
        """
        valids = []
        invalids = {}
        for creature in creatures:
            stem = creature.stem
            for record in SeqIO.parse(creature, "genbank"):
                title = record.description
                creature_name = record.annotations["organism"]
                if "complete" not in title:    # complete genome ではない
                    invalids[stem] = "not_complete_genome"
                elif " x " in creature_name:    # 〜 x 〜　の雑種
                    invalids[stem] = "mongrel"
                elif len(re.findall(r"[atgc]", str(record.seq.lower()))) == 0:    #配列がない
                    invalids[stem] = "have_no_seq"
                else:
                    valids.append(creature)
        return valids, invalids

    def _get_creature_name(self, p: Path) -> str:
        for record in SeqIO.parse(p, "genbank"):
            return record.annotations["organism"]
        raise NotFoundOrganismError(f"Not found organism name in {str(p)}")


class NotFoundOrganismError(Exception):
    pass


def fetch_taxon_from_NCBI(taxon_ids: list,
                          include_with_norank: bool = False) -> Iterator[dict]:
    """NCBIから生物種名と分類情報の情報を取得する

    Args:
        taxon_ids (list): TaxonIDの配列
        include_with_norank (bool): NCBI Taxonomy上でno rankとされるものも返すか. Defaults to False.

    Yields:
        Iterator[dict]: 取得した情報
    """
    result = {}
    with Entrez.efetch(db="Taxonomy", id=taxon_ids) as efetch_handle:
        efetch_record = Entrez.read(efetch_handle)
    for record in efetch_record:
        if include_with_norank:
            for index, lineage in enumerate(record["LineageEx"]):
                result[index] = {
                    'ScientificName': lineage["ScientificName"],
                    'Rank': lineage["Rank"]
                }
        else:
            result = {
                lineage["Rank"]: lineage["ScientificName"]
                for lineage in record["LineageEx"] if lineage["Rank"] != "no rank"
            }
        yield result


def fetch_taxon_from_GNR(names: list, priority: List[int] = None) -> Iterator[dict]:
    """Globak Names ResolverのAPIを叩いて分類情報を取得する

    Args:
        names (list): 生物の学名の集合体
        priority (List[int]): 使用するsource databases. http://resolver.globalnames.org/data_sourcesを参照. Defaults to None.

    Yields:
        Iterator[dict]: 取得した情報とソースのデータベースをまとめたdict
    """
    url = "http://resolver.globalnames.org/name_resolvers.json"
    if priority is None:
        priority = [4, 179, 11, 1, 8]

    invalid_ranknames = {"", "no rank", "no rank - terminal"}    # familyとかclassとか
    invalid_name = {"Not assigned"}    # mammaliaとか入る方
    params = {
        "names": "|".join(names),
        "best_match_only": True,
        "preferred_data_sources": "|".join(map(str, priority))
    }
    r = requests.get(url, params=params)
    try:
        data = r.json()
    except Exception:
        print(f"Error: {r}")
        import sys
        sys.exit(1)
    for entry in data["data"]:
        result = {}
        for record in entry["preferred_results"]:
            individual_db = {}
            for taxon_name, taxon_rank in zip(
                    record["classification_path"].split("|"),
                    record["classification_path_ranks"].split("|")):
                if taxon_name in invalid_name or taxon_rank in invalid_ranknames:
                    continue
                individual_db[taxon_rank] = taxon_name
            result[record["data_source_title"]] = individual_db
        yield result


def parser() -> Namespace:
    """setting for argparser

    Returns:
        Namespace: args namespace.
    """
    usage = f"Usage: python {__file__} [-i inputs...] [-d destination]"
    argparser = ArgumentParser(usage=usage)
    argparser.add_argument("-i", "--inputs", nargs="*", help="Target gbk files.")
    argparser.add_argument("-d",
                           "--destination",
                           help="Destination of taxonomy json data.")

    args = argparser.parse_args()
    return args


if __name__ == "__main__":
    args = parser()
    project_dir = Path(__file__).parents[2]

    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    inputs = args.inputs or list(Path(config["gbk_destination"]).glob("**/*.gbk"))
    inputs = list(map(lambda x: Path(x).resolve(), inputs))
    dst = Path(args.destination or config["taxoninfo_destination"]).resolve()

    cliant = TaxonFetchCliant()
    cliant.fetch_taxon_infos(inputs, dst, config["invalid_creatures"])
