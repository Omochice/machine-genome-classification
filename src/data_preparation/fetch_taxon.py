from Bio import Entrez
from collections import deque
import yaml
from argparse import Namespace, ArgumentParser
import requests
from seqtools.seq_tools import gbk_utils, utils
from dataclasses import dataclass, InitVar
from os import PathLike
from typing import Iterator, List, Iterable, Dict
from pathlib import Path


@dataclass
class TaxonFetchCliant:
    email: InitVar[str]
    n_once: int = 100

    def __post_init__(self, email) -> None:
        Entrez.email = email

    def fetch(self, pathes: Iterable[PathLike]) -> Iterator[dict]:
        """Fetch taxon information via NCBI and Global Names Resolver.

        Args:
            pathes (Iterable[PathLike]): Iterable object to gbk file's path.

        Yields:
            Iterator[dict]: Dict accession to taxon information.
        """
        # fetch taxon data NCBI and GlobalNameResolver
        for use_gbks in utils.split_per_n(pathes, n=self.n_once):
            use_gbks = tuple(use_gbks)    # it is used  three times.
            taxon_ids = map(gbk_utils.get_taxonID, use_gbks)    # it is used once.
            names = tuple(map(gbk_utils.get_creature_name,
                              use_gbks))    # it is used twice.
            results = {}

            for gbk, binomial_name, res_ncbi, res_gnr in zip(
                    use_gbks, names, self._fetch_taxon_from_NCBI(taxon_ids),
                    self._fetch_taxon_from_GNR(names)):
                taxon = {**{"NCBI Taxonomy": res_ncbi.copy()}, **res_gnr}
                # copyじゃないと次のループで書き換えられる
                # res_ncbiがポインタで渡されてるっぽい
                acc = gbk.stem
                results[acc] = {
                    "binmomial_name": binomial_name,
                    "accession": acc,
                    "taxon": taxon
                }
            yield results

    def _fetch_taxon_from_NCBI(self,
                               taxon_ids: Iterable,
                               include_with_norank: bool = False) -> Iterator[dict]:
        """fetch information from NCBI taxonomy DB.

        Args:
            taxon_ids (Iterable): Iterable object that have accession numbers.
            include_with_norank (bool): Leave ambigous taxon informations?.

        Yields:
            Iterator[dict]: [description]
        """
        result = {}
        with Entrez.efetch(db="Taxonomy", id=taxon_ids) as handle:
            records = Entrez.read(handle)
        for record in records:
            if include_with_norank:
                for index, lineage in enumerate(record["LineageEx"]):
                    result[index] = {
                        "ScientificName": lineage["ScientificName"],
                        "Rank": lineage["Rank"].lower()
                    }    # too long to use list comprehensions
            else:
                for lineage in record["LineageEx"]:
                    if lineage["Rank"] != "no rank":
                        result[lineage["Rank"].lower()] = lineage["ScientificName"]
            yield result

    def _fetch_taxon_from_GNR(self,
                              names: Iterable,
                              priority: List[int] = None) -> Iterator[dict]:
        """_fetch_taxon_fron_GNR.

        Args:
            names (Iterable): Iterable object that have biological names.
            priority (List[int]): Priority of database to use.
                                  See http://resolver.globalnames.org/data_sources

        Yields:
            Iterator[dict]:
        """
        names = tuple(names)
        url = "http://resolver.globalnames.org/name_resolvers.json"
        if priority is None:
            priority = [4, 3, 179, 11, 1, 8]    # decide by heuristic method

        # 正常なレスポンスが帰ってこなかったら半分でもう一度APIを叩く
        queue = deque()
        queue.append(names)
        while (queue):
            names_query = queue.popleft()
            params = self._params_for_gnr(names_query, priority)
            r = requests.get(url, params=params)
            try:
                data = r.json()
            except Exception:  # gnrが正常なデータを返さなかったら
                if len(names_query) == 1:
                    yield {}
                else:
                    n_queries = len(names_query)
                    queue.append(names_query[:n_queries//2])
                    queue.append(names_query[n_queries//2:])
            else:  # 正常なデータのとき(exceptに入らないとき)
                for entry in data["data"]:
                    yield self._parse_gnr_response(entry)

    def _parse_gnr_response(self, entry: list) -> Dict[str, str]:
        """Parse gnr response.

        Args:
            entry (list): Individual organism's response data entry.

        Returns:
            Dict[str, str]: Dict to rank name to rank. ex) class: Mammalia.
        """
        invalid_ranks_and_names = {"", "no rank", "no rank - terminal", "Not assingned"}
        results = {}
        for record in entry.get("preferred_results", []):
            ind_db = {}
            for rank, name in zip(record["classification_path_ranks"].split("|"),
                                  record["classification_path"].split("|")):
                if rank in invalid_ranks_and_names or name in invalid_ranks_and_names:
                    continue
                else:
                    ind_db[rank] = name
            results[record["data_source_title"]] = ind_db
        return results

    def _params_for_gnr(self, names: Iterable[str], priority: List[int]) -> dict:
        """Construct params for GNR.

        Args:
            names (Iterable[str]): Biological names used for requests.
            priority (List[int]): Priority of database.

        Returns:
            dict: Param for requests.
        """
        params = {
            "names": "|".join(names),
            "best_match_only": True,
            "preferred_data_sources": "|".join(map(str, priority))
        }
        return params


def parser() -> Namespace:
    """setting for argparser

    Returns:
        Namespace: args namespace.
    """
    usage = f"Usage: python {__file__} [-s sources...] [-d destination]"
    argparser = ArgumentParser(usage=usage)
    argparser.add_argument("-s",
                           "--sources",
                           nargs="*",
                           help="source gbk files.(optional)")
    argparser.add_argument("-d",
                           "--destination",
                           help="Destination of output data.(optional)")

    args = argparser.parse_args()
    return args


def main(sources: Iterable[Path], dst: Path, email: str):
    cliant = TaxonFetchCliant(email)
    cliant.fetch(sources, dst)


if __name__ == "__main__":
    args = parser()
    project_dir = Path(__file__).parents[2]

    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    sources = args.sources or (Path(config["destination"]) / "gbk" /
                               "valid").glob("**/*.gbk")
    sources = list(map(lambda x: Path(x).resolve(), sources))
    dst = Path(args.destination or config["destination"]).resolve() / "json"

    main(sources, dst, config["email"])
