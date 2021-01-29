import json
import pandas as pd
from collections import Counter
import sys
from pathlib import Path

import yaml


def filtering_by_n_class(df: pd.DataFrame, n_class: int) -> pd.DataFrame:
    df.query("complete and not shotgun and not chromosome and not is_mongrel",
             inplace=True)
    counter = Counter(df["class"])
    class_geq_n_class = {cl for cl in counter.keys() if counter[cl] >= n_class}
    # df = df.query("class in @class_geq_n_class") # not work
    df = df[df["class"].isin(class_geq_n_class)]
    return df.drop(["complete", "shotgun", "chromosome", "is_mongrel"], axis=1)


def make_roster(df: pd.DataFrame, taxon: str, n_split: int = 4) -> dict:
    # TODO
    # change method huristic to algorizmic
    # counter = Counter(df[taxon])
    all_taxon = set(df[taxon])
    no_1 = {"Actionopteri"}
    no_2 = {"Insecta"}
    no_3 = {"Mammalia", "Aves", "Malcostraca"}
    no_4 = all_taxon - no_1 - no_2 - no_3
    d = {
        "all": list(all_taxon),
        "no1": list(no_1),
        "no2": list(no_2),
        "no3": list(no_3),
        "no4": list(no_4),
    }
    return d


if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parents[2]

    with open(Path(project_dir) / "setting.yml") as f:
        config = yaml.safe_load(f)

    filename = sys.argv[1]
    df = pd.read_csv(project_dir / filename, index_col=0)

    df = filtering_by_n_class(df, config["use_limit"])

    dst = Path(config["destination"]) / "csv" / "formatted.csv"
    df.to_csv(dst)
    with open(Path(config["destination"]) / "json" / "roster.json", "w") as f:
        json.dump(make_roster(df, config["focus_rank"]), f, indent=2)
