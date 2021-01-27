import pandas as pd
from collections import Counter
import sys
from pathlib import Path


def filtering_by_n_class(df: pd.DataFrame, n_class: int) -> pd.DataFrame:
    df.query("complete and not shotgun and not chromosome and not is_mongrel",
             inplace=True)
    counter = Counter(df["class"])
    class_geq_n_class = {cl for cl in counter.keys() if counter[cl] >= n_class}
    # df = df.query("class in @class_geq_n_class") # not work
    df = df[df["class"].isin(class_geq_n_class)]
    return df.drop(["complete", "shotgun", "chromosome", "is_mongrel"], axis=1)


if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parents[2]
    filename = sys.argv[1]
    df = pd.read_csv(project_dir / filename, index_col=0)

    df = filtering_by_n_class(df, 5)

    dst = Path(filename).parent / "formatted.csv"
    df.to_csv(dst)
