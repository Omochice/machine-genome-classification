from Bio import SeqIO
from collections import Counter
from pathlib import Path
import yaml
import pandas as pd
from glob import glob


def base_rate(seq) -> dict:
    counter = Counter(seq)
    rst = {b: counter[b] / len(seq) for b in "ATGC"}
    return rst


if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parents[1]
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)
    graph = config["graph_destination"]
    print(graph)
    gbk_dir = Path(config["gbk_destination"]) / "valid"
    arr = []
    for cl in gbk_dir.glob("*"):
        cl_counter = {}
        for gbk in cl.glob("*.gbk"):
            seq = ""
            for record in SeqIO.parse(gbk, "genbank"):
                seq += record.seq
            rate = base_rate(seq)
            image_path = Path(graph) / cl.name / (gbk.stem + ".png")
            if not image_path.exists():
                print(image_path, "is not exist.")
                continue  # 保険、全種画像もあるはず
            row = [
                gbk.stem, cl.name,
                str(image_path),
                len(seq), rate["A"], rate["T"], rate["G"], rate["C"]
            ]
            arr.append(row)
    df = pd.DataFrame(arr,
                      columns=[
                          "accesion", "classname", "imagepath", "seq_len", "rate_A",
                          "rate_T", "rate_G", "rate_C"
                      ])
    df.to_csv("features.csv", index=False)
