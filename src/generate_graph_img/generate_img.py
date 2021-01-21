import json
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from Bio import Seq, SeqIO
from seqtools.seq_tools import calc_entropy, gbk_utils

matplotlib.use("Agg")


def generate_image(sequense: Seq.Seq,
                   weight: dict,
                   pix_of_a_side: int = 192) -> matplotlib.figure.Figure:
    """配列をグラフ化する 

    Args:
        sequense (Seq.Seq): 塩基配列
        weight (dict): 重み
        pix_of_a_side (int, optional): 生成画像の一辺のピクセル数. Defaults to 192.

    Returns:
        matplotlib.figure.Figure: 生成された画像
    """
    x_coo, y_coo = calc_coordinates(sequense, weight)
    dpi = 100
    figsize = (pix_of_a_side / dpi, pix_of_a_side / dpi)
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.plot(x_coo, y_coo, color="Black", lw=1)
    # 先にそれぞれのcooをソートしてからmin,maxかけるのとどっちが早いんだろう
    format_graph_img(ax, min(x_coo), max(x_coo), min(y_coo), max(y_coo))
    return fig


def format_graph_img(ax: matplotlib.axes._axes.Axes, xmin: np.float64, xmax: np.float64,
                     ymin: np.float64, ymax: np.float64) -> None:
    """縦横比を揃えたりする

    Args:
        fig (matplotlib.figure.Figure): グラフ画像

    Returns:
        matplotlib.figure.Figure: [description]
    """
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.axis("off")


def calc_coordinates(seq: Seq.Seq, weight: dict) -> Tuple[list, list]:
    """入力された配列を数値化する

    Args:
        seq (Seq.Seq): 配列
        weight (dict): 重み

    Returns:
        Tuple[list, list]: x座標, y座標
    """
    VECTORS = {"A": (1, 1), "T": (-1, 1), "G": (-1, -1), "C": (1, -1)}
    x_coo, y_coo = [0], [0]
    for triplet in gbk_utils.window_serach(seq, overhang="before"):
        x_coo.append(x_coo[-1] + VECTORS[triplet[-1]][0] * weight.get(triplet, 1))
        y_coo.append(y_coo[-1] + VECTORS[triplet[-1]][1] * weight.get(triplet, 1))
    return x_coo, y_coo


if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parents[2]
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    df = pd.read_csv(Path(config["destination"]) / "csv" / "creatures.csv")

    query_string = "not is_mongrel and complete and not shotgun and not chromosome"
    accs = df.query(query_string)["accession"]

    acc_to_seq: Dict[str, Seq.Seq] = {}

    # construct list
    gbk_dir = Path(config["destination"]) / "gbk"
    for gbkfile in map(lambda acc: gbk_dir / f"{acc}.gbk", accs):
        for record in SeqIO.parse(gbkfile, "genbank"):
            seq = gbk_utils.get_seq(record,
                                    recursive=True,
                                    search_gbk_root=(gbk_dir / "contigs"))
            acc_to_seq[record.name] = seq

    # generate_weight
    weight = calc_entropy.calc_weights(acc_to_seq.values(), atgc_only=True)
    with open(Path(config["destination"]) / "json" / "weights.json", "w") as f:
        json.dump(weight, f)

    # plot each sequeces
    img_dst = Path(config["destination"]) / "img"
    img_dst.mkdir(parents=True, exist_ok=True)

    for acc, sequence in acc_to_seq.items():
        fig = generate_image(sequence, weight, config["graph_pix"])
        dst = img_dst / f"{acc}.png"
        plt.savefig(dst)
        plt.close()    # clear figure on each species
