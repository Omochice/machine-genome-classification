from pathlib import Path
from typing import Tuple
import json
import yaml
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from Bio import Seq, SeqIO
from argparse import Namespace, ArgumentParser

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


def parser() -> Namespace:
    """setting for argparser

    Returns:
        Namespace: args namespace.
    """
    usage = f"Usage: python {__file__} [-w weight] [-wo weight_output] [-i inputfiles] [-o output] "
    argparser = ArgumentParser(usage=usage)
    argparser.add_argument("-w",
                           "--weight",
                           help="If you have already weight data, write the path.")
    argparser.add_argument(
        "-wo",
        "--weight_output",
        help="if not use weight file, Destination of computed weight path.")
    argparser.add_argument("-i",
                           "--inputs",
                           nargs="*",
                           help="Input file pathes",
                           required=True)
    argparser.add_argument(
        "-o",
        "--output_dir",
        help=
        f"Destination of generated_image. if don't specify, use {str(Path(__file__).resolve() / 'setting.yml')}"
    )

    args = argparser.parse_args()
    return args


if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parents[2]
    args = parser()
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    inputs = list(map(lambda x: Path(x).resolve(), args.inputs))
    if args.weight is None:
        weight = calc_entropy.calc_self_entropies(inputs)
        if args.weight_output is not None:
            with open(Path(args.weight_output).resolve()) as f:
                json.dump(weight, f, indent=4)

    else:
        with open(Path(args.weight).resolve()) as f:
            weight = json.load(f)

    if args.output_dir is None:
        dst = config["graph_destination"]
    else:
        dst = Path(args.output_dir).resolve()

    dst.mkdir(exist_ok=True)

    for input_gbk in inputs:
        for record in SeqIO.parse(str(input_gbk), "genbank"):
            acc = record.name
            print(acc)
            fig = generate_image(record.seq, weight)
            plt.savefig(dst / f"{acc}.png")
            plt.close()
