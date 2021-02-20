import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image

import itertools
import matplotlib
from os import PathLike
import os
from sklearn.metrics import confusion_matrix

import numpy as np
import pandas as pd
import seaborn as sns

from typing import Iterable, List
sns.set()
import japanize_matplotlib

matplotlib.use("Agg")


def visualize_history(history: dict, title: str = "", dst: PathLike = "") -> None:
    fig, (axL, axR) = plt.subplots(ncols=2, figsize=(10, 4))

    axL.plot(history.get("acc", history.get("accuracy", None)),
             "o-",
             label="Train accuracy")
    axL.plot(history.get("val_acc", history.get("val_accuracy")),
             "o-",
             label="Validation accuracy")
    if "f1score" in history:
        axL.plot(history["f1score"], "*-", label="F1 score")
        axL.set_title("Accuracy and F1 score")
        axL.set_ylabel("value of score")
    else:
        axL.set_title("Accuracy")
        axL.set_ylabel("Accuracy")
    axL.set_xlabel("Epoch")
    axL.set_ylim(0, 1)
    axL.grid(True)
    # axL.legend(bbox_to_anchor=(0, 0), loc="lower left", borderaxespad=0)
    axL.legend(bbox_to_ancher=(1, -0.1), loc="upper right", borderaxespad=0)
    # axL.legend(loc="best")

    axR.plot(history["loss"], "o-", label="Train loss")
    axR.plot(history["val_loss"], "o-", label="Validation loss")
    axR.set_title("Loss")
    axR.set_xlabel("Epoch")
    axR.set_ylabel("Loss")
    axR.grid(True)
    # axR.legend(bbox_to_anchor=(0, 0), loc="lower left", borderaxespad=0)
    axL.legend(bbox_to_ancher=(1, -0.1), loc="upper right", borderaxespad=0)
    # axR.legend(loc="best")

    # TODO
    # 3つ目のグラフの作成

    if title == "":
        title = "model_history"
    fig.suptitle(title, fontsize=20)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    fig.savefig(os.path.join(dst, f"{title}.png"), bbox_inches="tight", pad_inches=0.05)
    fig.clf()


def plot_cmx(y_true: list, y_pred: list, labels: List[str], title: str, dst: PathLike):
    cmx = confusion_matrix(y_true, y_pred)
    normalized_cmx = [row / np.sum(row) for row in cmx]

    df_cmx = pd.DataFrame(cmx, index=labels, columns=labels)

    plt.figure(figsize=(15, 13))
    sns.heatmap(normalized_cmx,
                xticklabels=labels,
                yticklabels=labels,
                annot=False,
                square=True,
                cmap="Blues",
                vmin=0,
                vmax=1.0)
    plt.title("正解の綱と予測した綱の対応")
    plt.xlabel("予測した綱")
    plt.ylabel("正解の綱")

    # write value in cells
    data = df_cmx.values
    for y, x in itertools.product(range(data.shape[0]), range(data.shape[1])):
        if x == y:
            color = "red"
        elif normalized_cmx[y][x] > 0.5:
            color = "white"
        else:
            color = "black"
        plt.text(x + 0.5,
                 y + 0.5,
                 data[y, x],
                 horizontalalignment="center",
                 verticalalignment="center",
                 color=color)
    plt.savefig(Path(dst) / (title + ".png"), bbox_inches="tight", pad_inches=0.05)
    plt.clf()


def visualize_all_cmxs(rst: list, cmx_pathes: Iterable[PathLike],
                       dst: PathLike) -> None:
    images = []
    dst = Path(dst)
    # make table
    d = {}
    for row in rst:
        for k, v in row.items():
            d[k] = {kk: f"{vv:.3f}" for kk, vv in v.items()}
    df = pd.DataFrame.from_dict(d, orient="index")

    # concate table and cmxs
    for i, img_path in enumerate(cmx_pathes, 1):
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 4))
        ax1.axis("off")
        tb = ax1.table(cellText=df.values,
                       colLabels=df.columns,
                       rowLabels=df.index,
                       loc="center",
                       bbox=[0, 0, 1, 1])
        for j in range(len(df.columns)):
            tb[0, j].set_facecolor('#363636')
            tb[0, j].set_text_props(color='w')
            tb[i, j].set_facecolor("tomato")

        ax2.axis("off")
        img = plt.imread(img_path)
        ax2.imshow(img)
        filedst = dst / f".tmp_{i}.png"
        fig.savefig(filedst, bbox_inches="tight", pad_inches=0.05)
        images.append(Image.open(filedst))
        os.remove(filedst)

    # make gif animation
    images[0].save(dst / "log.gif",
                   format="gif",
                   save_all=True,
                   append_images=images[1:],
                   duration=1000,
                   loop=0)


if __name__ == "__main__":
    import sys
    import json
    filename = sys.argv[1]
    with open(filename) as f:
        d = json.load(f)
    project_dit = Path(__file__).resolve().parents[1]
    title = sys.argv[2]
    visualize_history(d, title, project_dit)
