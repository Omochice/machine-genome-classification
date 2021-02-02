import matplotlib.pyplot as plt
from pathlib import Path

import itertools
import japanize_matplotlib
import matplotlib
from os import PathLike
import os
from sklearn.metrics import confusion_matrix

import numpy as np
import pandas as pd
import seaborn as sns

from typing import List

matplotlib.use("Agg")


def visualize_history(history: dict, title: str = "", dst: PathLike = "") -> None:
    fig, (axL, axR) = plt.subplots(ncols=2, figsize=(10, 4))

    axL.plot(history["acc"], "o-", label="Train accuracy")
    axL.plot(history["val_acc"], "o-", label="Validation accuracy")
    axL.set_title("Accuracy")
    axL.set_xlabel("Epoch")
    axL.set_ylabel("Accuracy")
    axL.set_ylim(0, 1)
    axL.grid(True)
    axL.legend(bbox_to_anchor=(0, 0), loc="lower left", borderaxespad=0)

    axR.plot(history["loss"], "o-", label="Train loss")
    axR.plot(history["val_loss"], "o-", label="Validation loss")
    axR.set_title("Loss")
    axR.set_xlabel("Epoch")
    axR.set_ylabel("Loss")
    axR.grid(True)
    axR.legend(bbox_to_anchor=(0, 0), loc="lower left", borderaxespad=0)

    if title == "":
        title = "model_history"
    fig.suptitle(title, fontsize=20)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    fig.savefig(os.path.join(dst, f"{title}.png"))
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
    plt.title("真のラベルと予測ラベルの対応")
    plt.xlabel("予測されたラベル")
    plt.ylabel("真のラベル")

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
    plt.savefig(Path(dst) / (title + ".png"))
    plt.clf()


if __name__ == "__main__":
    import sys
    import json
    filename = sys.argv[1]
    with open(filename) as f:
        d = json.load(f)
    project_dit = Path(__file__).resolve().parents[1]
    title = sys.argv[2]
    visualize_history(d, title, project_dit)
