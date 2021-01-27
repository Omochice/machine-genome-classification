# 実際に機械学習を実行するプログラムを書く

import sys
from PIL import Image
from . import model, visualize
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow.keras.preprocessing.image import load_img
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import CSVLogger, History
from sklearn.utils.class_weight import compute_class_weight
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

from typing import Iterable, Optional
from os import PathLike


def calc_label_weights(labels: np.ndarray, option: Optional[str] = None) -> dict:
    """calc weight for imbalanced data.

    Args:
        labels (np.ndarray): labels. ex) [0, 1, 1, 1, ...]
        option (Optional[str]): select type option. in {None, "log"}
    Returns:
        dict: number to weight. ex) {0: 30, 1: 5, 2: 0.5}
    """
    if option not in {None, "log"}:
        print("Argment 'option' must be in {None, 'log'}")
        sys.exit(1)

    weights = compute_class_weight("balanced", np.unique(labels), labels)
    n_label_type = len(weights)
    if option == "log":
        return {
            i: weight
            for i, weight in enumerate(-np.log(1 / (n_label_type * weights)))
        }
    else:
        return {i: weight for i, weight in enumerate(weights)}


def load_data(filename: PathLike):
    """csvから特徴量を読み込む.

    Args:
        filename (PathLike): path to csv file.
    """
    table = pd.read_csv(filename, delimiter=",").to_numpy()
    names, classes, pathes, n_len, rate_a, rate_t, rate_g, rate_c = np.hsplit(
        table, table.shape[1])

    def _to1dim(arr: np.ndarray) -> np.ndarray:
        return arr.reshape(-1)

    return tuple(
        _to1dim(x)
        for x in (names, classes, pathes, n_len, rate_a, rate_t, rate_g, rate_t))


def to_categorical(labels: np.ndarray) -> np.ndarray:
    """labelをone-hotなラベルに変更する
    labelの数値はラベルの出現数の昇順になる
    ex) input  ['A', 'A', 'B', 'C', 'C', 'C']
        output [[0,1,0],[0,1,0],[1,0,0],[0,0,1],[0,0,1],[0,0,1]]

    Args:
        labels (np.ndarray): one-hotに変換するnumpy配列

    Returns:
        np.ndarray:
    """
    def _get_sorted_class(classes: np.ndarray) -> np.ndarray:
        """昇順でソートしたユニークな配列を返す

        Args:
            classes (np.ndarray): 対象のnumpy配列

        Returns:
            np.ndarray:
        """
        cl, counts = np.unique(classes, return_counts=True)
        cl_with_counts = np.stack([cl, counts], axis=1)
        sorted_classes = cl_with_counts[np.argsort(cl_with_counts[:, 1])]
        return sorted_classes[:, 0]

    def _label_to_number(labels: np.ndarray, uniq: np.ndarray) -> np.ndarray:
        """ユニークな配列を基準に数値の配列を返す

        Args:
            labels (np.ndarray): 数値に変換するnumpy配列
            uniq (np.ndarray): 基準となるユニークな配列

        Returns:
            np.ndarray:
        """
        numbers = np.array([np.where(ll == uniq) for ll in labels])
        return numbers

    uniq = _get_sorted_class(labels)
    numerical_labels = _label_to_number(labels, uniq)
    if not numerical_labels.ndim == 1:
        numerical_labels = numerical_labels.flatten()
    one_hot = np.identity(np.max(numerical_labels) + 1)[numerical_labels]
    return one_hot


def load_images(pathes: Iterable[PathLike]) -> np.ndarray:
    """画像のパスが格納されたnumpy配列から画素値の配列を得る

    Args:
        pathes (np.ndarray): 画像へのパスが格納されたnumpt配列

    Returns:
        np.ndarray:
    """
    return np.array([(np.array(Image.open(p).convert("L")).astype("float32")) / 255
                     for p in pathes])


if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parents[3]
    source_features = Path(sys.argv[1])

    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    df = pd.read_csv(source_features)

    labels = to_categorical(df["class"].values)
    weights = calc_label_weights(np.argmax(labels, axis=1))
    images = load_images(
        map(lambda x: project_dir / "data" / "img" / f"{x}.png", df["accession"]))
    seq_lens = df["seq_len"].values.reshape([-1, 1])
    n_class = len(set(np.argmax(labels, axis=1)))

    ml_model = model.construct_model(n_class)

    model.show_model(ml_model, "model.png")

    (train_images, test_images, train_seq_lens, test_seq_lens, train_labels,
     test_labels) = train_test_split(images, seq_lens, labels)

    csv_log = CSVLogger("logger.csv")
    history = History()
    ml_model.compile(loss="categorical_crossentropy", optimizer=Adam(), metrics=["acc"])
    history = ml_model.fit([train_images, train_seq_lens],
                           train_labels,
                           epochs=50,
                           class_weight=weights)
    visualize.visualize_history(history.history, "test", project_dir)
