# 実際に機械学習を実行するプログラムを書く

from PIL import Image
from . import model
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow.keras.preprocessing.image import load_img
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

from os import PathLike


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


def load_images(pathes: np.ndarray) -> np.ndarray:
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
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)
    fname = Path(config["data_destionation"]) / "csv" / "features.csv"
    accs, labels, pathes, n_len, a, t, g, c = load_data(fname)
    ss = StandardScaler()
    n_len = ss.fit_transform(n_len.reshape(n_len.shape[0], 1)).flatten()
    print(n_len.shape, a.shape)
    liner_fetures = np.stack([a, t, g, c, n_len], axis=1).astype("float32")
    print(liner_fetures)
    print(f"mean: {ss.mean_}, var: {ss.var_}")
    labels = to_categorical(labels)
    images = load_images(pathes)

    ml_model = model.construct_model(n_class=len(labels[0]))
    ml_model.compile(optimizer="adam",
                     loss="categorical_crossentropy",
                     metrics=["accuracy"])
    ml_model.fit({
        "input1": images,
        "input2": liner_fetures
    },
        labels,
        epochs=50,
        batch_size=64,
        validation_split=0.2)
