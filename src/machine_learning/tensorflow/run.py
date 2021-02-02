# 実際に機械学習を実行するプログラムを書く

import sys
from PIL import Image
from . import model, visualize
import json
import sklearn
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow.keras.preprocessing.image import load_img
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import CSVLogger, History, TensorBoard
from sklearn.utils.class_weight import compute_class_weight
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

from typing import Iterable, Optional
from os import PathLike

import datetime

project_dir = Path(__file__).resolve().parents[3]


def logdir(root: PathLike) -> Path:
    """Make logdir and return it.

    Args:

    Returns:
        Path: logdir
    """
    now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
    logdir = Path(root) / now
    logdir.mkdir(parents=True, exist_ok=True)
    return logdir


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

    # def _get_sorted_class(classes: np.ndarray) -> np.ndarray:
    #     """昇順でソートしたユニークな配列を返す
    #
    #     Args:
    #         classes (np.ndarray): 対象のnumpy配列
    #
    #     Returns:
    #         np.ndarray:
    #     """
    #     cl, counts = np.unique(classes, return_counts=True)
    #     cl_with_counts = np.stack([cl, counts], axis=1)
    #     sorted_classes = cl_with_counts[np.argsort(cl_with_counts[:, 1])]
    #     return sorted_classes[:, 0]

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

    uniq = get_sorted_class(labels)
    numerical_labels = _label_to_number(labels, uniq)
    if not numerical_labels.ndim == 1:
        numerical_labels = numerical_labels.flatten()
    one_hot = np.identity(np.max(numerical_labels) + 1)[numerical_labels]
    return one_hot


def get_sorted_class(classes: np.ndarray) -> np.ndarray:
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


def load_images(pathes: Iterable[PathLike]) -> np.ndarray:
    """画像のパスが格納されたnumpy配列から画素値の配列を得る

    Args:
        pathes (np.ndarray): 画像へのパスが格納されたnumpt配列

    Returns:
        np.ndarray:
    """
    return np.array([(np.array(Image.open(p).convert("L")).astype("float32")) / 255
                     for p in pathes])


def normalize(seq_lens: np.ndarray, sigma: int = 0.3) -> np.ndarray:
    average = np.mean(seq_lens / 1000)
    rst = np.vectorize(lambda x: np.exp(-1 * (((average - x)**2) / (2 * sigma**2))))(
        seq_lens / 1000)
    return rst


# def label_to_rester_number(labels: np.ndarray, roster: dict):
#     for i, label in enumerate(labels):
#         if label in roster["no1"]:
#             labels[i] = 0
#         elif label in roster["no2"]:
#             labels[i] = 1
#         elif label in roster["no3"]:
#             labels[i] = 2
#         else:
#             labels[i] = 3
#     return labels
#
#
# def extract_df(df, use_roster, roster) -> pd.DataFrame:
#     if use_roster == "all":
#         return df
#     elif use_roster == "to_4":
#         return df
#     elif use_roster == "no3":
#         return df[df["class"].isin(set(roster["no3"]))]
#     elif use_roster == "no4":
#         return df[df["class"].isin(set(roster["no4"]))]
#

if __name__ == "__main__":
    settings = {
        "use_weight_method": "log",
        "use_csv": sys.argv[1],
        "normalize_seq_len": True,
        "KFold": 5
    }

    source_features = Path(sys.argv[1])

    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    dst = Path(config["destination"])

    logdst = logdir(config["log_dst"])

    df = pd.read_csv(source_features)

    labels = to_categorical(df["class"].values)

    weights = calc_label_weights(np.argmax(labels, axis=1),
                                 settings["use_weight_method"])
    images = load_images(map(lambda x: dst / "img" / f"{x}.png", df["accession"]))
    seq_lens = df["seq_len"].values.reshape([-1, 1])
    if settings["normalize_seq_len"]:
        seq_lens = normalize(seq_lens)
    n_class = len(set(np.argmax(labels, axis=1)))

    study_log = {}
    skf = StratifiedKFold(settings["KFold"])
    for i, (data_index, test_index) in enumerate(skf.split(images, df["class"].values),
                                                 1):
        trial_dst = logdst / f"trial_{i}"
        trial_dst.mkdir()
        ml_model = model.construct_model(n_class)

        model.show_model(ml_model, trial_dst / "model.png")

        (train_images, test_images, train_seq_lens, test_seq_lens, train_labels,
         test_labels) = train_test_split(images[data_index],
                                         seq_lens[data_index],
                                         labels[data_index],
                                         test_size=0.2,
                                         stratify=labels[data_index])

        csv_log = CSVLogger(trial_dst / "logger.csv")
        tensor_board = TensorBoard(log_dir=trial_dst,
                                   write_graph=True,
                                   write_images=True,
                                   histogram_freq=1)
        history = History()
        ml_model.compile(loss="categorical_crossentropy",
                         optimizer=Adam(),
                         metrics=["acc"])
        history = ml_model.fit([train_images, train_seq_lens],
                               train_labels,
                               validation_data=([test_images,
                                                 test_seq_lens], test_labels),
                               epochs=50,
                               callbacks=[csv_log, tensor_board],
                               class_weight=weights)

        visualize.visualize_history(history.history, "study_log", trial_dst)

        loss, acc = ml_model.evaluate([test_images, test_seq_lens],
                                      test_labels,
                                      verbose=1)
        pred_labels = np.argmax(ml_model.predict([test_images, test_seq_lens]), axis=1)
        test_labels = np.argmax(test_labels, axis=1)

        visualize.plot_cmx(test_labels,
                           pred_labels,
                           get_sorted_class(df["class"].values),
                           title="cmx",
                           dst=trial_dst)

        with open(trial_dst / "weight.json", "w") as f:
            json.dump(
                {
                    str(k): weights[i]
                    for i, k in enumerate(get_sorted_class(df["class"].values))
                },
                f,
                indent=2)

    with open(logdst / "status.json", "w") as f:
        json.dump(settings, f, indent=2)
