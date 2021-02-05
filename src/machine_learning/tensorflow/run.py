# 実際に機械学習を実行するプログラムを書く

import sys
from PIL import Image
from . import model, visualize, calc_loss_weight, focal_loss
import json
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import CSVLogger, History, TensorBoard
from sklearn.metrics import classification_report
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
        pathes (np.ndarray): 画像へのパスが格納されたnumpy配列

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


if __name__ == "__main__":
    settings = {
        "use_weight_method": "log",
        "use_csv": sys.argv[1],
        "normalize_seq_len": True,
        "KFold": 5,
        "description": "seq_lenを標準化する"
    }

    source_features = Path(sys.argv[1])

    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    dst = Path(config["destination"])

    logdst = logdir(config["log_dst"])

    df = pd.read_csv(source_features)

    labels = to_categorical(df["class"].values)

    weights = calc_loss_weight.calc_class_weight(np.argmax(labels, axis=1),
                                                 settings["use_weight_method"])
    images = load_images(map(lambda x: dst / "img" / f"{x}.png", df["accession"]))
    seq_lens = df["seq_len"].values.reshape([-1, 1])
    if settings["normalize_seq_len"]:
        scaler = StandardScaler()
        scaler.fit_transform(seq_lens)
        settings["seq_len_mean"] = scaler.mesn_
        settings["seq_len_var"] = scaler.var_

    n_class = len(set(np.argmax(labels, axis=1)))

    study_log = {}
    skf = StratifiedKFold(settings["KFold"])
    row_labels = df["class"].values
    for i, (train_index, test_index) in enumerate(skf.split(images, row_labels), 1):
        trial_dst = logdst / f"trial_{i}"
        trial_dst.mkdir()
        ml_model = model.construct_model(n_class)

        model.show_model(ml_model, trial_dst / "model.png")

        train_images, test_images = images[train_index], images[test_index]
        train_labels, test_labels = labels[train_index], labels[test_index]
        train_seq_lens, test_seq_lens = seq_lens[train_index], seq_lens[test_index]
        test_row_labels = row_labels[test_index]

        csv_log = CSVLogger(trial_dst / "logger.csv")
        tensor_board = TensorBoard(log_dir=trial_dst,
                                   write_graph=False,
                                   write_images=True,
                                   histogram_freq=1)
        history = History()

        if n_class == 2:
            ml_model.compile(loss=[focal_loss.binary_focal_loss(alpha=.25, gamma=2)],
                             metrics=["accuracy"],
                             optimizer=Adam())
        else:
            ml_model.compile(loss=[
                focal_loss.categorical_focal_loss(alpha=[[.25 for _ in range(n_class)]],
                                                  gamma=2)
            ],
                             metrics=["accuracy"],
                             optimizer=Adam())

        history = ml_model.fit([train_images, train_seq_lens],
                               train_labels,
                               validation_data=([test_images,
                                                 test_seq_lens], test_labels),
                               epochs=50,
                               callbacks=[csv_log, tensor_board],
                               class_weight=weights)
        study_log[f"trial_{i}"] = history.history.copy()

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

        with open(trial_dst / "report.txt", "w") as f:
            print(classification_report(test_labels,
                                        pred_labels,
                                        target_names=get_sorted_class(
                                            df["class"].values),
                                        zero_division=0),
                  file=f)

    with open(logdst / "weight.json", "w") as f:
        json.dump(
            {
                str(k): weights[i]
                for i, k in enumerate(get_sorted_class(df["class"].values))
            },
            f,
            indent=2)

    with open(logdst / "status.json", "w") as f:
        json.dump(settings, f, indent=2, ensure_ascii=False)
