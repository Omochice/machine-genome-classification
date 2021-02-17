# 実際に機械学習を実行するプログラムを書く

import sys
from PIL import Image
from . import model, visualize, calc_loss_weight, focal_loss, utils, callbacks
import json
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler, minmax_scale
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import CSVLogger, History, TensorBoard, ModelCheckpoint
from sklearn.metrics import classification_report, f1_score
from tensorflow.keras.utils import model_to_dot
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


def main(settings: dict, log_dst: PathLike):
    source_features = Path(settings["use_csv"])

    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)

    dst = Path(config["destination"])

    logdst = logdir(log_dst)

    df = pd.read_csv(source_features)
    n_data = len(df)

    labels = to_categorical(df["class"].values)

    weights = calc_loss_weight.calc_class_weight(np.argmax(labels, axis=1),
                                                 settings["use_weight_method"])
    images = load_images(map(lambda x: dst / "img" / f"{x}.png", df["accession"]))
    seq_lens = df["seq_len"].values.reshape([-1, 1])
    at_gc_rates = df["at_gc_rate"].values.reshape([-1, 1])
    atgc = np.stack([df["A"].values, df["T"].values, df["G"].values, df["T"].values]).T

    if settings["seq_len"]["enable"]:
        if settings["seq_len"]["rescale"]:
            if settings["seq_len"]["way"] == "standardization":
                scaler = StandardScaler()
                scaler.fit_transform(seq_lens)
                settings["seq_len"]["mean"] = float(scaler.mean_)
                settings["seq_len"]["var"] = float(scaler.var_)
            elif settings["seq_len"]["way"] == "normalization":
                seq_lens = minmax_scale(seq_lens)
                settings["seq_len"]["max"] = int(max(seq_lens))
                settings["seq_len"]["min"] = int(min(seq_lens))
    else:
        seq_lens = np.zeros(seq_lens.shape)

    if settings["atgc_rate"]["enable"]:
        if settings["atgc_rate"]["rescale"]:
            if settings["atgc_rate"]["way"] == "standardization":
                scaler = StandardScaler()
                scaler.fit_transform(at_gc_rates)
                settings["atgc_rate"]["mean"] = float(scaler.mean_)
                settings["atgc_rate"]["var"] = float(scaler.var_)
            elif settings["atgc_rate"]["way"] == "normalization":
                at_gc_rates = minmax_scale(at_gc_rates)
                settings["atgc_rate"]["max"] = float(max(at_gc_rates))
                settings["atgc_rate"]["min"] = float(min(at_gc_rates))
    else:
        at_gc_rates = np.zeros(at_gc_rates.shape)

    if settings["atgc"]["enable"]:
        if settings["atgc"]["rescale"]:
            if settings["atgc"]["way"] == "standardization":
                scaler = StandardScaler()
                scaler.fit_transform(atgc)
                settings["atgc"]["mean"] = float(scaler.mean_)
                settings["atgc"]["var"] = float(scaler.var_)
            elif settings["atgc"]["way"] == "normalization":
                atgc = minmax_scale(atgc)
                settings["atgc"]["max"] = float(max(atgc))
                settings["atgc"]["min"] = float(min(atgc))
    else:
        atgc = np.zeros(atgc.shape)

    n_class = len(labels[0])

    skf = StratifiedKFold(settings["KFold"])
    raw_labels = df["class"].values
    study_log = {"acc": [], "f1": []}
    for i, (train_index,
            test_index) in enumerate(skf.split(images, df["raw_class"].values), 1):
        trial_dst = logdst / f"trial_{i}"
        trial_dst.mkdir()
        ml_model = model.construct_model(n_class)

        model.show_model(ml_model, trial_dst / "model.pdf")
        model_to_dot(ml_model, show_shapes=True).write(str(trial_dst / "model.svg"),
                                                       format="svg")

        train_images, test_images = images[train_index], images[test_index]
        train_labels, test_labels = labels[train_index], labels[test_index]
        train_seq_lens, test_seq_lens = seq_lens[train_index], seq_lens[test_index]
        train_at_gc_rates, test_at_gc_rates = at_gc_rates[train_index], at_gc_rates[
            test_index]
        train_atgc, test_atgc = atgc[train_index], atgc[test_index]
        test_raw_labels = raw_labels[test_index]

        csv_log = CSVLogger(trial_dst / "logger.csv")
        tensor_board = TensorBoard(log_dir=trial_dst,
                                   write_graph=False,
                                   write_images=True,
                                   histogram_freq=1)
        f1cb = callbacks.F1Callback_(
            ml_model, [test_images, test_seq_lens, test_atgc, test_at_gc_rates],
            test_labels)
        history = History()
        checkpoint = ModelCheckpoint(logdst.parent / ".tmp_models" /
                                     "model_{epoch:02d}.hdf5",
                                     monitor="val_loss",
                                     save_weights_only=True)

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
        # if n_class == 2:
        #     ml_model.compile(loss="binary_crossentropy",
        #                      metrics=["accuracy"],
        #                      optimizer=Adam())
        # else:
        #     ml_model.compile(loss="categorical_crossentropy",
        #                      metrics=["accuracy"],
        #                      optimizer=Adam())

        history = ml_model.fit(
            [train_images, train_seq_lens, train_atgc, train_at_gc_rates],
            train_labels,
            validation_data=([test_images, test_seq_lens, test_atgc,
                              test_at_gc_rates], test_labels),
            epochs=100,
            batch_size=settings["batch"],
            callbacks=[csv_log, tensor_board, f1cb, checkpoint],
            class_weight=weights)
        # study_log[f"trial_{i}"] = history.history.copy()

        # TODO
        # f1scoreもプロットする
        # f1s = f1cb.f1s
        history.history["f1score"] = f1cb.f1s

        visualize.visualize_history(history.history, "study_log", trial_dst)

        # load best f1 score model
        ml_model.load_weights(logdst.parent / ".tmp_models" /
                              f"model_{np.argmax(f1cb.f1s) + 1:02}.hdf5")

        loss, acc, *_ = ml_model.evaluate(
            [test_images, test_seq_lens, test_atgc, test_at_gc_rates],
            test_labels,
            verbose=1)
        pred_labels = np.argmax(ml_model.predict(
            [test_images, test_seq_lens, test_atgc, test_at_gc_rates]),
                                axis=1)
        test_labels = np.argmax(test_labels, axis=1)
        settings["results"].append({
            f"trial_{i}": {
                "Accuracy": float(acc),
                "F1 score": float(max(f1cb.f1s))
        # "micro_f1": float(f1_score(test_labels, pred_labels, average="micro"))
            }
        })
        study_log["acc"].append(float(acc))
        study_log["f1"].append(float(max(f1cb.f1s)))

        visualize.plot_cmx(test_labels,
                           pred_labels,
                           get_sorted_class(df["class"].values),
                           title="cmx",
                           dst=trial_dst)

        with open(trial_dst / "report.txt", "w") as f:
            print(classification_report(test_labels,
                                        pred_labels,
                                        target_names=get_sorted_class(raw_labels),
                                        zero_division=0),
                  file=f)

    with open(logdst / "weight.json", "w") as f:
        json.dump(
            {str(k): weights[i]
             for i, k in enumerate(get_sorted_class(raw_labels))},
            f,
            indent=2)

    settings["average"] = {
        "Accuracy": sum(study_log["acc"]) / len(study_log["acc"]),
        "F1 score": sum(study_log["f1"]) / len(study_log["f1"])
    }

    with open(logdst / "status.json", "w") as f:
        json.dump(settings, f, indent=2, ensure_ascii=False)

    visualize.visualize_all_cmxs(
        settings["results"],
        [logdst / f"trial_{i}" / "cmx.png"
         for i in range(1, settings["KFold"] + 1)], logdst)


if __name__ == "__main__":
    settings = {
        "use_weight_method": None,
        "use_csv": sys.argv[1],
        "seq_len": {
            "enable": False,
            "rescale": True,
            "way": "normalization"
        },
        "atgc_rate": {
            "enable": True,
            "rescale": True,
            "way": "normalization"
        },
        "atgc": {
            "enable": False,
            "rescale": False,
            "way": None
        },
        "KFold": 5,
        "batch": 32,
        "description": "v1 逆数focal AT/GCのみ",
        "results": [],
        "average": {}
    }
    with open(project_dir / "setting.yml") as f:
        config = yaml.safe_load(f)
    log_dst = logdir(config["log_dst"])

    main(settings, log_dst)