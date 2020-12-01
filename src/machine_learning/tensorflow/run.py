# 実際に機械学習を実行するプログラムを書く

from . import model
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from pathlib import Path
import numpy as np
import pandas as pd

from os import PathLike


def load_data(filename: PathLike):
    table = pd.read_csv(filename, delimiter=",").to_numpy()
    names, classes, pathes, n_len, rate_a, rate_t, rate_g, rate_c = np.hsplit(
        table, table.shape[1])

    def _to1dim(arr: np.ndarray) -> np.ndarray:
        return arr.reshape(-1)

    return tuple(
        _to1dim(x) for x in (names, classes, pathes, n_len, rate_a, rate_t, rate_g, rate_t))


if __name__ == "__main__":
    project_dir = ""
    fname = ""
    accs, labels, pathes, n_len, a, t, g, c = load_data(fname)
    base_rates = np.stack([a, t, g, c], axis=1)
    ss = StandardScaler()
    n_len = ss.fit_transform(n_len.reshape(n_len.shape[0], 1))
