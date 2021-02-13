import pandas as pd
import matplotlib.pyplot as plt
import japanize_matplotlib
import matplotlib as mpl
import seaborn as sns

from collections import Counter
import sys
from pathlib import Path
mpl.use("Agg")

if __name__ == "__main__":
    filename = Path(sys.argv[1])
    df = pd.read_csv(filename)
    df.rename(columns={"class": "cl"}, inplace=True)
    c = sorted(Counter(df["cl"]).items(), key=lambda x: -x[1])
    clnames = [x[0] for x in c]
    medians, means = [], []
    for clname in clnames:
        new_df = df.query("cl == @clname")
        medians.append(new_df.median()["at_gc_rate"])
        means.append(new_df.mean()["at_gc_rate"])

    holder = list(range(len(clnames)))
    sns.set()
    fig, (axL, axR) = plt.subplots(ncols=2, figsize=(10, 4))
    axL.barh(holder, means, align="center")
    axL.set_yticks(holder)
    axL.set_yticklabels(clnames, fontsize=3)
    axL.set_xlim([0, 4])
    axL.set_title("means")

    axR.barh(holder, medians, align="center")
    axR.set_yticks(holder)
    axR.set_yticklabels([""] * len(clnames), fontsize=3)
    axR.set_xlim([0, 4])
    axR.set_title("medians")

    fig.suptitle("AT/GC value of each class")
    fig.savefig(Path(__file__).resolve().parents[1] / "at_gc_rate.pdf",
                bbox_inches="tight",
                pad_inches=0.05)
