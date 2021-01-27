import json
from argparse import ArgumentParser, Namespace

import matplotlib.pyplot as plt
from Bio import SeqIO

from .generate_img import generate_image

# python -m generate_one.py --weight <weight path> --inputs *.gbk


def parser() -> Namespace:
    """setting for argparser

    Returns:
        Namespace: args namespace.
    """
    usage = f"Usage: python {__file__} [-w weight] [-i inputs]"
    argparser = ArgumentParser(usage=usage)
    argparser.add_argument("-w",
                           "--weight",
                           help="If you already have weight data, write the path.")
    argparser.add_argument("-i", "--inputs", nargs="*", help="Input file pathes")

    args = argparser.parse_args()
    return args


if __name__ == "__main__":
    args = parser()

    with open(args.weight) as f:
        weight = json.load(f)

    for gbkfile in args.input():
        for record in SeqIO.parse(gbkfile, "genbank"):
            acc = record.name
            seq = record.seq.upper()
            fig = generate_image(seq, weight)
            dst = f"{acc}.png"
            plt.savefig()
            plt.close()
            print(f"{acc} have imaged to {dst} .")
