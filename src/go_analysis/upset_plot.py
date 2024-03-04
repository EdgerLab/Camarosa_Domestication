__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd

from upsetplot import from_contents
from upsetplot import UpSet
import matplotlib.pyplot as plt


def read_go_enrichment_table(path):
    """
    Read in the preprocessed GO enrichment table
    """
    # TODO think about refactoring this to the find_abnormal_genes.py script
    df = pd.read_csv(path, sep="\t", header="infer")
    df = df.loc[:, ["GO_ID"]]
    df = df.drop_duplicates()
    return df


def read_go_enrichment_table(path):
    """
    Read in the preprocessed GO enrichment table
    """
    # TODO think about refactoring this to the find_abnormal_genes.py script
    df = pd.read_csv(path, sep="\t", header="infer")
    df = df.loc[:, ["GO_ID"]]
    df = df.drop_duplicates()
    return df


def decode_te_window_direction_str_go_file(path):
    # MAGIC: This is a hacky way to get the name of the TE window direction and
    # cutoff
    name = os.path.basename(args.DN_preprocessed_go_enrichment_table).split("_")[2:-2]
    # name[1] = name[1] + " BP "  # add 'BP' to the window value
    name.append("Percentile")
    return name


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "DN_preprocessed_go_enrichment_table",
        type=str,
        help="TODO",
    )

    parser.add_argument(
        "RR_preprocessed_go_enrichment_table",
        type=str,
        help="TODO",
    )

    parser.add_argument(
        "H4_preprocessed_go_enrichment_table",
        type=str,
        help="TODO",
    )

    parser.add_argument(
        "output_dir",
        type=str,
        help="parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.DN_preprocessed_go_enrichment_table = os.path.abspath(
        args.DN_preprocessed_go_enrichment_table
    )
    args.RR_preprocessed_go_enrichment_table = os.path.abspath(
        args.RR_preprocessed_go_enrichment_table
    )
    args.H4_preprocessed_go_enrichment_table = os.path.abspath(
        args.H4_preprocessed_go_enrichment_table
    )
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # ---------------------------
    # Read in the data
    DN_table = read_go_enrichment_table(args.DN_preprocessed_go_enrichment_table)
    RR_table = read_go_enrichment_table(args.RR_preprocessed_go_enrichment_table)
    H4_table = read_go_enrichment_table(args.H4_preprocessed_go_enrichment_table)

    # NOTE
    # Plot title is something like "Mutator 2500 Upstream Upper 95 Percentile"
    # Hopefully the user is only comparing DN and RR preprocessed GO enrichment
    # files that are of the same TE type, direction, and cutoff. As we are
    # basing the name off of the DN file, we can assume that the RR file is the
    # same type, direction, and cutoff.
    name = decode_te_window_direction_str_go_file(
        args.DN_preprocessed_go_enrichment_table
    )
    plot_title = " ".join(name)
    filename = "_".join(name)

    RR_terms = RR_table["GO_ID"].values
    DN_terms = DN_table["GO_ID"].values
    H4_terms = DN_table["GO_ID"].values

    data = from_contents({"RR": RR_terms, "DN": DN_terms, "H4": H4_terms})

    # Initialize the plot figure object
    fig = plt.figure()
    upset = UpSet(data, subset_size="count")
    upset.style_subsets(present="RR", absent=["DN", "H4"], facecolor="red")
    upset.style_subsets(present="DN", absent=["RR", "H4"], facecolor="blue")
    upset.style_subsets(present="H4", absent=["DN", "RR"], facecolor="green")
    upset.style_subsets(present=["RR", "DN"], hatch="xx", edgecolor="purple")
    upset.style_subsets(present=["H4", "DN"], hatch="xxx", edgecolor="yellow")
    upset.plot(fig=fig)
    plt.suptitle(plot_title)
    plt.show()
    raise ValueError("stop here")
    plt.savefig(os.path.join(args.output_dir, filename), bbox_inches="tight")
    plt.clf()
