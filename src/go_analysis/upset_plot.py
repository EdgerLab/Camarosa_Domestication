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


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "RR_preprocessed_go_enrichment_table",
        type=str,
        help="TODO",
    )
    parser.add_argument(
        "DN_preprocessed_go_enrichment_table",
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
    args.RR_preprocessed_go_enrichment_table = os.path.abspath(
        args.RR_preprocessed_go_enrichment_table
    )
    args.DN_preprocessed_go_enrichment_table = os.path.abspath(
        args.DN_preprocessed_go_enrichment_table
    )
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # ---------------------------
    # Read in the data
    RR_table = read_go_enrichment_table(args.RR_preprocessed_go_enrichment_table)
    DN_table = read_go_enrichment_table(args.DN_preprocessed_go_enrichment_table)

    RR_terms = RR_table["GO_ID"].values
    DN_terms = DN_table["GO_ID"].values

    data = from_contents({"RR": RR_terms, "DN": DN_terms})


def read_go_enrichment_table(path):
    """
    Read in the preprocessed GO enrichment table
    """
    # TODO think about refactoring this to the find_abnormal_genes.py script
    df = pd.read_csv(path, sep="\t", header="infer")
    df = df.loc[:, ["GO_ID"]]
    df = df.drop_duplicates()
    return df


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "RR_preprocessed_go_enrichment_table",
        type=str,
        help="TODO",
    )
    parser.add_argument(
        "DN_preprocessed_go_enrichment_table",
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
    args.RR_preprocessed_go_enrichment_table = os.path.abspath(
        args.RR_preprocessed_go_enrichment_table
    )
    args.DN_preprocessed_go_enrichment_table = os.path.abspath(
        args.DN_preprocessed_go_enrichment_table
    )
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # ---------------------------
    # Read in the data
    RR_table = read_go_enrichment_table(args.RR_preprocessed_go_enrichment_table)
    DN_table = read_go_enrichment_table(args.DN_preprocessed_go_enrichment_table)

    RR_terms = RR_table["GO_ID"].values
    DN_terms = DN_table["GO_ID"].values

    data = from_contents({"RR": RR_terms, "DN": DN_terms})

    # Initialize the plot figure object
    fig = plt.figure()
    ax_dict = UpSet(data, subset_size="count").plot(fig=fig)
    plt.show()
