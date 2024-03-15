__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd

from upsetplot import from_contents
from upsetplot import from_memberships
from upsetplot import UpSet
import matplotlib.pyplot as plt


def read_go_enrichment_table(path):
    """
    Read in the preprocessed GO enrichment table
    """
    # TODO think about refactoring this to the find_abnormal_genes.py script
    df = pd.read_csv(path, sep="\t", header="infer")
    return df


def get_nonredundant_terms(df):
    """Return as a set"""
    x = df.loc[:, ["GO_ID"]]

    # TODO check if we need to change any options in the drop duplicates
    # function
    x.drop_duplicates(inplace=True)

    return set(x["GO_ID"].values)


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
        "H4_preprocessed_go_enrichment_table",
        type=str,
        help="TODO",
    )

    parser.add_argument(
        "RR_preprocessed_go_enrichment_table",
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

    # TODO implement some sort of check to ensure that the input files have the
    # same string filename structure

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

    # print(args.DN_preprocessed_go_enrichment_table)
    # print(args.RR_preprocessed_go_enrichment_table)
    # print(args.H4_preprocessed_go_enrichment_table)
    # print(filename)
    # print()
    # print()
    # raise ValueError

    # Get the set of terms
    RR_terms = get_nonredundant_terms(RR_table)
    DN_terms = get_nonredundant_terms(DN_table)
    H4_terms = get_nonredundant_terms(H4_table)

    # Get the count of unique genes in the GO enrichment table
    RR_gene_count = RR_table["RR_Gene"].nunique()
    DN_gene_count = DN_table["DN_Gene"].nunique()
    H4_gene_count = H4_table["H4_Gene"].nunique()

    print(RR_table)
    print(len(RR_terms))
    print(RR_gene_count)
    print()

    # Generate the data structure for the UpSet plot
    data = {"RR": RR_terms, "DN": DN_terms, "H4": H4_terms}
    data = from_contents(data)

    # Initialize the plot figure object
    fig = plt.figure(figsize=(10, 9))
    upset = UpSet(data, subset_size="count", show_counts=True)
    upset.style_subsets(present="RR", absent=["DN", "H4"], facecolor="red")
    upset.style_subsets(present="DN", absent=["RR", "H4"], facecolor="blue")
    upset.style_subsets(present="H4", absent=["DN", "RR"], facecolor="green")
    # upset.style_subsets(present="RR", facecolor="red")
    # upset.style_subsets(present="DN", facecolor="blue")
    # upset.style_subsets(present="H4", facecolor="green")
    upset.style_subsets(present=["RR", "DN"], hatch="xx", facecolor="yellow")
    upset.style_subsets(present=["H4", "DN"], hatch="xx", facecolor="cyan")
    upset.style_subsets(present=["H4", "RR"], hatch="xx", facecolor="brown")
    upset.style_subsets(present=["DN", "H4", "RR"], hatch="xx", facecolor="black")
    upset.plot(fig=fig)
    plt.suptitle(plot_title)
    fig.supxlabel("Number of GO Terms")

    # Add the gene count to the plot
    # TODO work on this, it is getting too crowded
    gene_count_handles = []
    for i in [
        (RR_gene_count, "red", "RR"),
        (DN_gene_count, "blue", "DN"),
        (H4_gene_count, "green", "H4"),
    ]:
        gene_count_handles.append(
            plt.Line2D(
                [],
                [],
                marker="s",
                color=i[1],
                linestyle="",
                label=f"{i[2]} Unique Genes: {i[0]}",
            )
        )
    fig.legend(
        handles=gene_count_handles,
        loc="upper left",
        fontsize="small",
        bbox_to_anchor=(0.9, 0.9),
    )

    # plt.show()
    # plt.tight_layout()
    outfile = os.path.join(args.output_dir, filename)
    logger.info(f"Saving plot to {outfile}")
    plt.savefig(outfile, bbox_inches="tight")
    plt.clf()
