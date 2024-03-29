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

from src.go_analysis.upset_plot import (
    read_go_enrichment_table,
    get_nonredundant_terms,
    decode_te_window_direction_str_go_file,
)

"""
Code for generating a DN vs RR UpSet plot, this code takes from upset_plot.py
because I didn't have time to refactor to remove the H4 code.

This code will generate UpSet plots for the DN vs RR DIFFERNCE GO enrichment
files. That data is the set of DN-RR syntelog pairs that are biased towards one
genome or the other. So for example the file:
    'Overrepresented_Difference_Total_TE_Density_1000_Upstream_BiasedTowards_DN_95_density_percentile.tsv'
    is the file that has all of the genes that have a high TE density for the
    DN syntelog but almost nothing for the RR syntelog. These are the big TE
    PAV files.
    We anticipate that few GO terms will be shared within the UpSet plot
"""

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
        "output_dir",
        type=str,
        help="parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    parser.add_argument(
        "--syntelog", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.DN_preprocessed_go_enrichment_table = os.path.abspath(
        args.DN_preprocessed_go_enrichment_table
    )
    args.RR_preprocessed_go_enrichment_table = os.path.abspath(
        args.RR_preprocessed_go_enrichment_table
    )
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # ---------------------------
    # Read in the data
    DN_table = read_go_enrichment_table(args.DN_preprocessed_go_enrichment_table)
    RR_table = read_go_enrichment_table(args.RR_preprocessed_go_enrichment_table)

    # NOTE
    # Plot title is something like "Mutator 2500 Upstream Upper 95 Percentile"
    # Hopefully the user is only comparing DN and RR preprocessed GO enrichment
    # files that are of the same TE type, direction, and cutoff. As we are
    # basing the name off of the DN file, we can assume that the RR file is the
    # same type, direction, and cutoff.
    name = decode_te_window_direction_str_go_file(
        args.DN_preprocessed_go_enrichment_table
    )
    # Shorten the Name because we are performing this code on the DN vs RR
    # 'difference' files, and the string rule doesn't work out of the box with
    # those filenames

    # TODO the DIRECTION IS NOT IN THE FILENAME

    # Get the set of terms
    RR_terms = get_nonredundant_terms(RR_table)
    DN_terms = get_nonredundant_terms(DN_table)

    # Get the shared terms, there should be few, and write them to a file
    shared_terms = RR_terms.intersection(DN_terms)
    shared_table_subset = DN_table.loc[DN_table["GO_ID"].isin(shared_terms)]
    shared_table_subset = shared_table_subset.loc[:, ["GO_ID", "Term"]]
    shared_table_subset.drop_duplicates(inplace=True)

    if args.syntelog:
        name = name[:-2]
        filename = "_".join(name) + "_Syntelog"
        shared_table_out = os.path.join(
            args.output_dir,
            str("Shared_Terms_" + os.path.basename(filename) + ".tsv"),
        )
    else:
        name = name[:-3]
        filename = "_".join(name) + "_NonSyntelog"
        shared_table_out = os.path.join(
            args.output_dir,
            str("Shared_Terms_" + os.path.basename(filename) + ".tsv"),
        )
    filename += ".png"
    plot_title = " ".join(name)
    logger.info(f"Saving shared terms to {shared_table_out}")
    shared_table_subset.to_csv(shared_table_out, sep="\t", index=False)

    # Get the count of unique genes in the GO enrichment table
    RR_gene_count = RR_table["RR_Gene"].nunique()
    DN_gene_count = DN_table["DN_Gene"].nunique()

    # Generate the data structure for the UpSet plot
    data = {"RR": RR_terms, "DN": DN_terms}
    data = from_contents(data)

    # Initialize the plot figure object
    fig = plt.figure(figsize=(10, 9))
    upset = UpSet(data, subset_size="count", show_counts=True)
    upset.style_subsets(present="RR", absent=["DN"], facecolor="red")
    upset.style_subsets(present="DN", absent=["RR"], facecolor="blue")
    upset.style_subsets(present=["RR", "DN"], hatch="xx", facecolor="purple")
    upset.plot(fig=fig)
    plt.suptitle(plot_title)
    fig.supxlabel("Number of GO Terms")

    # Add the gene count to the plot
    # TODO work on this, it is getting too crowded
    gene_count_handles = []
    for i in [
        (RR_gene_count, "red", "RR"),
        (DN_gene_count, "blue", "DN"),
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
