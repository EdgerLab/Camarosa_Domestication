"""
1. Extracts AED scores from gene annotation files and saves them to a new file.
2. Removes genes that are not the primary gene model.

Output format is a tab-separated file with the following columns:
    - Gene_Name
    - Chromosome
    - Start
    - Stop
    - AED_Score
"""

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np
import argparse
import os
import logging
import coloredlogs

import matplotlib.pyplot as plt
import matplotlib.ticker

from src.orthologs.utils import (
    remove_str_from_val,
    drop_rows_with_bad_val_in_col,
    map_names,
)


def read_aed_output_table(filepath):
    return pd.read_csv(filepath, sep="\t", header="infer")


def plot_aed_score_distribution(aed_table, output_path, genome_name):
    counts, bins = np.histogram(
        aed_table["AED_Score"].to_numpy(np.float64), range=(0.0, 1.0), bins=25
    )

    fig1, ax1 = plt.subplots()
    ax1.stairs(counts, bins, fill=True)
    ax1.set_yscale("log")
    ax1.set_xlabel("AED Score")
    ax1.set_ylabel("Count of Genes")
    ax1.set_title("AED Score Distribution")
    plt.minorticks_off()

    if genome_name == "DN":
        ax1.set_yticks([100, 125, 250, 500, 1000, 5000, 10000, 15000])
    else:
        pass
        # No need to do anything special for RR
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    return plt, fig1


def get_AED(genes_input_path, genome_name, logger):
    """Import gene annotation."""

    col_names = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Score",
        "Strand",
        "Frame",
        "FullName",
    ]

    col_to_use = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Strand",
        "FullName",
    ]

    gene_pandaframe = pd.read_csv(
        genes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        dtype={"Stop": "float64", "Start": "float64", "Strand": str, "FullName": str},
        comment="#",
    )

    gene_pandaframe = gene_pandaframe[
        gene_pandaframe.Feature == "mRNA"
    ]  # drop non-gene rows in annotation

    # clean the names and set as the index (get row wrt name c.f. idx)
    gene_pandaframe["Gene_Name"] = gene_pandaframe["FullName"].str.extract(r"ID=(.*?);")

    # Define the AED score
    gene_pandaframe["AED_Score"] = gene_pandaframe["FullName"].str.extract(
        r"_AED=(.*?);"
    )

    # Remove genes that are not the primary gene model
    # Remove anything that doesn't end with '.1'

    if genome_name != "DN":
        # Only take the primary gene model
        gene_pandaframe = gene_pandaframe[
            gene_pandaframe["Gene_Name"].str.endswith(".1")
        ]
        # Reformat so that our gene names are nice
        # Remove the trailing '.1'
        gene_pandaframe["Gene_Name"] = gene_pandaframe["Gene_Name"].str.replace(
            ".1", "", regex=False
        )

    else:
        # Only take the primary gene model
        gene_pandaframe = gene_pandaframe[
            gene_pandaframe["Gene_Name"].str.endswith("-mRNA-1")
        ]
        # Reformat so that our gene names are nice
        # Remove the trailing -mRNA-1
        gene_pandaframe["Gene_Name"] = gene_pandaframe["Gene_Name"].str.replace(
            "-mRNA-1", "", regex=False
        )

    gene_pandaframe.set_index("Gene_Name", inplace=True)

    # Remove extaneous BS chromosomes
    for i in ["contig", "ptg", "scaf"]:
        gene_pandaframe = drop_rows_with_bad_val_in_col(
            gene_pandaframe, i, "Chromosome"
        )

    # Fix the chromosome names of the Del_Norte annotation so that they
    # correspond to the chromosomes of the TE annotation
    if genome_name == "DN":
        gene_pandaframe = remove_str_from_val(gene_pandaframe, "_RagTag", "Chromosome")
        gene_pandaframe = map_names(gene_pandaframe, "Chromosome")

    # Remove the Fvb prefix from the chr names, particularly in FNI, FVI, and H4
    gene_pandaframe = remove_str_from_val(gene_pandaframe, "Fvb", "Chromosome")

    # Remove the chr prefix from the chr names, particularly in RR, and FII
    gene_pandaframe = remove_str_from_val(gene_pandaframe, "chr_", "Chromosome")
    gene_pandaframe = remove_str_from_val(gene_pandaframe, "chr", "Chromosome")

    # Remove the _RagTag from the chromosome names
    gene_pandaframe = remove_str_from_val(gene_pandaframe, "_RagTag", "Chromosome")

    gene_pandaframe.drop(
        columns=["FullName", "Strand", "Feature", "Software", "Start", "Stop"],
        inplace=True,
    )

    gene_pandaframe.sort_index(inplace=True)

    return gene_pandaframe


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TODO")
    parser.add_argument(
        "gene_input_file", type=str, help="Parent path of gene annotation file"
    )
    parser.add_argument("genome_name", type=str)

    parser.add_argument(
        "output_graph",
        type=str,
        help="path for output graph",
    )
    parser.add_argument(
        "output_table",
        type=str,
        help="path for output file",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.gene_input_file = os.path.abspath(args.gene_input_file)
    args.output_table = os.path.abspath(args.output_table)
    args.output_graph = os.path.abspath(args.output_graph)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # ------------------------------------------------------------------

    cleaned_genes = get_AED(args.gene_input_file, args.genome_name, logger)
    logger.info(f"Saving AED table to {args.output_table}")
    cleaned_genes.to_csv(args.output_table, sep="\t", header=True, index=True)

    logger.info(f"Saving AED graph to {args.output_graph}")
    plt, fig = plot_aed_score_distribution(
        cleaned_genes, args.output_graph, args.genome_name
    )
    plt.tight_layout()
    plt.savefig(args.output_graph)
