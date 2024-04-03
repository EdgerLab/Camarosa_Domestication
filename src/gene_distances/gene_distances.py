"""
Calculate distances between genes in the genome.
    1. Calculate the average distance between genes in the genome.
    2. Create a bar graph of the distances between genes.
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs

import matplotlib.pyplot as plt

from transposon.import_filtered_genes import import_filtered_genes

# from transposon.gene_data import GeneData


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "cleaned_genes",
        type=str,
        help="Parent path of the cleaned gene annotation file",
    )
    parser.add_argument("genome_name", type=str, help="Name of the genome")
    parser.add_argument(
        "output_dir",
        type=str,
        help="Parent directory to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.cleaned_genes = os.path.abspath(args.cleaned_genes)
    args.output_dir = os.path.abspath(args.output_dir)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -------------------------------------------------------------- #
    # Read in the data
    cleaned_genes = import_filtered_genes(args.cleaned_genes, logger)

    # Sort the genes by start position
    cleaned_genes = cleaned_genes.sort_values(["Chromosome", "Start"])

    # For each gene, calculate the distance to the next gene, whichever is
    # closer of the upstream or downstream.
    cleaned_genes["distance_to_next_gene"] = (
        cleaned_genes.groupby("Chromosome")["Start"].shift(-1) - cleaned_genes["Stop"]
    )

    # Set anything that is negative (overlapping genes) to 0
    cleaned_genes.loc[
        cleaned_genes["distance_to_next_gene"] < 0, "distance_to_next_gene"
    ] = 0

    # Fill NaN values with the previous non-NaN value within each group
    cleaned_genes["distance_to_next_gene"] = cleaned_genes.groupby("Chromosome")[
        "distance_to_next_gene"
    ].fillna(method="ffill")

    # Scale the distances to KB
    cleaned_genes["distance_to_next_gene"] = (
        cleaned_genes["distance_to_next_gene"] / 1000
    )

    mean_distance = cleaned_genes["distance_to_next_gene"].mean()
    standard_deviation = cleaned_genes["distance_to_next_gene"].std(ddof=0)
    max_distance = cleaned_genes["distance_to_next_gene"].max()

    # Show me the gene and the genes around it that have the max distance
    # print(cleaned_genes.loc[cleaned_genes["distance_to_next_gene"] == max_distance])
    mean_distance_plotobj = plt.Line2D(
        [],
        [],
        color="r",
        marker="",
        linestyle="",
        label=f"Mean distance between all genes: {mean_distance:.3}",
    )
    standard_deviation_plotobj = plt.Line2D(
        [],
        [],
        color="r",
        marker="",
        linestyle="",
        label=f"Standard deviation of distance between all genes: {standard_deviation:.3}",
    )

    # We have already calculated the mean and standard deviation, let's trim
    # the dataset to only include distances less than 250KB, just so that it
    # plots more consistently with the bins. The outliers in the RR dataset
    # make the bins chunky. There are like 4 genes above 250 KB distane, one is
    # at 827KB...
    cleaned_genes = cleaned_genes.loc[cleaned_genes["distance_to_next_gene"] < 250]

    # Create a bar graph of the distances between genes
    values = cleaned_genes["distance_to_next_gene"].to_list()
    plt.hist(values, bins=100)
    plt.ylim(1, 100000)
    plt.xlim(0, 250)
    plt.yscale("log")
    plt.xlabel("Distance Between Genes (KB)")
    plt.ylabel("Number of Genes")
    plt.xlim(0, 175)
    plt.title(f"{args.genome_name}: Distribution of Distances Between Genes")
    plt.legend(
        handles=[mean_distance_plotobj, standard_deviation_plotobj],
        loc="upper right",
    )
    file_out = os.path.join(args.output_dir, args.genome_name + "_gene_distances.png")
    logger.info(f"Saving plot to: {file_out}")
    plt.tight_layout()
    plt.savefig(file_out)
