__author__ = "Scott Teresi"

"""
TODO
"""

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd

from transposon.import_filtered_genes import import_filtered_genes
from src.go_analysis.upset_plot import read_go_enrichment_table
from src.intersect_sweeps_w_enriched_terms import (
    read_domestication_sweep_table,
    intersect_genes_with_sweep_zones,
)

if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "cleaned_genes",
        type=str,
        help="Parent path of the cleaned gene annotation file",
    )
    parser.add_argument("domestication_sweep_table", type=str, help="TODO")
    parser.add_argument(
        "preprocessed_go_enrichment_table",
        type=str,
        help="TODO",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.cleaned_genes = os.path.abspath(args.cleaned_genes)
    args.domestication_sweep_table = os.path.abspath(args.domestication_sweep_table)
    args.preprocessed_go_enrichment_table = os.path.abspath(
        args.preprocessed_go_enrichment_table
    )
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -------------------------------------------------------------- #
    # Analysis
    cleaned_genes = import_filtered_genes(args.cleaned_genes, logger)
    cleaned_genes = cleaned_genes.sort_values(["Chromosome", "Start"])
    print(cleaned_genes["Length"].mean())
    length_of_all_genes = cleaned_genes["Length"].sum()

    sweep_table = read_domestication_sweep_table(args.domestication_sweep_table)
    sweep_table = sweep_table.loc[sweep_table["Breeding_Program"] == "UF"]
    sweep_table["Length"] = sweep_table["Stop"] - sweep_table["Start"]
    length_of_all_sweeps = sweep_table["Length"].sum()
    print(length_of_all_genes)
    print(length_of_all_sweeps)

    go_enrichment_table = read_go_enrichment_table(
        args.preprocessed_go_enrichment_table
    )

    print(cleaned_genes)
    cleaned_genes.reset_index(inplace=True)
    cleaned_genes.rename(
        columns={
            "Gene_Name": "RR_Gene",
            "Start": "RR_Start",
            "Stop": "RR_Stop",
            "Chromosome": "RR_Chromosome",
        },
        inplace=True,
    )
    intersected_annotation = intersect_genes_with_sweep_zones(
        sweep_table, cleaned_genes
    )
    print(intersected_annotation)

    intersected_go_enrichment_table = intersect_genes_with_sweep_zones(
        sweep_table, go_enrichment_table
    )
    print(intersected_go_enrichment_table)
    raise ValueError
    intersected_go_enrichment_table.to_csv(
        "intersected_go_enrichment_table.tsv", header=True, index=True, sep="\t"
    )
    number_unique_go_intersected_genes = len(
        intersected_go_enrichment_table.reset_index()["RR_Gene"].unique()
    )
    print(number_unique_go_intersected_genes)

    print(sweep_table)

    x = sweep_table.groupby(by=["Breeding_Program"]).sum()
    print(x)
