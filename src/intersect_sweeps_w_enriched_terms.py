__author__ = "Scott Teresi"

"""
Intersect the domestication sweeps with the enriched (RR) GO terms to see if any of the sweeps overlap with the enriched terms.

I have two paradigms for this:
    1. The syntelog RR biased files
    2. The super-dense RR files that don't have a DN gene.
"""

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd

from src.go_analysis.upset_plot import read_go_enrichment_table


def read_domestication_sweep_table(filename):
    data = pd.read_csv(
        filename,
        sep=",",
        names=[
            "Sweep_ID",
            "Chromosome_ID",
            "Number_of_Sig_Windows_in_Sweep",
            "Start",
            "Stop",
            "Max_XPLCR_Val",
            "Window_Start_For_Max_XPLCR_Val",
            "Window_Stop_For_Max_XPLCR_Val",
            "Size_of_Sweep",
            "Selection_Coefficient",
            "Breeding_Program",
        ],
        header=0,
        dtype={
            "Sweep_ID": str,
            "Chromosome_ID": str,
            "Number_of_Sig_Windows_in_Sweep": int,
            "Start": int,
            "Stop": int,
            "Max_XPLCR_Val": float,
            "Window_Start_For_Max_XPLCR_Val": int,
            "Window_Stop_For_Max_XPLCR_Val": int,
            "Size_of_Sweep": int,
            "Selection_Coefficient": float,
            "Breeding_Program": str,
        },
    )

    data = add_useful_IDs_to_sweep_table(data)
    return data


def add_useful_IDs_to_sweep_table(sweep_table):
    """
    The sweep table from the publication has integer sweep IDs that are not
    unique, so the only distinguising feature within the table are the start
    and stop positions, which is not ideal. This function adds a new column
    that is the Sweep ID and a greek letter name, it renames the original
    Sweep_ID columns to OLD_Sweep_ID
    """
    delim = "_"
    sweep_table["New_Sweep_ID"] = (
        sweep_table["Breeding_Program"]
        + delim
        + sweep_table["Chromosome_ID"]
        + delim
        + sweep_table["Sweep_ID"].astype(str)
    )
    sweep_table.set_index("New_Sweep_ID", inplace=True)

    # Remove 'chr_' string from the front end of the chromosome ID so that it
    # matches with my data
    sweep_table["Chromosome_ID"] = sweep_table["Chromosome_ID"].str.replace(
        "chr_", "", regex=True
    )

    sweep_table = reorder_sweep_table_columns(sweep_table)
    return sweep_table


def reorder_sweep_table_columns(sweep_table):
    columns_in_desired_order = [
        "Selection_Coefficient",
        "Start",
        "Stop",
        "Max_XPLCR_Val",
        "Window_Start_For_Max_XPLCR_Val",
        "Window_Stop_For_Max_XPLCR_Val",
        "Size_of_Sweep",
        "Number_of_Sig_Windows_in_Sweep",
        "Chromosome_ID",
        "Breeding_Program",
        "Sweep_ID",
    ]
    sweep_table = sweep_table[columns_in_desired_order]
    return sweep_table


def intersect_genes_with_sweep_zones(sweep_table, go_enrichment_table):
    """
    Subset the go_enrichment_table table so that only the genes that are within the
    sweep zones remain. This operates on the start and stop positions of the RR
    genes, and special attention must be paid to make sure the chromosomes
    match
    """
    go_enrichment_table.set_index(
        ["RR_Gene", "Arabidopsis_Gene", "GO_ID", "Term"], inplace=True
    )

    for sweep_id in sweep_table.index:
        go_enrichment_table[sweep_id] = "N"

    go_enrichment_table = go_enrichment_table.copy(deep=True)

    for sweep_id, sweep_row in sweep_table.iterrows():
        sweep_start = sweep_row["Start"]
        sweep_stop = sweep_row["Stop"]
        sweep_chromosome = sweep_row["Chromosome_ID"]

        go_enrichment_table.loc[
            (go_enrichment_table["RR_Start"] >= sweep_start)
            & (go_enrichment_table["RR_Stop"] <= sweep_stop)
            & (go_enrichment_table["RR_Chromosome"] == sweep_chromosome),
            sweep_id,
        ] = "Y"

    # Subset the table further so that only rows with at least one 'Y' remain
    mask = go_enrichment_table[sweep_table.index] == "Y"
    go_enrichment_table = go_enrichment_table.loc[mask.any(axis=1)]
    return go_enrichment_table


def save_sweep_table_with_enriched_terms(mod_table, output_filename):
    """
    Save the sweep table with the enriched terms as a new table
    """
    mod_table.to_csv(output_filename, header=True, index=True, sep="\t")
    return None


def describe_gene_of_interest(mod_table, sweep_table_cols, gene_of_interest):
    """
    Describe the gene of interest by printing out the RR start and stop
    positions, the total TE density 5000 bp upstream of the gene, and the
    enriched terms
    """
    # subset = mod_table.loc[gene_of_interest, :]  # OLD, unused?

    subset = mod_table.reset_index()
    subset = subset.loc[subset["RR_Gene"] == gene_of_interest]

    # TODO this may change when we have the DN syntelog code
    subset = subset.set_index(["RR_Gene", "Arabidopsis_Gene", "GO_ID", "Term"])
    bool_df = subset.loc[:, sweep_table_cols] == "Y"
    any_masked = bool_df.any(axis="index")
    filtered = subset[any_masked.index[any_masked]]

    return filtered


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "preprocessed_go_enrichment_table",
        type=str,
        help="TODO",
    )
    parser.add_argument(
        "intersected_output",
        type=str,
        help="TODO",
    )
    parser.add_argument("domestication_sweep_table", type=str, help="TODO")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.preprocessed_go_enrichment_table = os.path.abspath(
        args.preprocessed_go_enrichment_table
    )
    args.domestication_sweep_table = os.path.abspath(args.domestication_sweep_table)
    args.intersected_output = os.path.abspath(args.intersected_output)
    go_enrichment_table = read_go_enrichment_table(
        args.preprocessed_go_enrichment_table
    )
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    sweep_table = read_domestication_sweep_table(args.domestication_sweep_table)
    mod_table = intersect_genes_with_sweep_zones(sweep_table, go_enrichment_table)

    # NOTE MAGIC, if the user wants to print their gene of interest
    gene_of_interest = "Fxa1Cg101834"
    result = describe_gene_of_interest(
        mod_table, sweep_table.index.to_list(), gene_of_interest
    )
    # print(result)

    # Save the master table
    logger.info(f"Saving the modified sweep table to {args.intersected_output}")
    save_sweep_table_with_enriched_terms(mod_table, args.intersected_output)
