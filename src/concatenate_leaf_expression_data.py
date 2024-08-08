#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs
from functools import reduce

"""
- Import unclean gene expression data previously downloaded from Edger Lab
    (RR Strawberry Leaf Match ATACnMNASEseq)
- Filter and combine unclean gene expression data and save to disk
- Provide a function to read the filtered output file from disk
    (read_filtered_expression_data)
"""


def read_filtered_expression_data(input_path_str):
    data = pd.read_csv(
        input_path_str,
        sep="\t",
        header=0,
        names=["Gene_Name", "TPKM_Rep1", "TPKM_Rep2", "TPKM_Rep3", "Avg_Expression"],
        dtype={
            "Gene_Name": str,
            "TPKM_Rep1": np.float64,
            "TPKM_Rep2": np.float64,
            "TPKM_Rep3": np.float64,
            "Avg_Expression": np.float64,
        },
    )
    return data


def import_unclean_expression_RR(expression_input_file, rep_id):
    """
    Import and filter unclean expression data from RR (Royal Royce).

    Args:
        expression_input_file (str): Path to unclean expression data

    Returns:
        TODO edit this!
        expression_data (pandas.DataFrame): Filtered expression data
            Columns:
                Gene (str): Gene name
                RR2_C (int): Expression value for rep 2 (cold)
                RR2_W (int): Expression value for rep 2 (warm)
                RR3_C (int): Expression value for rep 3 (cold)
                RR3_W (int): Expression value for rep 3 (warm)
                RR4_C (int): Expression value for rep 4 (cold)
                RR4_W (int): Expression value for rep 4 (warm)
    """
    # Import unclean gene expression data
    TPKM_colname_w_rep_id = "TPKM_" + rep_id
    data = pd.read_csv(
        expression_input_file,
        sep="\t",
        header=0,
        names=[
            "Gene_Name",
            "Length",
            "EffectiveLength",
            "TPM",
            "NumReads",
            TPKM_colname_w_rep_id,
        ],
        usecols=["Gene_Name", TPKM_colname_w_rep_id],
        dtype={"Gene_Name": str, TPKM_colname_w_rep_id: np.float64},
    )

    # NOTE, need to trim the exon number from the gene names and drop any
    # duplicate genes that get through because of the messy data that I was
    # given, not going to deal with the non-primary transcripts
    data = data.loc[data["Gene_Name"].str.contains(".1", regex=False)]
    data["Gene_Name"] = data["Gene_Name"].str.replace(".1", "", regex=False)
    return data


if __name__ == "__main__":
    # Parse command line arguments
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument("expression_data_rep1", type=str)
    parser.add_argument("expression_data_rep2", type=str)
    parser.add_argument("expression_data_rep3", type=str)
    parser.add_argument(
        "output_file",
        type=str,
        help="Path to output the concatenated expression data to disk",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.expression_data_rep1 = os.path.abspath(args.expression_data_rep1)
    args.expression_data_rep2 = os.path.abspath(args.expression_data_rep2)
    args.expression_data_rep3 = os.path.abspath(args.expression_data_rep3)
    args.output_file = os.path.abspath(args.output_file)

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # ------------------------------------------------------------------------------------------------

    # Import and filter unclean gene expression data from Royal Royce
    rep1 = import_unclean_expression_RR(args.expression_data_rep1, "Rep1")
    rep2 = import_unclean_expression_RR(args.expression_data_rep2, "Rep2")
    rep3 = import_unclean_expression_RR(args.expression_data_rep3, "Rep3")

    merged_expression = reduce(
        lambda left, right: pd.merge(left, right, on=["Gene_Name"], how="outer"),
        [rep1, rep2, rep3],
    )

    # NOTE calculating the average expression ignores NaN values
    merged_expression["Avg_Expression"] = merged_expression[
        ["TPKM_Rep1", "TPKM_Rep2", "TPKM_Rep3"]
    ].mean(axis=1)

    logger.info(f"Writing filtered output data to {args.output_file}")
    merged_expression.to_csv(args.output_file, sep="\t", index=False, header=True)
