#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs

"""
- Import unclean gene expression data previously downloaded from Edger Lab
- Filter unclean gene expression data and save to disk
- Provide helper reader function to read the cleaned data from disk
"""


def import_unclean_expression_RR(expression_input_file):
    """
    Import and filter unclean expression data from RR (Royal Royce).

    Args:
        expression_input_file (str): Path to unclean expression data

    Returns:
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
    # TODO go back and ask Jordan why the numbers start with 2, and not 1.

    # Import unclean gene expression data
    expression_data = pd.read_csv(expression_input_file, sep=",", header="infer")
    expression_data.rename(columns={"gene_id": "Gene"}, inplace=True)

    # Calculate the mean of the cold treatments, and the mean of the warm treatments
    # Takes a mean for each gene across the three replicates
    expression_data["RR_Avg_Cold"] = expression_data[["RR2_C", "RR3_C", "RR4_C"]].mean(
        axis=1
    )

    expression_data["RR_Avg_Warm"] = expression_data[["RR2_W", "RR3_W", "RR4_W"]].mean(
        axis=1
    )

    # Calculate the log2 fold change
    # Takes the log2 of the ratio of the warm mean to the cold mean
    # This will give a divide by zero warning, but that's okay as it replaces
    # with a NaN value
    expression_data["log2_fold_change"] = np.log2(
        expression_data["RR_Avg_Warm"] / expression_data["RR_Avg_Cold"]
    )

    # Replace the np.nan values with 'NA'
    # expression_data.replace(np.nan, "NA", inplace=True)
    return expression_data


def save_clean_expression(expression_data, output_file):
    """
    Save filtered gene expression data to disk
    """
    expression_data.to_csv(output_file, index=False, header=True, sep="\t")


if __name__ == "__main__":
    # Parse command line arguments
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(
        description="""Import unclean gene expression data previously
        downloaded from Edger Lab"""
    )

    parser.add_argument(
        "expression_input_file", type=str, help="Path to unclean gene expression data"
    )
    parser.add_argument(
        "output_path",
        type=str,
        help="Path to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.expression_input_file = os.path.abspath(args.expression_input_file)
    args.output_path = os.path.abspath(args.output_path)

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Import and filter unclean gene expression data from Royal Royce
    unclean_expression_RR = import_unclean_expression_RR(args.expression_input_file)

    # NOTE MAGIC filename here
    # Save the filtered gene expression data to disk
    save_clean_expression(
        unclean_expression_RR,
        os.path.join(args.output_path, "Cleaned_Expression_Data_RR.tsv"),
    )

    # ------------------------------------------------------------
    # TODO import and filter unclean gene expression data from Del Norte
    # TODO change this arg.
    # unclean_expression_DN = import_unclean_expression_DN(args.expression_input_file)
    # save_clean_expression(unclean_expression_DN, args.output_path)
