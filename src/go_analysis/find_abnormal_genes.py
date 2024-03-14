#!/usr/bin/env/python

"""
Given a table of pre-filtered TE Density of strawberry genes, this script
generates more tables by calculating the top X and bottom Y percentiles of data
and saving them to disk. The plan is to use these tables to run TopGO and
determine which genes are over or under represented in the datasets.
This script also calculates the top and bottom percentiles of the 'Difference'
column which is was previously calculated and is the Royal Royce (RR) TE Density value
substracted from the Del Norte (DN) TE Density value.

An example of an output table is the top 1% most LTR-dense genes in the RR
genome for the 5KB upstream window.
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd
from collections import namedtuple

from src.syntelog_differences.bargraphs import decode_te_window_direction_str


def calculate_cutoff_value(arr, cutoff_int):
    return np.percentile(arr, cutoff_int)


def perform_upper_cutoff_and_subset(
    table, column_name, upper_percentile_cutoff_int, logger
):
    """
    TODO
    """
    table = table.copy(deep=True)

    # Make sure nonzero values are used
    # TODO verify this is the right thing to do, do we want to remove non-zeros
    # before we calculate the cutoff?
    # But most genes have a 0 value for TE density, so it makes sense to test
    # those who do have a TE near them
    table = table.loc[table[column_name] != 0.0]

    upper_cutoff_val = calculate_cutoff_value(
        table[column_name], upper_percentile_cutoff_int
    )

    # Check the cutoff value and raise a warning if the value is below 10%
    if upper_cutoff_val < 0.1:
        logger.warning(
            f"Upper Cutoff Value is below 10% density: it is {upper_cutoff_val} for {column_name}"
        )

    # This column is used to store the cutoff value for each gene and may
    # not actually be needed
    table["Upper_Cutoff_Value"] = upper_cutoff_val

    # Subset the table by the cutoff value
    table = table.loc[table[column_name] >= upper_cutoff_val]
    return table


def perform_lower_cutoff_and_subset(
    table, column_name, lower_percentile_cutoff_int, logger
):
    """
    This is only needed for the difference column?
    """
    table = table.copy(deep=True)

    # Make sure nonzero values are used
    table = table.loc[table[column_name] != 0.0]

    lower_cutoff_val = calculate_cutoff_value(
        table[column_name], lower_percentile_cutoff_int
    )

    # Check the cutoff value and raise a warning if the value is below 10%
    # TODO revist if this warning is needed for the lower cutoff
    if lower_cutoff_val > -0.1:
        logger.warning(
            f"Lower Cutoff Value is below 10%: {lower_cutoff_val} for {column_name}"
        )

    # This column is used to store the cutoff value for each gene and may
    # not actually be needed
    table["Lower_Cutoff_Value"] = lower_cutoff_val

    # Subset the table by the cutoff value
    # NOTE this is like the only thing different from the upper cutoff
    # function... TODO think about how to refactor
    table = table.loc[table[column_name] <= lower_cutoff_val]
    return table


def subset_by_arabidopsis_presence(table):
    """
    Return a table where each entry has an Arabidopsis gene
    """
    return table.loc[~table["Arabidopsis_Gene"].isna()]


def save_table_to_disk(table, output_dir, out_filename, logger):
    out_filepath = os.path.join(output_dir, out_filename)
    logger.info(f"Writing to: {out_filepath}")
    table.to_csv(out_filepath, sep="\t", header=True, index=False)


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "preprocessed_density_table",
        type=str,
        help="TODO",
    )

    parser.add_argument(
        "upper_percentile_cutoff_int",
        type=int,
        help="TODO",
    )
    parser.add_argument(
        "lower_percentile_cutoff_int",
        type=int,
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
    args.preprocessed_density_table = os.path.abspath(args.preprocessed_density_table)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -------------------------------------------------------------------------

    # Load in the pre-filtered TE Density data
    base_table = pd.read_csv(args.preprocessed_density_table, header="infer", sep="\t")
    te_type, window, direction = decode_te_window_direction_str(
        os.path.basename(args.preprocessed_density_table)
    )

    # Create a named tuple that is the cutoff function and its percentile in
    # integer format
    cutoff_function_w_percentile = namedtuple(
        "cutoff_function_w_percentile",
        ["cutoff_function", "percentile", "upper_or_lower_str"],
    )
    upper = cutoff_function_w_percentile(
        perform_upper_cutoff_and_subset, args.upper_percentile_cutoff_int, "Upper"
    )
    lower = cutoff_function_w_percentile(
        perform_lower_cutoff_and_subset, args.lower_percentile_cutoff_int, "Lower"
    )

    # This calculates the top X% of the genes for each genome, and then the
    # top and bottom for the difference.
    # MAGIC genome names
    for genome in ["H4", "DN", "RR"]:
        te_col = f"{genome}_{te_type}_{window}_{direction}"

        # Use the functions as first class objects to iterate over them, and
        # calculate the upper and lower cutoffs
        # TODO consider making these a named tuple
        for i in [upper]:
            if genome == "H4":
                # TODO it is possible for H4 to have duplicate gene names before we run
                # TopGO.
                # H4 has NaN so we need to remove those before calculations
                mod_table = base_table.copy(deep=True)
                mod_table.dropna(axis=0, subset=["H4_Gene"], inplace=True)
                mod_table = i.cutoff_function(mod_table, te_col, i.percentile, logger)
            else:
                # Perform the cutoffs and subset
                mod_table = i.cutoff_function(base_table, te_col, i.percentile, logger)

            # Subset the table so that each entry MUST have an Arabidopsis gene
            mod_table = subset_by_arabidopsis_presence(mod_table)
            out_filename = f"{te_col}_{i.upper_or_lower_str}_{str(i.percentile)}_density_percentile.tsv"
            save_table_to_disk(
                mod_table,
                args.output_dir,
                out_filename,
                logger,
            )

    # ----------------------------------------
    # Calculate the 'Difference' column
    column_name = "Difference"
    col_to_save = f"{column_name}_{te_type}_{window}_{direction}"
    for genome, i in [("DN", upper), ("RR", lower)]:
        # We don't have a 'Difference' column for H4, the difference
        # column is for DN vs RR.
        # NOTE MAGIC, if positive it is a DEL NORTE biased gene pair
        # NOTE MAGIC, if negative it is a Royal Royce biased gene pair, hence
        # we will apply the lower cutoff function to the RR dataset

        # Perform the cutoffs and subset
        mod_table = i.cutoff_function(base_table, column_name, i.percentile, logger)

        mod_table = subset_by_arabidopsis_presence(mod_table)
        out_filename = f"{col_to_save}_BiasedTowards_{genome}_{str(i.percentile)}_density_percentile.tsv"
        save_table_to_disk(
            mod_table,
            args.output_dir,
            out_filename,
            logger,
        )
