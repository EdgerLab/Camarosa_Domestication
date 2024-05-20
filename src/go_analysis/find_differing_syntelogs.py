#!/usr/bin/env/python

"""
TODO
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
from src.extract_AED_score import read_aed_output_table
from src.go_analysis.find_abnormal_genes import (
    perform_upper_cutoff_and_subset,
    perform_lower_cutoff_and_subset,
    subset_by_arabidopsis_presence,
    save_table_to_disk,
    subset_by_aed_score,
)

if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)

    parser = argparse.ArgumentParser(description="TODO")
    parser.add_argument(
        "preprocessed_density_table",
        type=str,
        help="TODO",
    )
    parser.add_argument("dn_aed_score_table", type=str)
    parser.add_argument("rr_aed_score_table", type=str)
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
    parser.add_argument("strawberry_ortholog_table", type=str)
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
    args.strawberry_ortholog_table = os.path.abspath(args.strawberry_ortholog_table)
    args.dn_aed_score_table = os.path.abspath(args.dn_aed_score_table)
    args.rr_aed_score_table = os.path.abspath(args.rr_aed_score_table)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -------------------------------------------------------------------------

    # Load in the pre-filtered TE Density data
    base_table = pd.read_csv(args.preprocessed_density_table, header="infer", sep="\t")

    # print(base_table)
    # print(base_table.columns)
    # raise ValueError

    # Load in the pre-filtered ortholog information
    ortholog_table = pd.read_csv(
        args.strawberry_ortholog_table, header="infer", sep="\t", low_memory=False
    )

    te_type, window, direction = decode_te_window_direction_str(
        os.path.basename(args.preprocessed_density_table)
    )

    dn_aed_scores = read_aed_output_table(args.dn_aed_score_table)
    rr_aed_scores = read_aed_output_table(args.rr_aed_score_table)
    for aed_scores, genome in zip([dn_aed_scores, rr_aed_scores], ["DN", "RR"]):
        aed_scores.drop(columns=["Chromosome"], inplace=True)
        aed_scores.rename(columns={"Gene_Name": f"{genome}_Gene"}, inplace=True)
    # MAGIC MAGIC MAGIC, anything below this is a well-supported gene model
    good_aed_score = 0.75

    # Create a named tuple that is the cutoff function and its percentile in
    # integer format, to avoid needing to retype this all the time
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

        # Apply an AED score cutoff to the dense gene table
        for aed_scores, genome in zip([dn_aed_scores, rr_aed_scores], ["DN", "RR"]):
            mod_table = subset_by_aed_score(
                mod_table, aed_scores, genome, good_aed_score
            )
            # Quick and dirty way to avoid trying to add in the column twice,
            # we don't need it after the subset, would fail on second iteration
            # otherwise
            mod_table.drop(columns=["AED_Score"], inplace=True)

        # Subset the table so that each entry MUST have an Arabidopsis gene
        mod_table = subset_by_arabidopsis_presence(mod_table)
        out_filename = f"{col_to_save}_Biased_Towards_{genome}_{str(i.percentile)}_density_percentile.tsv"
        save_table_to_disk(
            mod_table,
            args.output_dir,
            out_filename,
            logger,
        )
