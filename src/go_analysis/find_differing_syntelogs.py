#!/usr/bin/env/python

"""
Generate a table of syntelogs that have a large TE density bias towards either
DN or RR. This script will take in a preprocessed density table and a table of
AED scores. Perform a secondary cutoff of AED scores to make sure the syntelogs
are well supported.

Unlike `find_abnormal_genes.py` this script does not calculate the top X% of
TE-dense genes, it merely identifies the genes that have an absolute difference
of 0.75 or greater in TE density between DN and RR. This value is hard-coded as
`good_aed_score`
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


class SyntelogDifferences:
    def __init__(self, dn_aed_scores, rr_aed_scores, te_table):
        self.dn_aed_scores = dn_aed_scores
        self.rr_aed_scores = rr_aed_scores
        self.te_table = te_table

        self.biased_towards_dn = None  # initialize to None, b/c set later
        self.biased_towards_rr = None  # initialize to None, b/c set later

    def subset_by_te_difference(self, cutoff=0.75):
        """
        Subset the table so that we only have syntelogs with extreme TE
        difference
        """
        self.biased_towards_dn = self.te_table.loc[
            self.te_table["Difference"] >= cutoff
        ]
        self.biased_towards_rr = self.te_table.loc[
            self.te_table["Difference"] <= -cutoff
        ]

    def subset_by_arabidopsis(self):
        """
        Subset the table so that each entry MUST have an Arabidopsis gene
        """
        self.biased_towards_dn = subset_by_arabidopsis_presence(self.biased_towards_dn)
        self.biased_towards_rr = subset_by_arabidopsis_presence(self.biased_towards_rr)

    def reorder_columns(self):
        """
        Reorder the columns so that the 'Difference' column is towards the front
        """
        data.biased_towards_dn.insert(
            2, "Difference", data.biased_towards_dn.pop("Difference")
        )
        data.biased_towards_rr.insert(
            2, "Difference", data.biased_towards_rr.pop("Difference")
        )


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "preprocessed_density_table",
        type=str,
        help="TODO",
    )
    parser.add_argument("dn_aed_score_table", type=str)
    parser.add_argument("rr_aed_score_table", type=str)
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
    args.dn_aed_score_table = os.path.abspath(args.dn_aed_score_table)
    args.rr_aed_score_table = os.path.abspath(args.rr_aed_score_table)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -------------------------------------------------------------------------

    # Load in the pre-filtered TE Density data
    base_table = pd.read_csv(args.preprocessed_density_table, header="infer", sep="\t")

    # MAGIC
    # Subset the base table so that we only have gene pairs between RR and DN
    # that are syntenic
    base_table = base_table.loc[base_table["DN_RR_Point_of_Origin"] == "Synteny"]

    te_type, window, direction = decode_te_window_direction_str(
        os.path.basename(args.preprocessed_density_table)
    )

    dn_aed_scores = read_aed_output_table(args.dn_aed_score_table)
    rr_aed_scores = read_aed_output_table(args.rr_aed_score_table)
    for table, name in zip([dn_aed_scores, rr_aed_scores], ["DN", "RR"]):
        table.drop(columns=["Chromosome"], inplace=True)
        table.rename(columns={"Gene_Name": f"{name}_Gene"}, inplace=True)

    data = SyntelogDifferences(dn_aed_scores, rr_aed_scores, base_table)

    # MAGIC MAGIC MAGIC, anything below this is a well-supported gene model
    good_aed_score = 0.75

    # ----------------------------------------
    # Calculate the 'Difference' column
    column_name = "Difference"
    col_to_save = f"{column_name}_{te_type}_{window}_{direction}"

    # NOTE MAGIC, if positive it is a DEL NORTE biased gene pair
    # NOTE MAGIC, if negative it is a Royal Royce biased gene pair, hence
    # we will apply the lower cutoff function to the RR dataset
    # Perform the cutoffs and subset, MAGIC
    data.subset_by_te_difference()

    # Apply an AED score cutoff to the dense gene table
    # Repeating myself here, bad code.
    data.biased_towards_dn = subset_by_aed_score(
        data.biased_towards_dn, data.dn_aed_scores, "DN", good_aed_score
    )
    data.biased_towards_rr = subset_by_aed_score(
        data.biased_towards_rr, data.rr_aed_scores, "RR", good_aed_score
    )

    # Subset the table so that each entry MUST have an Arabidopsis gene
    data.subset_by_arabidopsis()

    # Reorder the columns so that the 'Difference' column is
    # towards the front
    data.reorder_columns()

    # Save the DN biased table to disk
    dn_out_filename = f"{col_to_save}_Syntelogs_Biased_Towards_DN.tsv"
    rr_out_filename = f"{col_to_save}_Syntelogs_Biased_Towards_RR.tsv"

    for filename, table in zip(
        [dn_out_filename, rr_out_filename],
        [data.biased_towards_dn, data.biased_towards_rr],
    ):
        save_table_to_disk(
            table,
            args.output_dir,
            filename,
            logger,
        )
