#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs

from src.orthologs.syntelogs import read_cleaned_syntelogs
from src.orthologs.homologs import read_cleaned_homologs

"""
- Create the table of orthologs (homologs and syntelogs) for the DN and RR genomes
"""


def merge_homologs_and_syntelogs(homologs, syntelogs):
    """
    Merge the homolog and syntelog dataframes
    """

    # Sort the syntelogs in alphabetical order by Royal Royce gene name
    syntelogs.sort_values(by=["RR_Gene"], inplace=True)

    # Are there any Royal Royce syntelogs that are not unique?
    # Yes indeed there are quite a lot...
    # print(syntelogs.loc[syntelogs.duplicated(subset=["RR_Gene"], keep=False)])

    # Let's merge the homologs and syntelogs
    # First I will rename the columns in the homolog file to make it easier to
    # merge
    homologs.rename(
        columns={
            "Del_Norte": "DN_Gene",
            "Royal_Royce": "RR_Gene",
            "E_Value": "BLAST_E_Value",
        },
        inplace=True,
    )
    # I will also rename the E_Value column in the syntelog file
    syntelogs.rename(columns={"E_Value": "Syntelog_E_Value"}, inplace=True)

    # Now I will merge the two dataframes
    merged_all = pd.concat([homologs, syntelogs], axis=0, join="outer")

    # Sort the merged dataframe by Royal Royce gene name and Point_of_Origin,
    # then by Evalues
    merged_all.sort_values(
        by=["RR_Gene", "Point_of_Origin", "Syntelog_E_Value", "BLAST_E_Value"],
        ascending=[True, False, True, True],
        inplace=True,
    )

    return merged_all


if __name__ == "__main__":

    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument("cleaned_syntelog_input_file", type=str, help="TODO")
    parser.add_argument("cleaned_homolog_input_file", type=str, help="TODO")
    parser.add_argument(
        "output_file",
        type=str,
        help="Path and filename to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.cleaned_syntelog_input_file = os.path.abspath(args.cleaned_syntelog_input_file)
    args.cleaned_homolog_input_file = os.path.abspath(args.cleaned_homolog_input_file)
    args.output_file = os.path.abspath(args.output_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Begin work, read in the data
    syntelogs = read_cleaned_syntelogs(args.cleaned_syntelog_input_file)
    homologs = read_cleaned_homologs(args.cleaned_homolog_input_file)

    orthologs = merge_homologs_and_syntelogs(homologs, syntelogs)

    logger.info(f"Writing orthologs to {args.output_file}")
    orthologs.to_csv(args.output_file, sep="\t", index=False, header=True)
