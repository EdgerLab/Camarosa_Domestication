#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs

"""
- Import unclean syntelog data previously downloaded from CoGe
- Filter unclean syntelog data and save to disk
- Provide helper reader function to read the cleaned data from disk
"""


def import_unclean_syntelogs(syntelog_input_file):
    """
    Import the syntelogs from the raw file and manage data filtration
    """
    col_names = [
        "OrgA_Chromosome",
        "OrgA_Gene_Region",
        "OrgA_Start",
        "OrgA_Stop",
        "OrgB_Chromosome",
        "OrgB_Gene_Region",
        "OrgB_Start",
        "OrgB_Stop",
        "E_Value",
        "Diagonal_Score",
        "Web_Link",
    ]

    col_to_use = [
        "OrgA_Chromosome",
        "OrgA_Gene_Region",
        "OrgB_Chromosome",
        "OrgB_Gene_Region",
        "E_Value",
        "Diagonal_Score",
    ]

    syntelogs = pd.read_csv(
        syntelog_input_file,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
        dtype={
            "OrgA_Chromosome": str,
            "OrgA_Gene_Region": str,
            "OrgB_Chromosome": str,
            "OrgB_Gene_Region": str,
            "E_Value": np.float64,
            "Diagonal_Score": np.int32,
        },
    )

    return syntelogs


def filter_syntelogs(syntelogs):
    """
    Filter the syntelogs to rename the columns so that they make more sense.
    Also, remove syntelogs with E-values greater than 0.05.
    Parse the string names for the genes and chromosomes to remove the
    extraneous information

    Args:
        syntelogs (pd.DataFrame): unclean syntelog data

    Returns:
        syntelogs (pd.DataFrame): clean syntelog data
    """

    # Get the correct name for the genes
    # MAGIC
    syntelogs["OrgA_Gene_Region"] = (
        syntelogs["OrgA_Gene_Region"].str.split("\|\|").str[3]
    )
    # MAGIC
    syntelogs["OrgB_Gene_Region"] = (
        syntelogs["OrgB_Gene_Region"].str.split("\|\|").str[3]
    )

    # Get the correct name for the chromosome
    # MAGIC splits to remove the nonsense info from SynMap
    syntelogs["OrgA_Chromosome"] = (
        syntelogs["OrgA_Chromosome"].str.split("_", n=1).str[1]
    )
    syntelogs["OrgB_Chromosome"] = (
        syntelogs["OrgB_Chromosome"].str.split("_", n=1).str[1]
    )

    # Fix Royal Royce genes in A
    syntelogs["OrgA_Gene_Region"] = syntelogs["OrgA_Gene_Region"].str.split("\.").str[0]

    # Fix Del Norte genes in B
    syntelogs["OrgB_Gene_Region"] = (
        syntelogs["OrgB_Gene_Region"].str.split("-mRNA-1").str[0]
    )

    # Rename columns so that they make more sense
    syntelogs.rename(
        columns={
            "OrgA_Gene_Region": "RR_Gene",
            "OrgB_Gene_Region": "DN_Gene",
            "OrgA_Chromosome": "RR_Chromosome",
            "OrgB_Chromosome": "DN_Chromosome",
        },
        inplace=True,
    )

    # MAGIC trim E-values less than 0.05
    syntelogs = syntelogs.loc[syntelogs["E_Value"] < 0.05]

    # Add column with identifier so we can later see what source we derived the
    # gene pair from
    syntelogs["Point_of_Origin"] = "Synteny"

    syntelogs.drop(columns=["Diagonal_Score"], inplace=True)

    return syntelogs


def save_clean_syntelogs(syntelogs, output_path, logger):
    """
    Save the clean syntelog data to disk

    Args:
        syntelogs (pd.DataFrame): clean syntelog data
        output_path (str): path to save the clean syntelog data
        logger (logging.Logger): logger object
    """
    logger.info(f"Saving clean syntelog file to disk at: {output_path}")
    syntelogs.to_csv(output_path, header=True, sep="\t", index=False)


def read_cleaned_syntelogs(syntelog_input_file):
    """
    Import the clean syntelogs from the pre-filtered file
    """
    return pd.read_csv(syntelog_input_file, sep="\t", header="infer")


if __name__ == "__main__":

    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "syntelog_input_file",
        type=str,
        help="parent path of syntelog file",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path and filename to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.syntelog_input_file = os.path.abspath(args.syntelog_input_file)
    args.output_file = os.path.abspath(args.output_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    unclean_syntelogs = import_unclean_syntelogs(args.syntelog_input_file)
    clean_syntelogs = filter_syntelogs(unclean_syntelogs)

    save_clean_syntelogs(clean_syntelogs, args.output_file, logger)
