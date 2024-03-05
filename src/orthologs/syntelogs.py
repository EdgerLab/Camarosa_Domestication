#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs

from src.orthologs.utils import (
    map_chromosomes,
    reformat_chromosomes_from_SynMap,
    reformat_gene_names_from_SynMap,
    reformat_gene_names_with_period,
    drop_rows_with_bad_val_in_col,
    remove_str_from_val,
)

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
        },
    )

    return syntelogs


# maybe make a custom one for each genome, since the rules are different for
# the H4/DN/RR gene names, and the outputs have different column names
def filter_syntelogs_DN(syntelogs):
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
    # MAGIC splits to remove the nonsense info from SynMap
    for i in ["OrgA_Gene_Region", "OrgB_Gene_Region"]:
        syntelogs = reformat_gene_names_from_SynMap(syntelogs, i)

    # Get the correct name for the chromosome
    # MAGIC splits to remove the nonsense info from SynMap
    for i in ["OrgA_Chromosome", "OrgB_Chromosome"]:
        syntelogs = reformat_chromosomes_from_SynMap(syntelogs, i)

    # Fix Royal Royce genes in A
    syntelogs = reformat_gene_names_with_period(syntelogs, "OrgA_Gene_Region")

    # Fix Del Norte genes in B
    # TODO check with Pat should I eleminate the -mRNA-1? or 2?
    # This is not an issue with H4-RR
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
    syntelogs = remove_str_from_val(syntelogs, "_RagTag", "DN_Chromosome")

    # Drop the rows where the chromosome starts with 'contig'
    for i in ["DN_Chromosome", "RR_Chromosome"]:
        syntelogs = drop_rows_with_bad_val_in_col(syntelogs, "contig", i)

    # Whitelist the genes we want to keep, some chromosomes shouldn't have
    # syntelogs between one another. This is probably biologically real, but
    # for our purposes we aren't trying to concern ourselves with the ancient
    # polyploidy events.
    # MAGIC took this list from the Hardigan paper
    # Renaming the Del Norte chromosomes to match the Royal Royce chromosomes
    syntelogs = map_chromosomes(syntelogs, "DN_Chromosome")

    # Drop the rows where the chromosomes do not match
    syntelogs = syntelogs.loc[syntelogs["DN_Chromosome"] == syntelogs["RR_Chromosome"]]

    # MAGIC trim E-values less than 0.05
    syntelogs = syntelogs.loc[syntelogs["E_Value"] < 0.05]

    # Rename the E-value column to be more descriptive
    syntelogs.rename(columns={"E_Value": "Synteny_E_Value"}, inplace=True)

    # Add column with identifier so we can later see what source we derived the
    # gene pair from
    syntelogs["Point_of_Origin"] = "Synteny"

    return syntelogs


def filter_syntelogs_H4(syntelogs):
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
    # MAGIC splits to remove the nonsense info from SynMap
    for i in ["OrgA_Gene_Region", "OrgB_Gene_Region"]:
        syntelogs = reformat_gene_names_from_SynMap(syntelogs, i)

    # Get the correct name for the chromosome
    # MAGIC splits to remove the nonsense info from SynMap
    for i in ["OrgA_Chromosome", "OrgB_Chromosome"]:
        syntelogs = reformat_chromosomes_from_SynMap(syntelogs, i)

    # Fix H4 and RR genes in A and B
    for i in ["OrgA_Gene_Region", "OrgB_Gene_Region"]:
        syntelogs = reformat_gene_names_with_period(syntelogs, i)

    # Rename columns so that they make more sense
    syntelogs.rename(
        columns={
            "OrgA_Gene_Region": "H4_Gene",
            "OrgB_Gene_Region": "RR_Gene",
            "OrgA_Chromosome": "H4_Chromosome",
            "OrgB_Chromosome": "RR_Chromosome",
        },
        inplace=True,
    )

    # Remove the prefix 'Fvb' from the H4 chromosome names
    syntelogs = remove_str_from_val(syntelogs, "Fvb", "H4_Chromosome")

    # Drop the rows where the chromosome starts with 'contig'
    for i in ["H4_Chromosome", "RR_Chromosome"]:
        syntelogs = drop_rows_with_bad_val_in_col(syntelogs, "contig", i)

    # Whitelist the genes we want to keep, some chromosomes shouldn't have
    # syntelogs between one another. This is probably biologically real, but
    # for our purposes we aren't trying to concern ourselves with the ancient
    # polyploidy events.
    # Drop the rows where the chromosomes do not match
    # TODO CHECK THIS ONE MORE TIME, check with Pat.
    # MAGIC to compare a 7D vs 7
    syntelogs = syntelogs.loc[
        syntelogs["H4_Chromosome"] == syntelogs["RR_Chromosome"].str[0]
    ]

    # MAGIC trim E-values less than 0.05
    syntelogs = syntelogs.loc[syntelogs["E_Value"] < 0.05]

    # Rename the E-value column to be more descriptive
    syntelogs.rename(columns={"E_Value": "Synteny_E_Value"}, inplace=True)

    # Add column with identifier so we can later see what source we derived the
    # gene pair from
    syntelogs["Point_of_Origin"] = "Synteny"

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
    return pd.read_csv(
        syntelog_input_file, sep="\t", header="infer", dtype={"H4_Chromosome": str}
    )


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
        "genome_name",
        type=str,
        choices=["DN", "H4"],
        help="name of genome to be used in proper function call",
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

    if args.genome_name == "DN":
        clean_syntelogs = filter_syntelogs_DN(unclean_syntelogs)
    if args.genome_name == "H4":
        clean_syntelogs = filter_syntelogs_H4(unclean_syntelogs)

    save_clean_syntelogs(clean_syntelogs, args.output_file, logger)
