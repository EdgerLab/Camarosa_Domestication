#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs


"""
- The Del Norte protein file that was used for BLAST has gene names
    that are not the same as the gene names in the final gene annotation
- This script replaces the gene names in the Del Norte protein file.
- This step ought to performed before filtering the results for identifying all
    syntelogs and homologs
- This is performed with the BLAST output and a file that maps the old gene, a
    decoder ring
"""


def import_unclean_homologs(homolog_input_file):
    """
    Import the homologs from the raw file
    """
    col_names = [
        "Query",
        "Subject",
        "Percent ID",
        "Alignment Length",
        "Mismatches",
        "Gap Openings",
        "Q_Start",
        "Q_Stop",
        "Subject_Start",
        "Subject_Stop",
        "E_Value",
        "Bit_Score",
    ]

    col_to_use = [
        "Query",
        "Subject",
        "E_Value",
    ]

    homolog_pd = pd.read_csv(
        homolog_input_file,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
    )

    # Make sure E-Value is float64
    homolog_pd.E_Value = homolog_pd.E_Value.astype("float64")
    # Trim E-values less than 0.05
    homolog_pd = homolog_pd.loc[homolog_pd["E_Value"] < 0.05]

    # Sort by Name and E-Values
    homolog_pd.sort_values(
        by=["Query", "E_Value"], ascending=(True, True), inplace=True
    )

    # Need to take first occurrence of a duplicated gene, the one with the smallest
    # E-Value
    homolog_pd = homolog_pd.drop_duplicates(subset=["Query"], keep="first")

    # Rename columns
    homolog_pd.rename(
        columns={"Query": "Del_Norte", "Subject": "Royal_Royce"}, inplace=True
    )

    # Get the correct name for the Del_Norte genes
    homolog_pd["Del_Norte"] = homolog_pd["Del_Norte"].str.split("-mRNA-1").str[0]

    # Get the correct name for the Royal_Royce genes
    homolog_pd["Royal_Royce"] = homolog_pd["Royal_Royce"].str.split(".").str[0]

    # Add column with identifier so we can later see what source we derived the
    # gene pair from
    homolog_pd["Point_of_Origin"] = "BLAST"

    return homolog_pd


def import_decoder_ring(decoder_ring_input_file):
    decoder_ring = pd.read_csv(
        decoder_ring_input_file, sep="\t", usecols=[3, 4], names=["Old", "New"]
    )
    decoder_ring = dict(zip(decoder_ring["Old"], decoder_ring["New"]))
    return decoder_ring


def replace_names(homolog_pd, decoder_ring):
    homolog_pd["Del_Norte"] = homolog_pd["Del_Norte"].replace(decoder_ring)
    return homolog_pd


if __name__ == "__main__":

    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "homolog_input_file",
        type=str,
        help="parent path of homolog file",
    )
    parser.add_argument(
        "decoder_ring_input_file",
        type=str,
        help="parent path of decoder ring file",
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
    args.homolog_input_file = os.path.abspath(args.homolog_input_file)
    args.decoder_ring_input_file = os.path.abspath(args.decoder_ring_input_file)
    args.output_file = os.path.abspath(args.output_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    logger.info(f"Importing decoder ring...")
    decoder_ring = import_decoder_ring(args.decoder_ring_input_file)

    logger.info(f"Importing unclean homologs...")
    unclean_homologs = import_unclean_homologs(args.homolog_input_file)

    logger.info(f"Replacing names...")
    name_replaced_homologs = replace_names(unclean_homologs, decoder_ring)

    logger.info(f"Saving results to disk at: {args.output_file}")
    name_replaced_homologs.to_csv(args.output_file, sep="\t", header=True, index=False)
