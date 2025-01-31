#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs

from transposon.import_filtered_genes import import_filtered_genes


"""
Reformat the BLAST results for the H4 and Royal_Royce genes.
Perform an E Value filter and remove any genes that are not in the gene data
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

    # NOTE take first occurrence of a duplicated gene, the one with the smallest
    # E-Value
    homolog_pd = homolog_pd.drop_duplicates(subset=["Query"], keep="first")

    # Rename columns
    homolog_pd.rename(columns={"Query": "H4", "Subject": "Royal_Royce"}, inplace=True)

    # Get the correct name for the H4 genes
    homolog_pd["H4"] = homolog_pd["H4"].str.split(".").str[0]

    # Get the correct name for the Royal_Royce genes
    homolog_pd["Royal_Royce"] = homolog_pd["Royal_Royce"].str.split(".").str[0]

    # Add column with identifier so we can later see what source we derived the
    # gene pair from
    homolog_pd["Point_of_Origin"] = "BLAST"

    return homolog_pd


if __name__ == "__main__":

    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "homolog_input_file",
        type=str,
        help="parent path of homolog file",
    )

    parser.add_argument(
        "RR_gene_data",
        type=str,
        help="gene data file that is the input to TE Density for RR",
    )

    parser.add_argument(
        "H4_gene_data",
        type=str,
        help="gene data file that is the input to TE Density for H4",
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
    args.output_file = os.path.abspath(args.output_file)

    # Gene Data files for RR and H4, to be used in removing genes that don't
    # match my gene annotation, so we can have an easy merge with the TE
    # Density results later
    args.RR_gene_data = os.path.abspath(args.RR_gene_data)
    args.H4_gene_data = os.path.abspath(args.H4_gene_data)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    logger.info(f"Importing unclean homologs...")
    clean_homologs = import_unclean_homologs(args.homolog_input_file)

    # Rename the columns
    clean_homologs.rename(
        columns={
            "H4": "H4_Gene",
            "Royal_Royce": "RR_Gene",
            "E_Value": "BLAST_E_Value",
        },
        inplace=True,
    )

    # Remove a gene from the homolog table if is not in the cleaned gene file
    cleaned_RR_genes = import_filtered_genes(args.RR_gene_data, logger)
    cleaned_H4_genes = import_filtered_genes(args.H4_gene_data, logger)
    for i in zip(("RR_Gene", "H4_Gene"), (cleaned_RR_genes, cleaned_H4_genes)):
        clean_homologs = clean_homologs.loc[clean_homologs[i[0]].isin(i[1].index)]

    # First merge the data in and verify that they have the same chromosome
    for i in zip((cleaned_RR_genes, cleaned_H4_genes), ("RR_Gene", "H4_Gene")):
        genes = i[0].reset_index()
        genes.rename(
            columns={
                "Gene_Name": i[1],
                "Chromosome": i[1].split("_")[0] + "_Chromosome",
            },
            inplace=True,
        )
        clean_homologs = clean_homologs.merge(genes, on=i[1], how="left")
        clean_homologs.drop(
            columns=["Feature", "Start", "Stop", "Strand", "Length"], inplace=True
        )

    # Perform a chromosome check similar to the syntelog data, don't want BLAST
    # results from different chromosomes
    clean_homologs = clean_homologs.loc[
        clean_homologs["H4_Chromosome"] == clean_homologs["RR_Chromosome"].str[0]
    ]

    logger.info(f"Saving results to disk at: {args.output_file}")
    clean_homologs.to_csv(args.output_file, sep="\t", header=True, index=False)
