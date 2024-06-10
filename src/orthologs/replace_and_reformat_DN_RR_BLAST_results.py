#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs

from transposon.import_filtered_genes import import_filtered_genes

from src.orthologs.utils import map_names


"""
Take the vanilla output from a BLAST run between Del Norte (DN) and Royal
Royce (RR) and filter it to reformat the gene and chromosome names.

The gene names from the vanilla BLAST output do not match the "Updated" gene
names that are similar to the Arabidopsis gene naming scheme.
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

    # Get the correct name for the Del_Norte genes by removing the -mRNA-
    # TODO need to check for duplicates here? Or later?
    homolog_pd["Del_Norte"] = homolog_pd["Del_Norte"].str.split("-mRNA-").str[0]

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


def blacklist_if_no_new_name(homolog_pd, decoder_ring, column="Del_Norte"):
    # Remove a Del Norte gene if it is not within the decoder ring
    homolog_pd = homolog_pd.loc[~homolog_pd[column].isin(decoder_ring)]
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
        "decoder_ring_input_file",
        type=str,
        help="parent path of decoder ring file",
    )
    parser.add_argument(
        "RR_gene_data",
        type=str,
        help="gene data file that is the input to TE Density for RR",
    )
    parser.add_argument(
        "DN_gene_data",
        type=str,
        help="gene data file that is the input to TE Density for DN",
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
    args.DN_gene_data = os.path.abspath(args.DN_gene_data)
    args.RR_gene_data = os.path.abspath(args.RR_gene_data)
    args.output_file = os.path.abspath(args.output_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    logger.info(f"Importing decoder ring...")
    decoder_ring = import_decoder_ring(args.decoder_ring_input_file)

    logger.info(f"Importing unclean homologs...")
    unclean_homologs = import_unclean_homologs(args.homolog_input_file)

    logger.info(f"Replacing names...")
    name_replaced_homologs = map_names(unclean_homologs, "Del_Norte", decoder_ring)
    name_replaced_homologs = blacklist_if_no_new_name(
        name_replaced_homologs, decoder_ring
    )

    # Rename the columns
    name_replaced_homologs.rename(
        columns={
            "Del_Norte": "DN_Gene",
            "Royal_Royce": "RR_Gene",
            "E_Value": "BLAST_E_Value",
        },
        inplace=True,
    )

    # Remove a gene from the BLAST results if it is not within the Cleaned Gene
    # Annotation
    # First read the gene annotation file
    cleaned_DN_genes = import_filtered_genes(args.DN_gene_data, logger)
    cleaned_RR_genes = import_filtered_genes(args.RR_gene_data, logger)
    # MAGIC tuple indexing, the 0th element is the column name, the first is
    # the reference to the pandaframe
    # print(name_replaced_homologs["DN_Gene"])
    for i in zip(("RR_Gene", "DN_Gene"), (cleaned_RR_genes, cleaned_DN_genes)):
        name_replaced_homologs = name_replaced_homologs.loc[
            name_replaced_homologs[i[0]].isin(i[1].index)
        ]

    # Add the chromosome IDs to the files, from the gene data
    # MAGIC tuple indexing, the 0th element is the chrom column name,
    # the first is the gene column name, and the second is the gene data
    for i in zip(
        ("DN_Chromosome", "RR_Chromosome"),
        ("DN_Gene", "RR_Gene"),
        (cleaned_DN_genes, cleaned_RR_genes),
    ):
        name_replaced_homologs[i[0]] = name_replaced_homologs[i[1]].map(
            i[2]["Chromosome"]
        )

    # Drop the rows where the chromosomes do not match
    name_replaced_homologs = name_replaced_homologs.loc[
        name_replaced_homologs["DN_Chromosome"]
        == name_replaced_homologs["RR_Chromosome"]
    ]

    logger.info(f"Saving results to disk at: {args.output_file}")
    name_replaced_homologs.to_csv(args.output_file, sep="\t", header=True, index=False)
