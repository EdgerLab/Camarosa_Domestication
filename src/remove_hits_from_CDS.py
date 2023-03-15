#!/usr/bin/env python3

"""
Name: Filter CDS FASTA for EDTA

Problem: 1. EDTA expects gene names to be short. Most gene names are not short.
         2. Via a previously performed BLAST search, I have genes that have
            strong similarity to TEs. They must be removed from the CDS FASTA
            before I run EDTA to reduce "contamination", and improve annotation
            quality.
Solution:
    Arbitrarily rename the gene names to comply with EDTA. Remove the TE-like
    genes from the CDS FASTA. Write a new FASTA with both of these elements
    included. Create a simplified list of the TE-like genes in order to
    simplify later analysis.

Things to consider: I have written several variations of this script over the
    past year. Each script involves some element of reading the FASTA file
    within the Bio.SeqIO package and trimming/eliminating sequence IDs (genes).
    Every time I need to tailor some aspect of this to the string format of the
    genes. Future work might consider moving this to a class, to reduce the
    amount of copy-paste code between projects. Although, this is heavily
    tailored to using EDTA software, so it may become obsolete once a different
    annotation software is chosen.

NOTE: 4/19/2023 @ 2:17 pm, validated that the correct number of genes are being
dropped from Del Norte
"""

__author__ = "Scott Teresi"

import argparse
import os
import csv
import logging
import coloredlogs

from Bio import SeqIO
import pandas as pd
import numpy as np
import re


def load_blast_hits(blast_hit_file, genome_name):
    """
    Read the output from BLAST of Ning's Tpases file against the Royal Royce
    genome

    Args:
        blast_hit_file (str): Filepath to the output from a BLAST calculation,
        format 8

        genome_name (str): String of the genome name, needed for custom
        filtering

    Returns:
        blast_hit_pandas (Pandas.Data.Frame): Pandas object

            Index:
                RangeIndex
            Columns: More information found within BLAST documentation
                Query_ID: TE name, dtype: str
                Subject_ID: Gene name, dtype: str
                Percent_Match: percent match, dtype: float32
                Alignment_Length: the length of the alignment, dtype: int
                Mismatches: the number of mismatches, dtype: int
                Gap_Openings: dtype: int
                Query_Start: position where the TE starts, dtype: int
                Query_End: position where the TE ends, dtype: int
                Subject_Start: position where the gene starts, dtype: int
                Subject_End: position where the gene ends, dtype: int
                E_Value: e value, dtype: float
                Bit_Score: bit score, dtype: int
            Shape: For RR (309129 x 12)

    """
    blast_hit_pandas = pd.read_csv(
        blast_hit_file,
        header=None,
        sep="\t",
        names=[
            "Query_ID",
            "Subject_ID",
            "Percent_Match",
            "Alignment_Length",
            "Mismatches",
            "Gap_Openings",
            "Query_Start",
            "Query_End",
            "Subject_Start",
            "Subject_End",
            "E_Value",
            "Bit_Score",
        ],
        dtype={
            "Query_ID": str,
            "Subject_ID": str,
            "Percent_Match": float,
            "Alignment_Length": int,
            "Mismatches": int,
            "Gap_Openings": int,
            "Query_Start": int,
            "Query_End": int,
            "Subject_Start": int,
            "Subject_End": int,
            "E_Value": float,
            "Bit_Score": float,
        },
    )
    return blast_hit_pandas


def reformat_cds_seq_iq(
    input_fasta_filepath,
    new_fasta_filepath,
    blacklist_genes,
    genome_name,
    output_dir,
    logger,
):
    """
    Reformat a CDS FASTA file to remove genes that match with Ning's tpase
    file. ALSO reformat the genes to comply with EDTA's requirement that gene
    names aren't 13 characters or greater.

    Args:
        input_fasta (str): String path to input fasta file
        new_fasta_filepath (str): String path to output fasta file
        blacklist_genes (list of str): List of string gene names

        logger (logging.Logger): Object to log information to

    Returns:
        None, just saves the edited FASTA file to disk. Also writes a
        conversion table to disk for the old names and their new name
        counterparts
    """
    # MAGIC file suffixes
    name_key = os.path.join(
        output_dir,
        (genome_name + "_CDS_Seq_ID_Conversion.txt"),
    )
    shortened_gene_name_dict = {}
    # NB this is used to write the conversion key later for
    # clarity

    count = 0
    with open(input_fasta_filepath, "r") as input_fasta:
        with open(new_fasta_filepath, "w") as new_fasta_output:
            for s_record in SeqIO.parse(input_fasta, "fasta"):
                # NB the s_record.id and s_record.description combined contain
                # all the information for each entry following the '>' character
                # in the fasta
                if s_record.id in blacklist_genes:
                    continue  # Just start the next iteration, doesn't write
                # the gene

                shortened_gene_name_dict[s_record.id] = "gene_" + str(count)
                s_record.id = "gene_" + str(count)
                count += 1

                s_record.description = ""  # NB edit the description so that when
                # we rewrite we don't have the extraneous info
                SeqIO.write(s_record, new_fasta_output, "fasta")
    logger.info("Finished writing new fasta to: %s" % new_fasta_filepath)

    # Write the conversion table for record-keeping.
    header = ["Original_Name", "New_Name"]
    with open(name_key, "w", newline="") as output:
        writer = csv.writer(output, lineterminator=os.linesep)
        writer.writerow(header)
        for key, val in shortened_gene_name_dict.items():
            writer.writerow([key, val])
    logger.info(
        "Finished writing name conversion table to: %s"
        % os.path.join(output_dir, name_key)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove TE-like genes from CDS FASTA")

    parser.add_argument("fasta_input_file", type=str, help="parent path of fasta file")
    parser.add_argument(
        "blast_input_file",
        type=str,
        help="parent path of blast table file, format form 8",
    )
    parser.add_argument("new_cds_file", type=str, help="parent path of fasta file")
    parser.add_argument("genome_name", type=str, help="string of genome name")
    parser.add_argument(
        "output_dir",
        type=str,
        help="Parent directory to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.fasta_input_file = os.path.abspath(args.fasta_input_file)
    args.blast_input_file = os.path.abspath(args.blast_input_file)
    args.new_cds_file = os.path.abspath(args.new_cds_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    blast_pandas = load_blast_hits(args.blast_input_file, args.genome_name)
    blacklist_genes = blast_pandas["Subject_ID"].unique().tolist()
    num_unique_genes = len(blast_pandas["Subject_ID"].unique())

    with open(
        os.path.join(args.output_dir, args.genome_name + "_Blacklist.tsv"), "w"
    ) as f_in:
        for gene in blacklist_genes:
            f_in.write("%s\n" % gene)
        logger.info("Wrote unique blacklisted genes to: %s" % f_in)

    print(f"Genome Name: {args.genome_name}")
    print(f"Number of unique blacklisted genes: {num_unique_genes}")
    print()

    reformat_cds_seq_iq(
        args.fasta_input_file,
        args.new_cds_file,
        blacklist_genes,
        args.genome_name,
        args.output_dir,
        logger,
    )
