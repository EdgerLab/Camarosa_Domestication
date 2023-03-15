#!/usr/bin/env python3

"""
Name: Trim (edit) gene names in CDS for BLAST and downstream analyses

Problem: The gene names within the various strawberry CDS FASTA files are
    rather long and sometimes improperly formatted. E.g the default Del Norte file
    causes BLAST to fail because it can't extract the right ID for each sequence
    entry. Additionally, messy CDS FASTA files complicates downstream analysis.

Solution: Meaningfully shorten the gene names within the CDS FASTA without losing
        useful information.
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


def fix_sequence_ID_names(
    input_fasta_filepath, output_fasta_filepath, genome_name, logger
):
    """
    input_fasta_filepath (str): String path to input fasta file
    output_fasta_filepath (str): String path to output fasta file
    genome_name (str): String representing genome name, so that we can easily
        apply a regex rule for different genomes
    """
    with open(input_fasta_filepath, "r") as input_fasta_filepath:
        with open(output_fasta_filepath, "w") as output_fasta_filepath:
            for s_record in SeqIO.parse(input_fasta_filepath, "fasta"):

                # Rules for each genome below
                if genome_name == "DN":
                    # NOTE, this has been validated 4/18/2023 @ 1:30 pm. It
                    # works and faithfully trims all the excess string
                    # information from the DN genes. It leaves you with:
                    # Fvb7-4g20850-mRNA-1

                    # NB, the -mRNA- suffix can have more than one digit at the
                    # end
                    regex_rule = r"\|\|.+?\|\|.+?\|\|.+?\|\|(.+?-mRNA-\d+)"
                    match = re.search(regex_rule, s_record.description)
                    if match:
                        s_record.id = match.group(1)

                if genome_name in ["RR", "H4", "FVI", "FNI", "FII"]:
                    # NOTE, this has been validated 4/18/2023 @ 2:00 pm. It
                    # works and faithfully trims all the excess string
                    # information. It leaves in the isoforms of the genes

                    # NOTE, it doesn't actually do anything for FII since that
                    # genome is already minimally formatted
                    regex_rule = r"^([^\s]+)"
                    match = re.search(regex_rule, s_record.description)
                    if match:
                        s_record.id = match.group(1)

                s_record.description = ""  # NB edit the description so that when
                # we rewrite we don't have the extraneous info
                SeqIO.write(s_record, output_fasta_filepath, "fasta")
    logger.info("Finished writing new fasta to: %s" % output_fasta_filepath)


if __name__ == "__main__":
    # TODO add general description
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument("fasta_input_file", type=str, help="parent path of fasta file")
    parser.add_argument("new_cds_file", type=str, help="parent path of fasta file")
    parser.add_argument("genome_name", type=str, help="TODO")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.fasta_input_file = os.path.abspath(args.fasta_input_file)
    args.new_cds_file = os.path.abspath(args.new_cds_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    fix_sequence_ID_names(
        args.fasta_input_file,
        args.new_cds_file,
        args.genome_name,
        logger,
    )
