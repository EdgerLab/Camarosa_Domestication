#!/usr/bin/env python3

__author__ = "Scott Teresi"

from Bio import SeqIO
from Bio.Seq import Seq

import os
import argparse
import logging
import coloredlogs

"""
Translate a nucleotide FASTA file to a protein FASTA file
"""


def translate(input_file, output_file):
    # Open the input and output files
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Iterate through each sequence in the input file
        for record in SeqIO.parse(infile, "fasta"):

            # Create a Bio.Seq object from the nucleotide sequence
            nucleotide_seq = Seq(record.seq)

            # Translate the nucleotide sequence to protein sequence
            protein_seq = nucleotide_seq.translate()

            # Write the protein sequence to the output file in FASTA format
            outfile.write(f">{record.id}\n{protein_seq}\n")


if __name__ == "__main__":
    # Parse command line arguments
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(
        description="""Import unclean gene expression data previously
        downloaded from Edger Lab"""
    )

    parser.add_argument("input_nucleotides", type=str, help="TODO")
    parser.add_argument(
        "output_nucleotides",
        type=str,
        help="TODO",
    )
    args = parser.parse_args()
    args.input_nucleotides = os.path.abspath(args.input_nucleotides)
    args.output_nucleotides = os.path.abspath(args.output_nucleotides)

    translate(args.input_nucleotides, args.output_nucleotides)
