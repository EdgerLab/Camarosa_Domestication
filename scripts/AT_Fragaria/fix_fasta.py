#!/usr/bin/env python

"""
Fix the front matter of protein FASTA files to run BLASTP
"""

__author__ = "Scott Teresi"

import argparse
import os

# import logging
import argparse


def filter_fasta(fasta_file, genome_id, output_dir):
    """
    Filter the extraneous front matter from a fasta file

    fasta file (str): Path to peptide fasta file

    genome_id (str): Identifier for the genome of the fasta file

    output_dir (str): Path to output directory
    """
    fasta_file = open(fasta_file, "r")
    output_file = open(str("fasta_revised_" + genome_id + ".fasta"), "w")

    # with open(output_file, "w") as output_file:
    for a_line in fasta_file:
        if a_line.startswith(">"):
            test = a_line.split("||")
            a_line = test[4]  # MAGIC NUMBER
            a_line = ">" + a_line + "\n"
        output_file.write(a_line)
    fasta_file.close()
    output_file.close()


if __name__ == "__main__":

    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.join(dir_main, "..", "data_dir")  # NOTE change
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "strawberry_input_file",
        type=str,
        help="parent path to peptide fasta for strawberry file",
    )

    parser.add_argument("genome_id", type=str, help="identifier for our genome at hand")

    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=output_default,
        help="parent directory to output results",
    )

    args = parser.parse_args()
    args.strawberry_input_file = os.path.abspath(args.strawberry_input_file)
    args.output_dir = os.path.abspath(args.output_dir)

    # Execute
    filter_fasta(args.strawberry_input_file, args.genome_id, args.output_dir)
