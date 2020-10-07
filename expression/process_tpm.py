#!/usr/bin/env python3

"""
Generate and save an TPM table from gene annotation and count matrix
"""

__author__ = "Scott Teresi"

import argparse
import os
import numpy as np
import pandas as pd

from gene_lengths import import_genes
from count_matrix import import_count_matrix
from tpm import calc_tpm


def process(gene_annotation, count_matrix, output_dir, file_name):
    """Run the script, go from annotation and count file to TPM table.

    Args:
        gene_annotation (str): string of path to gene annotation file. GFF
        format.

        count_matrix (str): string of path to count matrix file. Shape:
            (N_Genes, N_Samples).

        output_dir (str): string of path for output directory to store data

        file_name (str): a string to be used as a prefix for the file_name,
        ought to represent whether the file is single or all haplotype in
        origin.

    Returns:
        None. Saves a table of TPM values to the output_dir
    """
    # import genes and their gene lengths
    ids_and_exon_lengths = import_genes(gene_annotation)

    # import the count matrix
    counts = import_count_matrix(count_matrix)

    # join the count matrix and gene/length matrix on the gene names column
    merged_data = pd.merge(ids_and_exon_lengths, counts, on="Gene_Name")

    # Parse out the data in prep for fpkm calculator
    merged_data.set_index("Gene_Name", inplace=True)
    lengths = merged_data.Total_Exon_Length
    merged_data.drop(columns=["Total_Exon_Length"], inplace=True)

    # Calcuate TPM
    tpm_vals = calc_tpm(merged_data, lengths)
    tpm_vals = pd.DataFrame(tpm_vals, columns=merged_data.columns)
    tpm_vals.set_index(merged_data.index, inplace=True)

    # To disk
    tpm_vals.to_csv(
        os.path.join(output_dir, f"Blueberry_TPM_{file_name}.tsv"), sep="\t"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculate gene lengths")
    path_main = os.path.abspath(__file__)
    parser.add_argument("genes_input_file", type=str, help="parent path of gene file")
    parser.add_argument("count_matrix", type=str, help="parent path of counts file")
    parser.add_argument("selection", type=str, help="all or single haplotype")
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=os.path.join(path_main, "../../../../", "Blueberry_Data/TPM"),
        help="parent directory to output results",
    )

    args = parser.parse_args()
    args.genes_input_file = os.path.abspath(args.genes_input_file)
    args.count_matrix = os.path.abspath(args.count_matrix)
    args.output_dir = os.path.abspath(args.output_dir)

    process(args.genes_input_file, args.count_matrix, args.output_dir, args.selection)
