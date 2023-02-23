"""
Filter a gene annotation file for the TE Density algorithm
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs


def write_cleaned_genes(gene_pandaframe, output_dir, old_filename, logger):
    file_name = os.path.join(
        output_dir,
        ("Cleaned_" + os.path.splitext(os.path.basename(old_filename))[0]) + ".tsv",
    )  # MAGIC to get proper extension

    logger.info("Writing cleaned gene annotation file to: %s" % file_name)
    gene_pandaframe.to_csv(file_name, sep="\t", header=True, index=True)


def import_genes(genes_input_path, logger):
    """Import gene annotation."""

    col_names = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Score",
        "Strand",
        "Frame",
        "FullName",
    ]

    col_to_use = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Strand",
        "FullName",
    ]

    gene_pandaframe = pd.read_csv(
        genes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        dtype={"Stop": "float64", "Start": "float64", "Strand": str, "FullName": str},
        comment="#",
    )

    gene_pandaframe = gene_pandaframe[
        gene_pandaframe.Feature == "gene"
    ]  # drop non-gene rows in annotation

    # clean the names and set as the index (get row wrt name c.f. idx)
    gene_pandaframe["Gene_Name"] = gene_pandaframe["FullName"].str.extract(r"ID=(.*?);")
    gene_pandaframe.set_index("Gene_Name", inplace=True)

    # Drop extraneous columns
    gene_pandaframe = gene_pandaframe.drop(columns=["FullName", "Software"])

    # Edit gene length
    gene_pandaframe["Length"] = gene_pandaframe.Stop - gene_pandaframe.Start + 1

    gene_pandaframe.sort_values(by=["Chromosome", "Start"], inplace=True)

    return gene_pandaframe


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat gene annotation file")

    parser.add_argument(
        "gene_input_file", type=str, help="Parent path of gene annotation file"
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Parent directory to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.gene_input_file = os.path.abspath(args.gene_input_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Execute
    cleaned_genes = import_genes(args.gene_input_file, logger)
    write_cleaned_genes(cleaned_genes, args.output_dir, args.gene_input_file, logger)
