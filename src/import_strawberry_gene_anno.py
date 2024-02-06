"""
Filter a gene annotation file for the TE Density algorithm
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs


def write_cleaned_genes(gene_pandaframe, out_file, logger):
    logger.info(f"Writing cleaned gene annotation file to: {out_file}")
    gene_pandaframe.to_csv(out_file, sep="\t", header=True, index=True)


def import_genes(genes_input_path, genome_name, logger):
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

    # NOTE, this works for H4, DN, and RR, still need to check FNI, FII, and
    # FVI
    gene_pandaframe = gene_pandaframe.loc[
        ~gene_pandaframe["Chromosome"].str.contains("contig")
    ]
    gene_pandaframe = gene_pandaframe.loc[
        ~gene_pandaframe["Chromosome"].str.contains("ptg")
    ]
    gene_pandaframe = gene_pandaframe.loc[
        ~gene_pandaframe["Chromosome"].str.contains("scaf")
    ]

    # Fix the chromosome names of the Del_Norte annotation so that they
    # correspond to the chromosomes of the TE annotation
    if genome_name == "DN":
        # Fix the '_RagTag' suffix and replace the hyphen between the
        # 'subchromosomes'
        gene_pandaframe["Chromosome"] = gene_pandaframe["Chromosome"].str.replace(
            "_RagTag", ""
        )
        gene_pandaframe["Chromosome"] = gene_pandaframe["Chromosome"].str.replace(
            "-", "_"
        )

    # Fix the chromosome names of the FVI/FNI annotation so that they
    # correspond to the chromosomes of the TE annotation
    if genome_name == "FVI" or genome_name == "FNI":
        # Fix the '_RagTag' suffix and replace the hyphen between the
        # 'subchromosomes'
        gene_pandaframe["Chromosome"] = gene_pandaframe["Chromosome"].str.replace(
            "_RagTag", ""
        )

    # print(gene_pandaframe["Chromosome"].unique())

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
    parser.add_argument("genome_name", type=str)
    parser.add_argument(
        "output_file",
        type=str,
        help="path for output file",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.gene_input_file = os.path.abspath(args.gene_input_file)
    args.output_file = os.path.abspath(args.output_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Execute
    cleaned_genes = import_genes(args.gene_input_file, args.genome_name, logger)
    write_cleaned_genes(cleaned_genes, args.output_file, logger)
