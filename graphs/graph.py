#!/usr/bin/env/python

"""
Execute graphing commands
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd
import h5py

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

from scripts.density_data import DensityData
from graphs.dotplots import *
from transposon.gene_data import GeneData


def test_data():
    """remove at later point"""
    # TEST DATA
    table = np.array(range(30))
    table = table.reshape(3, 2, 5)
    # Functionally equivalent
    # table[1][1][:]
    # table[1, 1, :]


def _split_wrt_chromosome(filtered_genes):
    """Segment data frames with respect to chromosome.

    Data is split wrt chromosome as each chromosome is processed indepenently.

    Args:
        filtered_genes(pandas.DataFrame): preprocessed genes
    Returns:
        list(pandas.DataFrame): gene frames
    """
    # NOTE copied from transposon.preprocess, will need to relocate later
    # perhaps when DensityData is refactored into the TE density code.
    # I dropped the TE portion of the code since it was doing it in tuples of
    # genes and TEs, but for the density data I only need the genes.
    # I also had to remove the validate split command

    group_key = "Chromosome"  # MAGIC NUMBER our convention
    gene_groups = filtered_genes.groupby(group_key)
    gene_list = [gene_groups.get_group(g) for g in gene_groups.groups]
    return gene_list


def plot_all_H4(h5_file_list, gene_pandaframe, logger, output_default):
    genes_split = _split_wrt_chromosome(gene_pandaframe)
    genes_split = [
        GeneData(geneframe, ("H4_" + geneframe.Chromosome.unique()[0]))
        for geneframe in genes_split
    ]
    chromosomes = [chromosome for chromosome in genes_split]  # gene data
    dd_objects = [
        DensityData(h5_file, chromosome, logger)
        for h5_file, chromosome, in zip(h5_file_list, chromosomes)
    ]
    for chromosome_of_density_data in dd_objects:
        plot_density_all(chromosome_of_density_data, "Superfamily", output_default)
        plot_density_all(chromosome_of_density_data, "Order", output_default)


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    # output_default = "/home/scott/Documents/Uni/Research/Projects/Domestication_Data/"
    output_default = os.path.abspath(
        os.path.join(dir_main, "../../", "Domestication_Data")
    )
    parser = argparse.ArgumentParser(description="generate graphs")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=output_default,
        help="parent directory to output results",
    )
    args = parser.parse_args()
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Below contains an option for the 1 chromosome set
    # test_file = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/ShortTest_Fvb1-1.h5"
    test_file1 = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/H4/5KB_Up/H4_Fvb1.h5"
    test_file2 = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/H4/5KB_Up/H4_Fvb2.h5"
    test_file3 = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/H4/5KB_Up/H4_Fvb3.h5"
    test_file4 = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/H4/5KB_Up/H4_Fvb4.h5"
    test_file5 = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/H4/5KB_Up/H4_Fvb5.h5"
    test_file6 = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/H4/5KB_Up/H4_Fvb6.h5"
    test_file7 = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/H4/5KB_Up/H4_Fvb7.h5"

    # NOTE I need to split the GeneData up by chromosome, I am pretty sure I
    # already have a function that does that...
    # Below contains an option for the 1 chromosome set
    all_genes = pd.read_csv(
        # "/home/scott/Documents/Uni/Research/Projects/TE_Data/filtered_input_data/Cleaned_ShortFvb1-1_Genes.tsv",
        "/home/scott/Documents/Uni/Research/Projects/TE_Data/filtered_input_data/Cleaned_H4_Genes.tsv",
        sep="\t",
        header="infer",
        dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
        index_col="Gene_Name",  # this is crucial
    )

    plot_all_H4(
        [
            test_file1,
            test_file2,
            test_file3,
            test_file4,
            test_file5,
            test_file6,
            test_file7,
        ],
        all_genes,
        logger,
        output_default,
    )
    # raise ValueError
