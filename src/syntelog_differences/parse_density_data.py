#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs

"""
- Take the TE Density data and merge it with the syntelog table
"""

from src.orthologs.pan_orthology_table import read_pan_orthology_table
from transposon.import_filtered_genes import import_filtered_genes
from transposon.gene_data import GeneData

from transposon.density_data import DensityData

# from src.fish_utils import (
#    add_hdf5_indices_to_gene_data_from_list_hdf5,
#    add_te_vals_to_gene_info_pandas_from_list_hdf5,
#    add_te_vals_to_syntelogs,
#    parse_analysis_config,
# )


def get_gene_data_as_list(cleaned_genes):
    """
    Take a cleaned genes annotation file from TE Density
    (import_filtered_genes) and break it into a list of GeneData objects by
    chromosome ID. This is used to initialize all of the DensityData objects in
    a list.

    Args:
        cleaned_genes (pandas.DataFrame)
            Index:
                Name: Gene_Name, strings of gene names
            Columns:
                Name: Chromosome, object
                Name: Feature, object
                Name: Start, float64
                Name: Stop, float64
                Name: Strand, object
                Name: Length, float64

    Returns:
        genedata_list (list of GeneData)
    """
    # MAGIC group by column Chromosome
    gene_dataframe_list = [
        dataframe for k, dataframe in cleaned_genes.groupby("Chromosome")
    ]

    # MAGIC initialize GeneData iteratively using the magic unique chromosome
    genedata_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in gene_dataframe_list
    ]
    return genedata_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "syntelog_file",
        type=str,
        help="syntelog file of all strawberry genes",
    )
    parser.add_argument(
        "DN_gene_data",
        type=str,
        help="gene data file that is the input to TE Density for DN",
    )
    parser.add_argument(
        "RR_gene_data",
        type=str,
        help="gene data file that is the input to TE Density for RR",
    )

    parser.add_argument(
        "DN_HDF5_dir",
        type=str,
        help=".h5 files corresponding to the DN genome",
    )
    parser.add_argument(
        "RR_HDF5_dir",
        type=str,
        help=".h5 files corresponding to the RR genome",
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Parent directory to output results",
    )
    logger = logging.getLogger(__name__)
    args = parser.parse_args()
    args.syntelog_file = os.path.abspath(args.syntelog_file)
    args.output_dir = os.path.abspath(args.output_dir)

    args.RR_HDF5_dir = os.path.abspath(args.RR_HDF5_dir)
    args.RR_gene_data = os.path.abspath(args.RR_gene_data)

    args.DN_gene_data = os.path.abspath(args.DN_gene_data)
    args.DN_HDF5_dir = os.path.abspath(args.DN_HDF5_dir)
    # -----------------------------------------------------------------------
    # Read and initialize various data

    # Read the universal syntelog table
    orthologs = read_pan_orthology_table(args.syntelog_file)

    # Read cleaned genes for the given genome as pandas
    cleaned_RR_genes = import_filtered_genes(args.RR_gene_data, logger)
    cleaned_DN_genes = import_filtered_genes(args.DN_gene_data, logger)

    print(orthologs)
    print(cleaned_RR_genes)
    print(cleaned_DN_genes)

    # Get list of GeneData for each genome to enable initialization of
    # DensityData
    genedata_DN_list = get_gene_data_as_list(cleaned_DN_genes)
    genedata_RR_list = get_gene_data_as_list(cleaned_RR_genes)

    # Initialize DensityData for each genome
    processed_DN_data = DensityData.from_list_gene_data_and_hdf5_dir(
        genedata_DN_list, args.DN_HDF5_dir, "DN_(.*?).h5", logger
    )
    processed_RR_data = DensityData.from_list_gene_data_and_hdf5_dir(
        genedata_RR_list, args.RR_HDF5_dir, "RR_(.*?).h5", logger
    )
    print(processed_DN_data[0])
    print(processed_RR_data[0])
    raise NotImplementedError("This is not finished")

    # Reset index to make it easier to add the HDF5 indices to a pandas frame
    cleaned_genes.reset_index(inplace=True)

    # Add HDF5 indices to a pandas dataframe to enable easier HDF5 TE value
    # access later
    # This is used for the special and regular gene set independently
    gene_frame_with_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_genes, processed_dd_data
    )

    # -----------------------------------------------------------------------
    # Begin analysis of regular genes

    # Do the loops to create graphs for the general set
    for window in windows:
        for direction in directions:
            for order in orders:

                # NOTE this is only for the orders. I could do it for the
                # superfamilies but I don't think that is wanted or needed
                cleaned_with_te_vals = add_te_vals_to_gene_info_pandas_from_list_hdf5(
                    gene_frame_with_indices,
                    processed_dd_data,
                    "Order",
                    order,
                    direction,
                    window,
                )

                syntelogs_w_te_vals = add_te_vals_to_syntelogs(
                    syntelogs, cleaned_with_te_vals, args.genome_name
                )

                # MAGIC string formatting
                # NOTE A genome MINUS B genome. Negative values means that the
                # B copy had more TE
                syntelogs_w_te_vals["Difference"] = (
                    syntelogs_w_te_vals[
                        order + "_" + str(window) + "_" + direction + "_A"
                    ]
                    - syntelogs_w_te_vals[
                        order + "_" + str(window) + "_" + direction + "_B"
                    ]
                )

                total_length = len(syntelogs_w_te_vals)

                # NB subset to have only rows with a nonzero difference

                syntelogs_w_te_vals = syntelogs_w_te_vals.loc[
                    syntelogs_w_te_vals["Difference"] != 0
                ]
                number_of_zeros = total_length - len(syntelogs_w_te_vals)

                graph_histogram_density_differences(
                    syntelogs_w_te_vals["Difference"].to_list(),
                    order,
                    window,
                    direction,
                    number_of_zeros,
                    args.genome_name,
                    args.genome_name + "_A",
                    args.genome_name + "_B",
                    args.output_dir,
                    logger,
                    display=False,
                )
