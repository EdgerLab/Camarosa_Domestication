#!/usr/bin/env python3

__author__ = "Scott Teresi"

"""
- Take the TE Density data of the RR and DN genomes and merge it with the
preconstructed syntelog table
"""

import pandas as pd
import numpy as np
import os
import argparse
import logging
import coloredlogs


from src.orthologs.pan_orthology_table import read_pan_orthology_table

# from src.syntelog_differences.bargraphs import graph_barplot_density_differences

from transposon.import_filtered_genes import import_filtered_genes
from transposon.gene_data import GeneData
from transposon.density_data import DensityData
from transposon.density_utils import (
    add_hdf5_indices_to_gene_data_from_list_hdf5,
    add_te_vals_to_gene_info_pandas_from_list_hdf5,
    add_te_vals_to_gene_info_pandas,
    get_specific_slice,
    add_hdf5_indices_to_gene_data,
    info_of_gene,
)


def make_table_for_te_type_and_direction(
    orthologs,
    cleaned_DN_genes,
    cleaned_RR_genes,
    processed_DN_data,
    processed_RR_data,
    major_group,
    te_type,
    direction,
    window,
):

    # Add HDF5 indices to a pandas dataframe to enable easier HDF5 TE value
    # access later
    # This is used for the special and regular gene set independently

    DN_gene_frame_with_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_DN_genes, processed_DN_data
    )
    DN_gene_frame_with_values = add_te_vals_to_gene_info_pandas_from_list_hdf5(
        DN_gene_frame_with_indices,
        processed_DN_data,
        major_group,
        te_type,
        direction,
        window,
    )

    RR_gene_frame_with_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_RR_genes, processed_RR_data
    )
    RR_gene_frame_with_values = add_te_vals_to_gene_info_pandas_from_list_hdf5(
        RR_gene_frame_with_indices,
        processed_RR_data,
        major_group,
        te_type,
        direction,
        window,
    )
    te_window_direction_str = f"{te_type}_{window}_{direction}"

    DN_gene_frame_with_values.rename(
        columns={
            te_window_direction_str: f"DN_{te_window_direction_str}",
            "Gene_Name": "DN_Gene",
        },
        inplace=True,
    )
    RR_gene_frame_with_values.rename(
        columns={
            te_window_direction_str: f"RR_{te_window_direction_str}",
            "Gene_Name": "RR_Gene",
        },
        inplace=True,
    )

    for dataframe in [DN_gene_frame_with_values, RR_gene_frame_with_values]:
        dataframe.drop(
            columns=[
                "Feature",
                "Start",
                "Stop",
                "Strand",
                "Length",
                "Index_Val",
                "Chromosome",
                "index",
            ],
            inplace=True,
        )
    big_merge = pd.merge(orthologs, RR_gene_frame_with_values, on="RR_Gene")
    bigga_merge = pd.merge(big_merge, DN_gene_frame_with_values, on="DN_Gene")

    # MAGIC
    # NOTE HARD CODED TO HAVE DN MINUS RR
    bigga_merge["Difference"] = (
        bigga_merge[f"DN_{te_window_direction_str}"]
        - bigga_merge[f"RR_{te_window_direction_str}"]
    )
    return bigga_merge


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
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.syntelog_file = os.path.abspath(args.syntelog_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    args.RR_HDF5_dir = os.path.abspath(args.RR_HDF5_dir)
    args.RR_gene_data = os.path.abspath(args.RR_gene_data)

    args.DN_gene_data = os.path.abspath(args.DN_gene_data)
    args.DN_HDF5_dir = os.path.abspath(args.DN_HDF5_dir)
    # -----------------------------------------------------------------------
    # Manually define our config
    # MAGIC
    orders = ["LTR", "TIR", "Total_TE_Density"]
    superfamilies = ["Mutator", "Copia", "Gypsy", "hAT", "Tc1-Mariner"]
    windows = [1000, 2500, 5000, 10000]
    directions = ["Upstream", "Downstream"]

    # -----------------------------------------------------------------------
    # Read and initialize various data

    # Read the universal syntelog table
    orthologs = read_pan_orthology_table(args.syntelog_file)

    # Read cleaned genes for the given genome as pandas
    cleaned_RR_genes = import_filtered_genes(args.RR_gene_data, logger)
    cleaned_DN_genes = import_filtered_genes(args.DN_gene_data, logger)

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

    # Reset index to make it easier to add the HDF5 indices to a pandas frame
    cleaned_DN_genes.reset_index(inplace=True)
    cleaned_RR_genes.reset_index(inplace=True)

    # Start looping to make the tables

    for window in windows:
        for direction in directions:
            for te_type in orders:
                major_group = "Order"
                table = make_table_for_te_type_and_direction(
                    orthologs,
                    cleaned_DN_genes,
                    cleaned_RR_genes,
                    processed_DN_data,
                    processed_RR_data,
                    major_group,
                    te_type,
                    direction,
                    window,
                )
                file_string = os.path.join(
                    args.output_dir, f"DN_minus_RR_{te_type}_{window}_{direction}.tsv"
                )
                logger.info(f"Writing to file: {file_string}")
                table.to_csv(file_string, header=True, index=False, sep="\t")

            # NOTE duplicate code...
            for te_type in superfamilies:
                major_group = "Superfamily"
                table = make_table_for_te_type_and_direction(
                    orthologs,
                    cleaned_DN_genes,
                    cleaned_RR_genes,
                    processed_DN_data,
                    processed_RR_data,
                    major_group,
                    te_type,
                    direction,
                    window,
                )
                file_string = (
                    f"{args.output_dir}/DN_minus_RR_{te_type}_{window}_{direction}.tsv"
                )
                logger.info(f"Writing to file: {file_string}")
                table.to_csv(file_string, header=True, index=False, sep="\t")
                # graph_barplot_density_differences(
                #    table["Difference"],
                #    te_type,
                #    window,
                #    direction,
                #    0,
                #    "home/scott",
                #    logger,
                #    display=True,
                # )
                # raise ValueError
