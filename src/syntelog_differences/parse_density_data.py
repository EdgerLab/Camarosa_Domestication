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
from src.syntelog_differences.Extract_Density import Strawberry_Specific_Density

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
        "H4_gene_data",
        type=str,
        help="gene data file that is the input to TE Density for H4",
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
        "H4_HDF5_dir",
        type=str,
        help=".h5 files corresponding to the H4 genome",
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

    args.H4_gene_data = os.path.abspath(args.H4_gene_data)
    args.H4_HDF5_dir = os.path.abspath(args.H4_HDF5_dir)
    # -----------------------------------------------------------------------
    # Manually define our config
    # MAGIC
    orders = ["LTR", "TIR", "Total_TE_Density"]
    superfamilies = ["Mutator", "Copia", "Gypsy", "hAT", "CACTA"]
    windows = [1000, 2500, 5000, 10000]
    directions = ["Upstream", "Downstream"]

    # -----------------------------------------------------------------------
    # Read and initialize various data
    # TODO there has to be a better way to do this, I am repeating the same
    # code for DN, RR, and H4. This is so ugly. Maybe I can make a function to
    # clean it up but then I am just recreating DensityData?

    # Read the universal syntelog table
    orthologs = read_pan_orthology_table(args.syntelog_file)

    # NB dropping the these columns now to avoid pandas.merge conflicts later
    orthologs.drop(
        labels=["RR_Chromosome", "DN_Chromosome", "H4_Chromosome"],
        axis=1,
        inplace=True,
    )

    # Read cleaned genes for the given genome as pandas
    cleaned_RR_genes = import_filtered_genes(args.RR_gene_data, logger)
    cleaned_DN_genes = import_filtered_genes(args.DN_gene_data, logger)
    cleaned_H4_genes = import_filtered_genes(args.H4_gene_data, logger)

    # Get list of GeneData for each genome to enable initialization of
    # DensityData
    genedata_DN_list = get_gene_data_as_list(cleaned_DN_genes)
    genedata_RR_list = get_gene_data_as_list(cleaned_RR_genes)
    genedata_H4_list = get_gene_data_as_list(cleaned_H4_genes)

    # Initialize DensityData for each genome
    processed_DN_data = DensityData.from_list_gene_data_and_hdf5_dir(
        genedata_DN_list, args.DN_HDF5_dir, "DN_(.*?).h5", logger
    )
    processed_RR_data = DensityData.from_list_gene_data_and_hdf5_dir(
        genedata_RR_list, args.RR_HDF5_dir, "RR_(.*?).h5", logger
    )
    processed_H4_data = DensityData.from_list_gene_data_and_hdf5_dir(
        genedata_H4_list, args.H4_HDF5_dir, "H4_(.*?).h5", logger
    )

    DN_gene_frame_with_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_DN_genes, processed_DN_data
    )
    RR_gene_frame_with_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_RR_genes, processed_RR_data
    )
    H4_gene_frame_with_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_H4_genes, processed_H4_data
    )
    logger.info("Starting to iterate and create the sub-tables...")

    # Initialize the Strawberry Dataclass to help cut down on duplicate code
    # Start looping to make the tables
    for window in windows:
        for direction in directions:
            for te_list, major_group in [
                (orders, "Order"),
                (superfamilies, "Superfamily"),
            ]:
                for te_type in te_list:

                    H4 = Strawberry_Specific_Density(
                        H4_gene_frame_with_indices,
                        processed_H4_data,
                        "H4",
                        major_group,
                        te_type,
                        direction,
                        window,
                    )
                    DN = Strawberry_Specific_Density(
                        DN_gene_frame_with_indices,
                        processed_DN_data,
                        "DN",
                        major_group,
                        te_type,
                        direction,
                        window,
                    )
                    RR = Strawberry_Specific_Density(
                        RR_gene_frame_with_indices,
                        processed_RR_data,
                        "RR",
                        major_group,
                        te_type,
                        direction,
                        window,
                    )
                    for density, name in [(H4, "H4"), (DN, "DN"), (RR, "RR")]:
                        density.table.rename(
                            {
                                "Chromosome": f"{name}_Chromosome",
                                "Start": f"{name}_Start",
                                "Stop": f"{name}_Stop",
                                "Strand": f"{name}_Strand",
                                "Length": f"{name}_Length",
                            },
                            axis=1,
                            inplace=True,
                        )

                    # This is the general table for each genome, i.e
                    # DN_1KB_Upstream_LTR etc...
                    for i in [H4, DN, RR]:
                        i.save_table_to_disk(args.output_dir)

                    # Merge the tables with the ortholog data and subset,
                    # because we are interested in the syntelog tables at this
                    # point... We will keep merging against the ortholog table
                    big_merge = pd.merge(orthologs, RR.table, on="RR_Gene", how="left")
                    bigga_merge = pd.merge(
                        big_merge, DN.table, on="DN_Gene", how="left"
                    )
                    biggest_merge = pd.merge(
                        bigga_merge, H4.table, on="H4_Gene", how="left"
                    )

                    # MAGIC
                    # NOTE HARD CODED TO HAVE DN MINUS RR
                    te_window_direction_str = f"{te_type}_{window}_{direction}"
                    biggest_merge["Difference"] = (
                        biggest_merge[f"DN_{te_window_direction_str}"]
                        - biggest_merge[f"RR_{te_window_direction_str}"]
                    )

                    # TODO save the merged table to disk
                    outfile = os.path.join(
                        args.output_dir,
                        f"DN_minus_RR_{te_type}_{window}_{direction}.tsv",
                    )
                    logger.info(f"Writing to file: {outfile}")
                    biggest_merge.to_csv(outfile, header=True, index=False, sep="\t")
