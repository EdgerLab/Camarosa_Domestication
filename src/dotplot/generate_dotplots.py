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
from collections import defaultdict, namedtuple
from configparser import ConfigParser
import ast
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

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

from src.syntelog_differences.parse_density_data import get_gene_data_as_list
from src.orthologs.pan_orthology_table import read_pan_orthology_table


def read_dotplot_table(filepath):
    data = pd.read_csv(filepath, sep="\t", header="infer", index_col=0)
    return data


def parse_dotplot_config(config_file):
    """Return parameters for dotplot generation from a config file"""
    parser = ConfigParser()
    parser.read(config_file)
    orders = ast.literal_eval(parser.get("dotplot_parameters", "orders"))
    superfamilies = ast.literal_eval(parser.get("dotplot_parameters", "superfamilies"))
    directions = ast.literal_eval(parser.get("dotplot_parameters", "directions"))
    window_start = parser.getint("dotplot_parameters", "first_window_size")
    window_step = parser.getint("dotplot_parameters", "window_delta")
    window_stop = parser.getint("dotplot_parameters", "last_window_size")
    dotplot_param = {
        "orders": orders,
        "superfamilies": superfamilies,
        "window_range": range(window_start, window_stop + 1, window_step),
        "directions": directions,
    }
    return dotplot_param


def generate_plotting_dict(table, te_groups, windows, directions):
    """
    Generate a dictionary that can be used to plot the data.
    Dictionary is of the form: {TE_Type_Direction: [Mean_TE_Density]}, i.e
    LTR_Upstream: [0.1, 0.2, 0.3], and values are in the order of the windows
    provided.

    Parameters:
        table: pandas dataframe, the table that contains the TE density data
        te_groups: list of strings, the TE types to be plotted
        windows: list of ints, the window sizes to be plotted
        directions: list of strings, the directions to be plotted
    Returns:
        plotting_dict: defaultdict, a dictionary that can be used to plot the
        data, described above
    """
    plotting_dict = defaultdict(list)
    for te_type in te_groups:
        for direction in directions:
            for window in windows:
                mean_val = table[te_type + "_" + str(window) + "_" + direction].mean()
                plotting_dict[te_type + "_" + direction].append(mean_val)
    return plotting_dict


def generate_plotting_table(
    te_groups,
    windows,
    directions,
    processed_data,
    gene_frame_with_indices,
    order_or_super,
):
    """
    Generate a table that can be used to plot the data.

    Parameters:
        te_groups: list of strings, the TE types to be plotted
        windows: list of ints, the window sizes to be plotted
        directions: list of strings, the directions to be plotted
        processed_data: list of DensityData
        gene_frame_with_indices: pandas dataframe, the gene data with indices
            Format:
                Gene_Name (str)
                Chromosome (str)
                Start (Float or Int)
                End (Float or Int)
                Strand (str)
                Length (Float or Int)
                TE_Index (int): index in the processed density data
        order_or_super: string, either "Order" or "Superfamily"


    Returns:
        gene_frame_with_values (pandas dataframe): Same format as
             gene_frame_with_indices, but we have added the TE columns to it,
             i.e LTR_1000_Upstream, LTR_2000_Upstream, etc. This does it for
             ALL TE types for a given Order or Superfamily
    """

    i = 0
    for te_type in te_groups:
        for direction in directions:
            for window in windows:
                if i == 0:
                    gene_frame_with_values = (
                        add_te_vals_to_gene_info_pandas_from_list_hdf5(
                            gene_frame_with_indices,
                            processed_data,
                            order_or_super,
                            te_type,
                            direction,
                            window,
                        )
                    )
                if i != 0:
                    gene_frame_with_values = (
                        add_te_vals_to_gene_info_pandas_from_list_hdf5(
                            gene_frame_with_values,
                            processed_data,
                            order_or_super,
                            te_type,
                            direction,
                            window,
                        )
                    )
                i += 1
    return gene_frame_with_values


def plot_dotplot(plotting_dict, plot_name, output_dir, logger):
    """
    Generate the dotplot and save to disk

    Parameters:
        plotting_dict: defaultdict, a dictionary that can be used to plot the
        data, described in `generate_plotting_dict()`
        plot_name: string, the name of the plot
        output_dir: string, the directory to save the plot
        logger: logging object
    """
    # Make a figure with 2 subplots, subplot 1 is upstream, 2 downstream
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey="col")
    fig.set_size_inches(16, 9.5)  # MAGIC, we want the figure to be big

    # Plot the data
    for key, val in plotting_dict.items():
        if "Upstream" in key:
            ax1.plot(
                windows,
                val,
                label=key,
                linestyle=(0, (3, 1, 1, 1)),
                marker="o",
            )

        # NOTE, a little hacky, using the downstream plot for our legend,
        # so I am going to modify the string to remove the direction
        if "Downstream" in key:
            shortened_label = key.replace("_Downstream", "")
            shortened_label = shortened_label.replace("_", " ")
            ax2.plot(
                windows,
                val,
                label=shortened_label,
                linestyle=(0, (3, 1, 1, 1)),
                marker="o",
            )

    # Set the y-axis to be a percentage
    for ax in [ax1, ax2]:
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))

    # Set the labels for the upstream and downstream subplots
    ax1.set(
        yticks=np.arange(0, 0.375, 0.025),  # MAGIC
        xlabel="BP Upstream",
        ylabel="TE Density",
        xlim=[max(windows), min(windows)],
        xticks=range(min(windows), (max(windows) + 1), 1000),
        ylim=[0.0, 0.35],  # MAGIC
    )
    ax2.set(
        yticks=np.arange(0, 0.375, 0.025),  # MAGIC
        xlabel="BP Downstream",
        xlim=[min(windows), max(windows)],
        xticks=range(min(windows), (max(windows) + 1), 1000),
        ylim=[0.0, 0.35],  # MAGIC
    )
    ax2.yaxis.tick_right()  # MAGIC put ticks on the right hand side
    ax2.legend(loc="upper right")
    fig.suptitle(
        "Average TE Density of All Genes as a Function of Window Size and Location"
    )

    file_out = os.path.join(args.output_dir, plot_name)
    logger.info(f"Saving graphic to: {file_out}")
    plt.savefig(file_out)
    plt.close()


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="Generate dotplots for TE density")

    parser.add_argument(
        "strawberry_hdf5_dir",
        type=str,
        help="parent path to a dir of TE density data",
    )
    parser.add_argument(
        "strawberry_gene_data",
        type=str,
        help="parent path to a strawberry genome's filtered gene data file",
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="parent directory to output results",
    )
    parser.add_argument("ortholog_table", type=str, help="path to the ortholog table")

    parser.add_argument("config_file", type=str, help="path to the config file")

    parser.add_argument(
        "hdf5_string",
        type=str,
        help="corresponding access string to parse the hdf5 file",
    )
    parser.add_argument(
        "genome", type=str, help="The name of the genome to be analyzed"
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.strawberry_hdf5_dir = os.path.abspath(args.strawberry_hdf5_dir)
    args.strawberry_gene_data = os.path.abspath(args.strawberry_gene_data)
    args.output_dir = os.path.abspath(args.output_dir)
    args.config_file = os.path.abspath(args.config_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # -----------------------------------------------------------------------
    # Read and initialize various data
    # Read cleaned genes for the given genome as pandas
    cleaned_genes = import_filtered_genes(args.strawberry_gene_data, logger)

    # Read the ortholog table
    orthologs = read_pan_orthology_table(args.ortholog_table)

    # Get list of GeneData for the genome to enable initialization of
    # DensityData
    genedata_list = get_gene_data_as_list(cleaned_genes)

    # Initialize DensityData for the genome
    processed_data = DensityData.from_list_gene_data_and_hdf5_dir(
        genedata_list, args.strawberry_hdf5_dir, args.hdf5_string, logger
    )

    # Add indices to the gene data, so we can get the TE values easier
    gene_frame_with_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_genes, processed_data
    )

    # Set our paramaters for the dotplot from the config file
    dotplot_parameters = parse_dotplot_config(args.config_file)

    # Redefine some of the parameters as a namedtuple to make it easier to keep
    # track of which ones are Orders and Superfamilies, makes references easier
    TE_Set = namedtuple("TE_Set", ["te_list", "te_type"])
    orders = TE_Set(dotplot_parameters["orders"], "Order")
    superfamilies = TE_Set(dotplot_parameters["superfamilies"], "Superfamily")
    windows = [i for i in dotplot_parameters["window_range"]]
    directions = dotplot_parameters["directions"]
    # -----------------------------------------------------------------------

    # Add superfamilies to this list if wanted
    for tes in [orders]:
        logger.info("Generating plotting table...")
        table_to_plot = generate_plotting_table(
            tes.te_list,
            windows,
            directions,
            processed_data,
            gene_frame_with_indices,
            tes.te_type,
        )
        logger.info("Generating plotting dictionary...")
        plotting_dict = generate_plotting_dict(
            table_to_plot, tes.te_list, windows, directions
        )
        plot_name = f"{args.genome}_Dotplot_Strawberry_AllGenes_{tes.te_type}.png"

        # Save the dict to disk so we can use it for other analyses
        data = pd.DataFrame.from_dict(plotting_dict)
        data.index = windows
        table_out = os.path.join(
            args.output_dir,
            f"{args.genome}_Dotplot_DataFrame_Strawberry_{tes.te_type}.tsv",
        )
        data.to_csv(table_out, sep="\t", header=True, index=True)

        # Print the last entry in the values of the plotting dict
        logger.info("All genes plotting dict:")
        for key, val in plotting_dict.items():
            logger.info(f"\tLast value of {key}: {val[-1]}")
        plot_dotplot(plotting_dict, plot_name, args.output_dir, logger)

        # -----------------------------------------------------------------------
        # Plot the dotplot, but this time only include genes that have an AT
        # ortholog
        logger.info("Starting to plot the dotplot for genes with orthologs...")
        orthologs_to_plot = table_to_plot.merge(
            orthologs, left_on="Gene_Name", right_on=f"{args.genome}_Gene", how="inner"
        )

        # Make it so that we only have rows where the Arabidopsis_Gene column
        # is not null
        orthologs_to_plot = orthologs_to_plot.dropna(subset=["Arabidopsis_Gene"])

        # Redefine the plotting dict
        plotting_dict = generate_plotting_dict(
            orthologs_to_plot, tes.te_list, windows, directions
        )
        logger.info("Ortholog plotting dict:")

        for key, val in plotting_dict.items():
            logger.info(f"\tLast value of {key}: {val[-1]}")
        plot_name = f"{args.genome}_Dotplot_Strawberry_Orthologs_{tes.te_type}.png"
        plot_dotplot(plotting_dict, plot_name, args.output_dir, logger)
