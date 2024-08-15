#!/usr/bin/env/python

"""
Generate a dotplot of TE Density for upstream/downstream for Del Norte and
Royal Royce in ONE figure, using Matplotlib subplots
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import ast
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from src.dotplot.generate_dotplots import read_dotplot_table


def plot_dual_dotplot(
    RR_plotting_table, DN_plotting_table, plot_name, output_dir, logger
):
    """
    NOTE a lot of this is similar to and copied from `generate_dotplot.py`
    Generate the dotplot and save to disk

    Parameters:
        plotting_dict: defaultdict, a dictionary that can be used to plot the
        data, described in `generate_plotting_dict()`
        plot_name: string, the name of the plot
        output_dir: string, the directory to save the plot
        logger: logging object
    """
    windows = RR_plotting_table.index.to_list()
    # Make a figure with 2 subplots, subplot 1 is upstream, 2 downstream
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey="col")
    fig.set_size_inches(16, 9.5)  # MAGIC, we want the figure to be big

    for table, name in zip([RR_plotting_table, DN_plotting_table], ["RR", "DN"]):
        plotting_dict = table.to_dict("list")

        # Plot the data
        for key, val in plotting_dict.items():
            if "Upstream" in key:
                ax1.plot(
                    windows,
                    val,
                    linestyle=(0, (3, 1, 1, 1)),
                    marker="o",
                )

            # NOTE, a little hacky, using the downstream plot for our legend,
            # so I am going to modify the string to remove the direction
            if "Downstream" in key:
                shortened_label = key.replace("_Downstream", "")
                shortened_label = shortened_label.replace("_", " ")
                shortened_label = name + " " + shortened_label
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
        ax2.legend(loc="upper left")
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
        "RR_dotplot_table", type=str, help="path to the dotplot table for Royal Royce"
    )
    parser.add_argument(
        "DN_dotplot_table", type=str, help="path to the dotplot table for Del Norte"
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="parent directory to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.RR_dotplot_table = os.path.abspath(args.RR_dotplot_table)
    args.DN_dotplot_table = os.path.abspath(args.DN_dotplot_table)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # -----------------------------------------------------------------------
    # Read and initialize various data
    RR = read_dotplot_table(args.RR_dotplot_table)
    DN = read_dotplot_table(args.DN_dotplot_table)

    print(RR)
    print(DN)

    plot_dual_dotplot(RR, DN, "TE_Density_Dual_Dotplot.png", args.output_dir, logger)
