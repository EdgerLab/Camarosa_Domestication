#!/usr/bin/env python3

"""
Generate dotplot
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_density_intra(te_type, output_dir, display=False):
    """
    Args:
        te_type (str): String representing the category of TE, should be either
        'Order' or 'Superfamily', raises ValueError if it does not match.

        output_dir (str): Path representing output directory.

        display (str): Defaults to False, if True it will show the graph
        (usually opens another window).

    """
    test_densities_ltr = [0.3]
    test_densities_line = [0.8]
    data = {"ltr": test_densities_ltr, "line": test_densities_line}

    for key, val in data.items():
        plt.scatter(
            0, val, label=str(key),
        )
        plt.legend()

    plt.yticks(np.arange(0, 1.01, 0.1))
    plt.xticks([0])
    plt.xlabel("Window Position in BP")
    plt.ylabel("Intronic TE Density")

    if te_type == "Order":
        plt.title("TE Orders")
    elif te_type == "Superfamily":
        plt.title("TE Superfamilies")
    else:
        raise ValueError("Please provide Order or Superfamily")
    plt.savefig(os.path.join(output_dir, str(te_type + "_Intra_Plot.png")))

    if display:
        plt.show()


def plot_density_all(te_type, output_dir, display=False):
    """
    Plot the mean density value of each TE category over each window for "left"
    and "right".

    Args:
        te_type (str): String representing the category of TE, should be either
        'Order' or 'Superfamily', raises ValueError if it does not match.

        output_dir (str): Path representing output directory.

        display (str): Defaults to False, if True it will show the graph
        (usually opens another window).


    Returns:
        Plot. Also saves the plot

    """
    # Get the mean density value for all genes for each TE type. E.g LTR
    # density at 500 distance, take the mean of that for all genes, that will
    # be the value plotted. Do this for Order and SuperFamily

    # TODO:
    # Funciton to check the density values being provided to graph are
    # the same length as the number of windows.

    # Should provide it a dictionary of lists
    test_densities_ltr_forward = [0.3, 0.6, 0.5, 0.2]
    test_densities_line_forward = [0.8, 0.9, 0.2, 0.3]
    test_densities_ltr_reverse = [0.3, 0.4, 0.42, 0.40]
    test_densities_line_reverse = [0.65, 0.9, 0.91, 0.93]
    test_densities_ltr_intra = [0.1]
    test_densities_line_intra = [0.4]
    test_windows = [500, 1000, 1500, 2000]

    data_intra = {"ltr": test_densities_ltr_intra, "line": test_densities_line_intra}

    data_forward = {
        "ltr": test_densities_ltr_forward,
        "line": test_densities_line_forward,
    }

    data_reverse = {
        "ltr": test_densities_ltr_reverse,
        "line": test_densities_line_reverse,
    }

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)

    if te_type == "Order":
        fig.suptitle("Mean TE Density of Genes by TE Type: " + str(te_type))
    elif te_type == "Superfamily":
        fig.suptitle("Mean TE Density of Genes by TE Type: " + str(te_type))
    else:
        raise ValueError("Please provide Order or Superfamily")

    for key, val in data_reverse.items():
        ax1.plot(
            test_windows,
            val,
            linestyle=(0, (3, 10, 1, 10)),
            marker="o",
            label=str(key),
        )
    ax1.axis([max(test_windows), min(test_windows), 0, 1])
    ax1.set(
        xlabel="BP Upstream",
        ylabel="TE Density",
        yticks=np.arange(0, 1.01, 0.1),
        xticks=range(min(test_windows), (max(test_windows) + 1), 500),
    )

    for key, val in data_intra.items():
        ax2.scatter(
            0, val, label=str(key),
        )

    ax2.set(xlabel="Intronic TEs", xticks=[])
    ax2.legend(loc="upper right")

    for key, val in data_forward.items():
        ax3.plot(
            test_windows,
            val,
            linestyle=(0, (3, 10, 1, 10)),
            marker="o",
            label=str(key),
        )
    ax3.axis([min(test_windows), max(test_windows), 0, 1])
    ax3.set(
        xlabel="BP Downstream",
        xticks=range(min(test_windows), (max(test_windows) + 1), 500),
    )

    plt.savefig(os.path.join(output_dir, str(te_type + "_Density_Plot.png")))

    if display:
        plt.show()


if __name__ == "__main__":
    """Command line interface to graph."""
    path_main = os.path.abspath(__file__)
    parser = argparse.ArgumentParser(description="graph TE density")

    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=os.path.join(path_main, "../../..", "TE_Data/graphs"),
        help="parent directory to output results",
    )
    args = parser.parse_args()
    args.output_dir = os.path.abspath(args.output_dir)

    # Give it a dictionary for upstream, downstream, and intronic TEs, where
    # that dictionary represents the Order or Superfamily as well.
    # Should give it an upstream and downstream dictionary
    plot_density_all("Order", args.output_dir, display=True)
