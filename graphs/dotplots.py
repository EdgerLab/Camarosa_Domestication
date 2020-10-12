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


def plot_density_outer(te_type, direction, output_dir, display=False):
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

    """
    # Get the mean density value for all genes for each TE type. E.g LTR
    # density at 500 distance, take the mean of that for all genes, that will
    # be the value plotted. Do this for Order and SuperFamily

    # TODO:
    # Funciton to check the density values being provided to graph are
    # the same length as the number of windows.

    # Should provide it a dictionary of lists
    test_densities_ltr = [0.3, 0.6, 0.5, 0.2]
    test_densities_line = [0.8, 0.9, 0.2, 0.3]
    test_windows = [500, 1000, 1500, 2000]

    data = {"ltr": test_densities_ltr, "line": test_densities_line}

    for key, val in data.items():
        plt.plot(
            test_windows,
            val,
            linestyle=(0, (3, 10, 1, 10)),
            marker="o",
            label=str(key),
        )
        plt.legend()

    # NOTE how to do it without iterating over a dictionary
    # plt.plot("windows", "ltr", linestyle=(0, (3, 10, 1, 10)), marker="o", data=data)

    if direction == "Upstream":
        plt.axis([max(test_windows), min(test_windows), 0, 1])
    elif direction == "Downstream":
        plt.axis([min(test_windows), max(test_windows), 0, 1])
    else:
        raise ValueError("Upstream or Downstream")
    plt.xticks(range(min(test_windows), (max(test_windows) + 1), 500))
    plt.xlabel("Window Position in BP")
    plt.ylabel("Density Value")

    if te_type == "Order":
        plt.title("TE Orders")
    elif te_type == "Superfamily":
        plt.title("TE Superfamilies")
    else:
        raise ValueError("Please provide Order or Superfamily")
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

    # Should try to provide a dictionary of lists, with the key being the TE,
    # and the list containing all of the density values for the windows.
    # Pre-filter ahead of time?
    # plot_density_intra("Order", args.output_dir, display=True)

    # Should give it an upstream and downstream dictionary
    plot_density_outer("Order", "Upstream", args.output_dir, display=True)
