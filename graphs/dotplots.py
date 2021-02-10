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
from collections import defaultdict


def plot_intra_density(dd_obj, order_or_super, output_dir, display=False):
    """
    Args:
        dd_obj (DensityData):

        order_or_super (str): String representing the category of TE, should be either
        'Order' or 'Superfamily', raises ValueError if it does not match.

        output_dir (str): Path representing output directory.

        display (str): Defaults to False, if True it will show the graph
        (usually opens another window).

    """
    plotting_dict = {}

    if order_or_super == "Order":
        te_index_dict = dd_obj.order_index_dict
        h5_frame = dd_obj.intra_orders
    elif order_or_super == "Superfamily":
        te_index_dict = dd_obj.super_index_dict
        h5_frame = dd_obj.intra_supers
    else:
        raise ValueError("Please provide Order or Superfamily")
    # NOTE
    # I could add a FOR loop here to go over each dd_obj and get the mean that
    # way, plotting that would be super easy and would only yield one line.
    for te_type, index_val in te_index_dict.items():
        plotting_dict[te_type] = np.mean(h5_frame[te_index_dict[te_type], :, :])

    for key, val in plotting_dict.items():
        plt.scatter(0, val, label=key)
        plt.legend()

    plt.yticks(np.arange(0, 1.01, 0.1))
    plt.xticks([0])
    plt.xlabel("Window Position in BP")
    plt.ylabel("Intronic TE Density")

    if order_or_super == "Order":
        plt.title("TE Orders")
    elif order_or_super == "Superfamily":
        plt.title("TE Superfamilies")
    else:
        raise ValueError("Please provide Order or Superfamily")

    plt.savefig(
        os.path.join(
            output_dir, str(order_or_super + "_" + dd_obj.genome_id + "_Intra_Plot.png")
        )
    )
    if display:
        plt.show()
    plt.close()


def plot_density_all(dd_obj, order_or_super, output_dir, display=False):
    """
    Plot the mean density value of each TE category over each window for "left"
    and "right".

    Args:
        dd_obj (DensityData):

        order_or_super (str): String representing the category of TE, should be either
        'Order' or 'Superfamily', raises ValueError if it does not match.

        output_dir (str): Path representing output directory.

        display (str): Defaults to False, if True it will show the graph
        (usually opens another window).

    """
    if order_or_super == "Order":
        te_index_dict = dd_obj.order_index_dict
        left_h5_frame = dd_obj.left_orders
        intra_h5_frame = dd_obj.intra_orders
        right_h5_frame = dd_obj.right_orders
    elif order_or_super == "Superfamily":
        te_index_dict = dd_obj.super_index_dict
        left_h5_frame = dd_obj.left_supers
        intra_h5_frame = dd_obj.intra_supers
        right_h5_frame = dd_obj.right_supers
    else:
        raise ValueError("Please provide Order or Superfamily")

    # NOTE rename this to upstream and downstream
    data_left_dict = defaultdict(list)
    data_intra_dict = defaultdict(list)
    data_right_dict = defaultdict(list)
    for te_type, index_val in te_index_dict.items():
        for window_idx in range(len(dd_obj.window_list)):
            data_left_dict[te_type].append(
                np.mean(left_h5_frame[te_index_dict[te_type], window_idx, :])
            )
            data_right_dict[te_type].append(
                np.mean(right_h5_frame[te_index_dict[te_type], window_idx, :])
            )
        data_intra_dict[te_type].append(
            np.mean(intra_h5_frame[te_index_dict[te_type], :, :])
        )

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey="col")
    fig.set_size_inches(16, 8.5)
    # define colors
    NUM_COLORS = sum(1 for te_type in te_index_dict.items())
    cm = plt.get_cmap("tab20")
    ax1.set_prop_cycle("color", [cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])
    ax2.set_prop_cycle("color", [cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])
    ax3.set_prop_cycle("color", [cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])

    for key, val in data_left_dict.items():
        ax1.plot(
            dd_obj.window_list, val, label=key, linestyle=(0, (3, 1, 1, 1)), marker="o",
        )
    ax1.axis([max(dd_obj.window_list), min(dd_obj.window_list), 0, 1])
    ax1.set(
        xlabel="BP Upstream",
        ylabel="TE Density",
        yticks=np.arange(0, 1.01, 0.1),
        xticks=range(min(dd_obj.window_list), (max(dd_obj.window_list) + 1), 1000),
    )

    for key, val in data_intra_dict.items():
        ax2.scatter(
            0, val, label=key,
        )

    ax2.set(xlabel="Intronic TEs", xticks=[])
    ax2.legend(loc="upper right")

    for key, val in data_right_dict.items():
        ax3.plot(
            dd_obj.window_list, val, label=key, linestyle=(0, (3, 1, 1, 1)), marker="o",
        )
    ax3.axis([min(dd_obj.window_list), max(dd_obj.window_list), 0, 1])
    ax3.set(
        xlabel="BP Downstream",
        yticks=np.arange(0, 1.01, 0.1),
        xticks=range(min(dd_obj.window_list), (max(dd_obj.window_list) + 1), 1000),
    )
    ax3.yaxis.tick_right()

    # for key, val in data_right_dict.items():
    # plt.scatter(dd_obj.window_list, val, label=key)
    # plt.legend()
    plt.savefig(
        os.path.join(
            output_dir,
            str(order_or_super + "_" + dd_obj.genome_id + "_Combined_Density_Plot.png"),
        ),
        bbox_inches="tight",
    )
    if display:
        plt.show()
