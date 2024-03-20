import pandas as pd
import numpy as np
import os
import argparse
import logging
import coloredlogs

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

import scipy.stats as stats
import statsmodels.api as sm

"""
Generate barplots of syntelog TE density differences
"""


def graph_barplot_density_differences(
    values,
    te_type,
    window_val,
    direction,
    number_of_zeros,
    output_dir,
    logger,
    display=False,
    align="left",
    upper_and_lower_cutoff_int=(95, 5),
):
    """
    Plot a histogram of TE density differences between syntelog pairs

    Args:
        values (list): A list of values representing the TE density differences
        between syntelog pairs

        te_type (str): String representing the TE type being plotted

        window_val (int): Integer representing the current window of which the
        data is being plotted

        direction (str): string representing whether or not the graphs are
        coming from upstream or downstream TE density data

        number_of_zeros ():

        logger (logging.Logger): Object to log information to

        display (boolean): Defaults to False, if True shows the plot upon
        generation with the plt.show() command
    """
    # MAGIC, the bins to group density values for the histogram AND the values
    # for the xticks on the xaxis
    tick_bins = [
        -1.0,
        -0.9,
        -0.8,
        -0.7,
        -0.6,
        -0.5,
        -0.4,
        -0.3,
        -0.2,
        -0.1,
        0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
    ]
    plt.figure(figsize=(8, 6))
    n, bins, patches = plt.hist(
        values, bins=tick_bins, facecolor="blue", ec="black", alpha=0.5
    )
    plt.rcParams["xtick.labelsize"] = 7  # MAGIC set size of axis ticks
    plt.ylabel("Number of Genes")
    plt.xlabel("Difference in TE Density Values")
    plt.title("Del Norte vs Royal Royce")  # MAGIC genome name order here
    N = mpatches.Patch(
        label="Total Plotted Genes: %s \nTE type: %s \nWindow: %s \nDirection: %s \nNo. 0 Differences: %s"
        % (len(values), te_type, window_val, direction, str(number_of_zeros))
    )
    # -----------------------------------------
    # Cutoffs
    upper_cutoff_val = np.percentile(values, upper_and_lower_cutoff_int[0])
    lower_cutoff_val = np.percentile(values, upper_and_lower_cutoff_int[1])
    plt.axvline(upper_cutoff_val, color="g", linestyle="dashed", linewidth=2)
    plt.axvline(lower_cutoff_val, color="g", linestyle="dashed", linewidth=2)

    upper_label = plt.Line2D(
        [],
        [],
        color="g",
        marker="",
        linestyle="",
        label=f"Upper {upper_and_lower_cutoff_int[0]} Percentile Cutoff: {upper_cutoff_val:.2f}",
    )
    lower_label = plt.Line2D(
        [],
        [],
        color="g",
        marker="",
        linestyle="",
        label=f"Bottom {upper_and_lower_cutoff_int[1]} Percentile Cutoff: {lower_cutoff_val:.2f}",
    )
    # -----------------------------------------
    # Wilcoxon signed-rank test
    # TODO possibly add this to the legend or print out
    # Wilxocon signed rank test is supposed be to carried out after the zero
    # difference values are filtered away
    # Because the data does not come from a normal distribution, I will use the
    # Wilcoxon signed rank test to test if the median of the differences is
    # significantly different than zero.
    # wilcoxon = stats.wilcoxon(values, alternative="two-sided")

    # I could also edit the alernative hypothesis to be one-sided, and specify
    # that DN is always higher than RR. This would help me get at whether or
    # not DN is more dense.
    wilcoxon = stats.wilcoxon(values, alternative="greater")

    significance = wilcoxon[1]
    significance = plt.Line2D(
        [],
        [],
        color="r",
        marker="",
        linestyle="",
        label=f"Wilcoxon signed-rank test p-value: {significance:.3e}",
    )

    # -----------------------------------------
    # Get the mean and plot as red line, create object for legend
    diff_mean = np.mean(values)
    mean_label = plt.Line2D(
        [], [], color="r", marker="", linestyle="", label=f"Mean: {diff_mean:.3f}"
    )
    plt.axvline(diff_mean, color="r", linestyle="dashed", linewidth=2)

    plt.xticks(tick_bins)
    plt.legend(
        handles=[N, mean_label, lower_label, upper_label, significance],
        loc="upper left",
    )
    path = os.path.join(
        output_dir,
        (te_type + "_" + str(window_val) + "_" + direction + "_DensityDifferences.png"),
    )

    logger.info("Saving graph to: %s" % path)
    plt.savefig(path)
    if display:
        plt.show()
    plt.close()


def decode_te_window_direction_str(filename):
    """
    Decode the TE type, window, and direction from a string, the filename of a
    preprocessed density table
    """
    # NOTE a lot of MAGIC here, the filename is expected to be in a certain
    # pattern
    if len(filename.split("_")) < 7:
        te_type, window = filename.split("_")[3:5]
        direction = filename.split("_")[5].split(".")[0]
    else:
        te_type = filename.split("_")[3:6]
        te_type = "_".join(te_type)
        window = filename.split("_")[6]
        direction = filename.split("_")[7].split(".")[0]
    return te_type, window, direction


def calculate_number_of_nonzero_differences(table):
    return len(get_table_of_nonzero_differences(table))


def get_table_of_nonzero_differences(table):
    return table.loc[table["Difference"] != 0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "preprocessed_density_table",
        type=str,
        help="TODO",
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Parent directory to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    logger = logging.getLogger(__name__)
    args = parser.parse_args()
    args.preprocessed_density_table = os.path.abspath(args.preprocessed_density_table)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    table = pd.read_csv(args.preprocessed_density_table, header="infer", sep="\t")
    te_type, window, direction = decode_te_window_direction_str(
        os.path.basename(args.preprocessed_density_table)
    )

    num_original_rows = len(table)
    num_non_zero = calculate_number_of_nonzero_differences(table)
    num_zero = num_original_rows - num_non_zero
    table_sans_zeros = get_table_of_nonzero_differences(table)

    # This is the code for a qqplot which tells us how well the data fits a
    # normal distribution
    # Based on the qq plot, the data does not appear to come from a normal
    # distribution
    # fig = sm.qqplot(table["Difference"], line="45")
    # plt.show()
    # plt.clf()

    # This is code to perform a Shapiro Wilkes test, more formal way to see if
    # normal distribution
    # The data does not appear to come from a normal distribution
    # Shapiro = stats.shapiro(table["Difference"])

    graph_barplot_density_differences(
        table_sans_zeros["Difference"],
        te_type,
        window,
        direction,
        num_zero,
        args.output_dir,
        logger,
        display=False,
    )
