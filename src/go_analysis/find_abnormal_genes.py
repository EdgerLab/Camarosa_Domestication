#!/usr/bin/env/python


"""
NOTE this script does too much and is a candidate for refactoring

Objectives:
1. Given a table of pre-filtered TE Density of strawberry genes, this script
generates more tables by calculating the top X and bottom Y percentiles of data
and saving them to disk. The plan is to use these tables to run TopGO and
determine which genes are over or under represented in the datasets.

An example of an output table is the top 1% most LTR-dense genes in the RR
genome for the 5KB upstream window.

2. This script also calculates the top and bottom percentiles of the 'Difference'
column which is was previously calculated and is the Royal Royce (RR) TE Density value
substracted from the Del Norte (DN) TE Density value.

3. This script also calculates the number of genes that have strawberry and AT
orthologs. This is done for a TE-dense gene set, and then the same number of
genes is randomly sampled from the genome N (usually 1000) times and the number
of AT orthologs is averaged. This is compared to the TE-dense gene set to see
if TE-dense genes have a reduced number of AT orthologs. A barplot and
histogram is generated for this purpose
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd
from collections import namedtuple

import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.api as sm
import matplotlib.patches as mpatches

from src.syntelog_differences.bargraphs import decode_te_window_direction_str
from src.extract_AED_score import read_aed_output_table


def subset_by_aed_score(mod_table, aed_scores, genome, good_aed_score):
    # Merge mod_table with the AED scores
    mod_table = mod_table.merge(
        aed_scores,
        on=f"{genome}_Gene",
        how="left",
    )
    # Subset mod_table by the good AED scores
    mod_table = mod_table.loc[mod_table["AED_Score"] <= good_aed_score]
    return mod_table


def plot_count_of_remaining_genes(
    count_of_genes_meeting_cutoff,
    count_of_genes_w_strawberry_orthologs,
    count_of_genes_w_arabidopsis_orthologs,
    count_of_reference_genes_meeting_cutoff,
    count_of_reference_w_strawberry,
    count_of_reference_w_arabidopsis,
    logger,
    output_filename,
    te_col,
    aed_scores,
):
    x = 0
    group_labels = ["Randomly Sampled Genes", "TE-Dense Genes"]

    bar_labels_and_vals = {
        "Starting Genes": [
            count_of_reference_genes_meeting_cutoff,
            count_of_genes_meeting_cutoff,
        ],
        "Count w/ Strawberry Orthologs": [
            count_of_reference_w_strawberry,
            count_of_genes_w_strawberry_orthologs,
        ],
        "Count w/ Arabidopsis Orthologs": [
            count_of_reference_w_arabidopsis,
            count_of_genes_w_arabidopsis_orthologs,
        ],
    }

    obsv_percent_reduction_to_arabidopsis = (
        count_of_genes_w_arabidopsis_orthologs / count_of_genes_meeting_cutoff
    )
    ref_percent_reduction_to_arabidopsis = (
        count_of_reference_genes_meeting_cutoff / count_of_reference_w_arabidopsis
    )

    # Get the X label locations
    x = np.arange(len(group_labels))
    width = 0.2  # MAGIC, to make the graph look nice
    multiplier = 0  # MAGIC, to make the graph look nice

    fig, ax = plt.subplots(layout="constrained", figsize=(10, 6))
    for attribute, measurement in bar_labels_and_vals.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute)

        # Displaying percentage of each value divided by 100 above each bar
        percentage_values = [
            ((value / count_of_genes_meeting_cutoff) * 100) for value in measurement
        ]
        ax.bar_label(
            rects,
            labels=[
                f"{base_val:.2f}\n({percent_val:.2f})%"
                for percent_val, base_val in zip(percentage_values, measurement)
            ],
            padding=3,
        )
        multiplier += 1

    ax.set_xticks(x + width, group_labels)
    ax.legend(loc="upper right", ncols=1, fontsize="small")
    ax.set_ylabel("Count of Genes")
    ax.set_xlabel("Groupings")

    # MAGIC set ylim to be a little taller so we can display the raw number and
    # percentage above the bars
    ax.set_ylim(
        0, count_of_genes_meeting_cutoff + (0.10 * count_of_genes_meeting_cutoff)
    )
    te_col = te_col.replace("_", " ")  # string format for plot title
    ax.set_title(
        f"""Difference in Ortholog Count Between Randomly Sampled and TE-Dense Genes: \n{te_col}"""
    )
    logger.info(f"Saving graph to: {output_filename}")
    plt.savefig(output_filename, bbox_inches="tight")
    plt.close()
    return


def generate_random_sample(
    base_table, ortholog_table, genome, aed_scores, good_aed_score, n
):
    """
    Randomly sample the base table (the table of ALL genes and their
    TE Density values for a specific window, TE type, direction dimension) N
    times and return the count of genes that have an ortholog in the random
    sample.

    Returns a list because we are randomly sampling and doing this 'n' times so
    we want to be able to average later

    Returns 2 things:
        random_strawberry (list of ints): count of genes with strawberry
            orthologs
        reference_arabidopsis_counts (list of ints): count of genes wit
            Arabidopsis orthologs
    """
    # Generate a comparable dataset from random sampling
    random_strawberry = []
    reference_arabidopsis_counts = []
    for i in range(n):
        random_sample = base_table.sample(count_of_genes_meeting_cutoff)
        random_merged = random_sample.merge(
            ortholog_table,
            left_on=[f"{genome}_Gene", f"{genome}_Chromosome"],
            right_on=[f"{genome}_Gene", f"{genome}_Chromosome"],
            how="inner",
        )

        random_merged = subset_by_aed_score(
            random_merged, aed_scores, genome, good_aed_score
        )

        # TODO write a loc command to make sure that the random strawberry is
        # checking for the OTHER genome to have a gene, not the same genome. I
        # think currently it works because the ortholog table is forced to have
        # something for DN <-> RR...
        random_strawberry.append(len(random_merged))
        reference_arabidopsis_counts.append(
            len(random_merged.loc[~random_merged["Arabidopsis_Gene"].isna()])
        )
    return random_strawberry, reference_arabidopsis_counts


def plot_random_distribution_vs_observed(
    reference_arabidopsis_counts,
    count_of_genes_w_arabidopsis_orthologs,
    std_dev_away,
    logger,
    output_filename,
    te_col,
):
    """
    Plot a histogram of the distribution of the random sample and surviving AT
    ortholog count, and plot the observed value (TE dense genes and the count
    of AT orthologs) as a red line.
    """
    muh_bins = 20  # MAGIC bincount to make the histogram look nice
    n, bins, hist_patch = plt.hist(
        reference_arabidopsis_counts,
        bins=muh_bins,
        label="Distribution of Randomly Sampled Genes w/ AT Orthologs",
        edgecolor="black",
        linewidth=0.4,
    )
    observation = plt.axvline(
        x=count_of_genes_w_arabidopsis_orthologs,
        color="red",
        linewidth=2,
        label=f"Count of TE-Dense Genes w/ AT Orthologs \nValue is {std_dev_away:.3f} Standard Deviations Away from the Mean",
    )
    plt.xlabel("Count of Genes with AT Orthologs")
    plt.ylabel("Frequency of Genes with AT Orthologs in Bin")
    plt.legend(loc="upper right")
    # Plot title is too big, may exclude because it is not needed for the paper
    # plt.title(
    #     "Distribution of Randomly Sampled Genes w/ AT Orthologs vs Count of AT Orthologs in TE-Dense Genes"
    # )

    te_col = te_col.replace("_", " ")  # string format for plot title
    plt.title(te_col)
    logger.info(f"Saving graph to: {output_filename}")
    plt.savefig(output_filename, bbox_inches="tight")
    plt.close()


def calculate_cutoff_value(dataframe, column_name, cutoff_int):
    cutoff_float = cutoff_int / 100
    return dataframe[column_name].quantile(cutoff_float)


def perform_upper_cutoff_and_subset(
    table, column_name, upper_percentile_cutoff_int, logger
):
    """
    TODO
    """
    table = table.copy(deep=True)

    # Make sure nonzero values are used?
    # TODO verify this is the right thing to do, do we want to remove non-zeros
    # before we calculate the cutoff?
    # But most genes have a 0 value for TE density, so it makes sense to test
    # those who do have a TE near them
    # table = table.loc[table[column_name] != 0.0]

    upper_cutoff_val = calculate_cutoff_value(
        table, column_name, upper_percentile_cutoff_int
    )

    # Check the cutoff value and raise a warning if the value is below 10%
    if upper_cutoff_val < 0.1:
        logger.warning(
            f"Upper Cutoff Value is below 10% density: it is {upper_cutoff_val} for {column_name}"
        )

    # This column is used to store the cutoff value for each gene and may
    # not actually be needed
    table["Upper_Cutoff_Value"] = upper_cutoff_val

    # Subset the table by the cutoff value
    table = table.loc[table[column_name] >= upper_cutoff_val]
    return table


def perform_lower_cutoff_and_subset(
    table, column_name, lower_percentile_cutoff_int, logger
):
    """
    This is only needed for the difference column?
    """
    table = table.copy(deep=True)

    # TODO verify if this is the approach that we want
    # Make sure nonzero values are used
    # table = table.loc[table[column_name] != 0.0]

    lower_cutoff_val = calculate_cutoff_value(
        table, column_name, lower_percentile_cutoff_int
    )

    # Check the cutoff value and raise a warning if the value is below 10%
    # TODO revist if this warning is needed for the lower cutoff
    if lower_cutoff_val > -0.1:
        logger.warning(
            f"Lower Cutoff Value is below 10%: {lower_cutoff_val} for {column_name}"
        )

    # This column is used to store the cutoff value for each gene and may
    # not actually be needed
    table["Lower_Cutoff_Value"] = lower_cutoff_val

    # Subset the table by the cutoff value
    # NOTE this is like the only thing different from the upper cutoff
    # function... TODO think about how to refactor
    table = table.loc[table[column_name] <= lower_cutoff_val]
    return table


def subset_by_arabidopsis_presence(table):
    """
    Return a table where each entry has an Arabidopsis gene
    """
    return table.loc[~table["Arabidopsis_Gene"].isna()]


def save_table_to_disk(table, output_dir, out_filename, logger):
    out_filepath = os.path.join(output_dir, out_filename)
    logger.info(f"Writing to: {out_filepath}")
    table.to_csv(out_filepath, sep="\t", header=True, index=False)


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "preprocessed_density_table",
        type=str,
        help="TODO",
    )
    parser.add_argument("aed_score_table", type=str, help="TODO")

    parser.add_argument(
        "upper_percentile_cutoff_int",
        type=int,
        help="TODO",
    )
    parser.add_argument(
        "lower_percentile_cutoff_int",
        type=int,
        help="TODO",
    )

    parser.add_argument("strawberry_ortholog_table", type=str, help="TODO")

    parser.add_argument(
        "output_dir",
        type=str,
        help="parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.preprocessed_density_table = os.path.abspath(args.preprocessed_density_table)
    args.strawberry_ortholog_table = os.path.abspath(args.strawberry_ortholog_table)
    args.aed_score_table = os.path.abspath(args.aed_score_table)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -------------------------------------------------------------------------

    # Load in the pre-filtered TE Density data
    base_table = pd.read_csv(args.preprocessed_density_table, header="infer", sep="\t")

    # Load in the pre-filtered ortholog information
    ortholog_table = pd.read_csv(
        args.strawberry_ortholog_table, header="infer", sep="\t", low_memory=False
    )

    # Magic remove the file extension after basename, won't work if multiple
    # '.' in the filename
    filename = os.path.splitext(os.path.basename(args.preprocessed_density_table))[0]
    if "minus" not in filename:
        if "Total" in filename:
            genome = filename.split("_")[0]
            te_type = filename.split("_")[1:4]
            te_type = "_".join(te_type)
            window = filename.split("_")[4]
            direction = filename.split("_")[5]
        else:
            pieces = filename.split("_")
            genome = pieces[0]
            te_type = pieces[1]
            window = pieces[2]
            direction = pieces[3]
        flag = True  # this is baaaaad to do it this way...
    else:
        # We are working with the DN minus RR table
        te_type, window, direction = decode_te_window_direction_str(
            os.path.basename(args.preprocessed_density_table)
        )
        # TODO NOTE the AED score stuff doesn't work for the DN minus RR table
        # because I would need two AED score tables
        flag = False

    # Load in the AED scores
    aed_scores = read_aed_output_table(args.aed_score_table)
    aed_scores.drop(columns=["Chromosome"], inplace=True)
    aed_scores.rename(columns={"Gene_Name": f"{genome}_Gene"}, inplace=True)
    # MAGIC MAGIC MAGIC, anything below this is a well-supported gene model
    good_aed_score = 0.75

    # Create a named tuple that is the cutoff function and its percentile in
    # integer format, to avoid needing to retype this all the time
    cutoff_function_w_percentile = namedtuple(
        "cutoff_function_w_percentile",
        ["cutoff_function", "percentile", "upper_or_lower_str"],
    )
    upper = cutoff_function_w_percentile(
        perform_upper_cutoff_and_subset, args.upper_percentile_cutoff_int, "Upper"
    )
    lower = cutoff_function_w_percentile(
        perform_lower_cutoff_and_subset, args.lower_percentile_cutoff_int, "Lower"
    )
    if flag:
        # We are starting with the individual genome files, and only doing the
        # 'upper' function
        # Perform the cutoffs and subset, we are not working with the file that
        # has the 'Difference' column or the ortholog information
        te_col = f"{genome}_{te_type}_{window}_{direction}"
        mod_table = upper.cutoff_function(base_table, te_col, upper.percentile, logger)

        # Apply an AED score cutoff to the dense gene table
        mod_table = subset_by_aed_score(mod_table, aed_scores, genome, good_aed_score)

        # Save the original table before we do any subsetting with Arabidopsis
        # presence
        out_filename = f"{te_col}_{upper.upper_or_lower_str}_{str(upper.percentile)}_no_Arabidopsis_density_percentile.tsv"

        # MAGIC filepath matches with Makefile
        special_out_dir = os.path.join(args.output_dir, "no_Arabidopsis")
        save_table_to_disk(
            mod_table,
            special_out_dir,
            out_filename,
            logger,
        )
        count_of_genes_meeting_cutoff = len(mod_table)

        # Merge in the ortholog table because the individual genome files do
        # not have the ortholog information in them
        # Subset the table so that each entry MUST have an Arabidopsis gene
        mod_table = mod_table.merge(
            ortholog_table,
            left_on=[f"{genome}_Gene", f"{genome}_Chromosome"],
            right_on=[f"{genome}_Gene", f"{genome}_Chromosome"],
            how="inner",
        )
        count_of_genes_w_strawberry_orthologs = len(mod_table)
        count_of_genes_w_arabidopsis_orthologs = len(
            mod_table.loc[~mod_table["Arabidopsis_Gene"].isna()]
        )

        # Calculate the random distro barplots and histogram
        # NOTE this is big ugly, should put the below into a function
        if genome != "H4":

            # TODO parametrize this so we have a choice to not do 1000 random
            # samples
            N = 1000  # MAGIC, number of random samples to generate
            ref_surviving_straw, ref_surviving_AT = generate_random_sample(
                base_table, ortholog_table, genome, aed_scores, good_aed_score, N
            )
            # Define additional reference values
            reference_strawberry_mean = np.mean(ref_surviving_straw)
            reference_arabidopsis_mean = np.mean(ref_surviving_AT)
            reference_arabidopsis_std = np.std(ref_surviving_AT)

            logger.info(f"Working on {te_col}")
            logger.info(
                f"My reference standard deviation is {reference_arabidopsis_std}"
            )
            logger.info(f"My reference mean is {reference_arabidopsis_mean}")
            logger.info(
                f"My observed sample value is {count_of_genes_w_arabidopsis_orthologs}"
            )
            std_dev_away = (
                count_of_genes_w_arabidopsis_orthologs - reference_arabidopsis_mean
            ) / reference_arabidopsis_std

            output_filename = os.path.join(
                args.output_dir,
                "ortholog_analysis",
                f"{te_col}_{upper.upper_or_lower_str}_{str(upper.percentile)}_ortholog_histogram.png",
            )
            plot_random_distribution_vs_observed(
                ref_surviving_AT,
                count_of_genes_w_arabidopsis_orthologs,
                std_dev_away,
                logger,
                output_filename,
                te_col,
            )

            output_filename = os.path.join(
                args.output_dir,
                "ortholog_analysis",
                f"{te_col}_{upper.upper_or_lower_str}_{str(upper.percentile)}_ortholog_barplot.png",
            )
            plot_count_of_remaining_genes(
                count_of_genes_meeting_cutoff,
                count_of_genes_w_strawberry_orthologs,
                count_of_genes_w_arabidopsis_orthologs,
                count_of_genes_meeting_cutoff,
                reference_strawberry_mean,
                reference_arabidopsis_mean,
                logger,
                output_filename,
                te_col,
                aed_scores,
            )

        out_filename = f"{te_col}_{upper.upper_or_lower_str}_{str(upper.percentile)}_density_percentile.tsv"
        save_table_to_disk(
            mod_table,
            args.output_dir,
            out_filename,
            logger,
        )

    else:

        # ----------------------------------------
        # Calculate the 'Difference' column
        column_name = "Difference"
        col_to_save = f"{column_name}_{te_type}_{window}_{direction}"
        for genome, i in [("DN", upper), ("RR", lower)]:
            # We don't have a 'Difference' column for H4, the difference
            # column is for DN vs RR.
            # NOTE MAGIC, if positive it is a DEL NORTE biased gene pair
            # NOTE MAGIC, if negative it is a Royal Royce biased gene pair, hence
            # we will apply the lower cutoff function to the RR dataset

            # Perform the cutoffs and subset
            mod_table = i.cutoff_function(base_table, column_name, i.percentile, logger)

            # Subset the table so that each entry MUST have an Arabidopsis gene
            mod_table = subset_by_arabidopsis_presence(mod_table)
            out_filename = f"{col_to_save}_Biased_Towards_{genome}_{str(i.percentile)}_density_percentile.tsv"
            save_table_to_disk(
                mod_table,
                args.output_dir,
                out_filename,
                logger,
            )
