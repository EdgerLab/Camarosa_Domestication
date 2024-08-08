#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs
from src.concatenate_leaf_expression_data import read_filtered_expression_data


import matplotlib.pyplot as plt
import matplotlib

# import matplotlib.axes.Axes.twiny as twiny

"""
- Import Pat A, B, C, and D syntelog table
"""


def import_unclean_syntelog_data(syntelog_input_file):
    """
    Import and filter Pat's Excel table of syntelogs for the A, B, C, and D
    homeologs for Royal Royce
    """
    data = pd.read_csv(
        syntelog_input_file,
        header=0,
        names=["Gene_A_Name", "Gene_B_Name", "Gene_C_Name", "Gene_D_Name", "Set"],
        dtype={
            "Gene_A_Name": str,
            "Gene_B_Name": str,
            "Gene_C_Name": str,
            "Gene_D_Name": str,
            "Set": str,
        },
        sep="\t",
        na_values=["---"],
    )

    # Let's first only look at rows where we have a syntelog in all four sets
    # (quartet)
    # data = data.loc[data["Set"] == "ABCD"]
    return data


def merge_expression_data(df, density_and_expression, gene_column):
    return (
        df.merge(
            density_and_expression,
            left_on=gene_column,
            right_on="Gene_Name",
            how="left",
        )
        .rename(
            columns={
                "Avg_Expression": f"{gene_column}_Avg_Expression",
                "Density": f"{gene_column}_Density",
            }
        )
        .drop(columns=["Gene_Name"])
    )


# Find the gene with the highest expression and density in each row
def find_highest_expression_and_density(row):
    genes = ["Gene_A_Name", "Gene_B_Name", "Gene_C_Name", "Gene_D_Name"]
    highest_expression_gene = max(
        genes,
        key=lambda x: row[f"{x}_Avg_Expression"]
        if pd.notna(row[f"{x}_Avg_Expression"])
        else -float("inf"),
    )
    highest_density_gene = max(
        genes,
        key=lambda x: row[f"{x}_Density"]
        if pd.notna(row[f"{x}_Density"])
        else -float("inf"),
    )
    return pd.Series(
        {
            "Highest_Expression_Gene": highest_expression_gene,
            "Highest_Expression_Value": row[
                f"{highest_expression_gene}_Avg_Expression"
            ],
            "Highest_Density_Gene": highest_density_gene,
            "Highest_Density_Value": row[f"{highest_density_gene}_Density"],
        }
    )


if __name__ == "__main__":
    # Parse command line arguments
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(
        description="""Import unclean gene expression data previously
        downloaded from Edger Lab"""
    )

    parser.add_argument(
        "syntelog_input_file",
        type=str,
        help="Path to Pat Excel A,B,C,D RR syntelog data",
    )
    parser.add_argument("expression_data", type=str)
    parser.add_argument("density_data", type=str)
    parser.add_argument(
        "output_path",
        type=str,
        help="Path to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.syntelog_input_file = os.path.abspath(args.syntelog_input_file)
    args.expression_data = os.path.abspath(args.expression_data)
    args.density_data = os.path.abspath(args.density_data)
    args.output_path = os.path.abspath(args.output_path)

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # ------------------------------------------------------------------------------------------------

    # Import and filter unclean gene expression data from Royal Royce
    syntelog_data = import_unclean_syntelog_data(args.syntelog_input_file)
    expression_data = read_filtered_expression_data(args.expression_data)

    density_data = pd.read_csv(args.density_data, sep="\t", header="infer")

    # Remove the reps from the expression data, we are just going to use the
    # average
    expression_data.drop(columns=["TPKM_Rep1", "TPKM_Rep2", "TPKM_Rep3"], inplace=True)

    # Lets keep the syntelogs only in quartet groups, and remove that column so
    # that things aren't complicated.
    # syntelog_data = syntelog_data.loc[syntelog_data["Set"] == "ABCD"]
    syntelog_data.drop(columns=["Set"], inplace=True)

    # Let's rename the density value to density
    density_data.rename(
        columns={"RR_Total_TE_Density_5000_Upstream": "Density"}, inplace=True
    )
    density_data.rename(columns={"RR_Gene": "Gene_Name"}, inplace=True)

    # Remove the other columns that we don't need
    density_data.drop(
        columns=["RR_Start", "RR_Stop", "RR_Strand", "RR_Length", "RR_Chromosome"],
        inplace=True,
    )

    # Merge the expression and density data
    density_and_expression = expression_data.merge(
        density_data, on="Gene_Name", how="left"
    )

    # Generate a bin plot of the density and expression data
    density_and_expression.set_index("Gene_Name", inplace=True)
    bins = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    density_and_expression["Density_Bin"] = pd.cut(
        density_and_expression["Density"], bins=bins, include_lowest=True
    )
    grouped = (
        density_and_expression.groupby("Density_Bin")
        .agg(Gene_Count=("Density", "size"), Avg_Expression=("Avg_Expression", "mean"))
        .reset_index()
    )

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.bar(
        grouped["Density_Bin"].astype(str),
        grouped["Gene_Count"],
        color="lightblue",
        label="Gene Count",
    )
    ax2 = ax1.twinx()
    ax2.plot(grouped["Density_Bin"].astype(str), grouped["Avg_Expression"], color="red")
    ax2.set_yscale("log")
    plt.show()

    # Check out the homeologs for each quartet, and checkout which genes have
    # high expression and high TE density

    # NOTE, may be worth forcing it so that I am only considering A, B, C, D
    # retained homeologs
    merged_A = merge_expression_data(
        syntelog_data, density_and_expression, "Gene_A_Name"
    )
    merged_B = merge_expression_data(merged_A, density_and_expression, "Gene_B_Name")
    merged_C = merge_expression_data(merged_B, density_and_expression, "Gene_C_Name")
    merged_D = merge_expression_data(merged_C, density_and_expression, "Gene_D_Name")

    result = merged_D.apply(find_highest_expression_and_density, axis=1)

    final_result = pd.concat([merged_D, result], axis=1)

    desired_output = final_result.loc[
        (final_result["Highest_Expression_Gene"] == "Gene_B_Name")
        & (final_result["Highest_Density_Gene"] == "Gene_B_Name")
        & (final_result["Highest_Density_Value"] > 0.50)
    ]

    # Pat also wants an example of 3 genes with high TE load, low expression,
    # and one gene has no TE and high expression.

    # NOTE, the correlation between expression and TE density is very close to
    # 0 which is odd.
    # print(
    #     final_result["Gene_A_Name_Avg_Expression"].corr(
    #         final_result["Gene_A_Name_Density"]
    #     )
    # )
    # print(
    #     final_result["Gene_D_Name_Avg_Expression"].corr(
    #         final_result["Gene_D_Name_Density"]
    #     )
    # )
    # raise ValueError

    desired_output.to_csv("test.tsv", header=True, sep="\t", index=False)

    # This one is probably the best because the other homeologs have (similar)
    # very low levels of TE presence, and similar expression
    # NOTE Fxa4Bg103345 has much higher expression than all other homeologs, and
    # has a CACTA with high identity right in front of it. This is a good example.

    # NOTe Fxa5Bg100284 is a good example as it has crazy expression and has a big
    # piece of LTR, LINE, and some other upstream. Other homeologs have good
    # expression, 2 of the homeologs have 0 TE, and one has 60\% TE but has lowest
    # expression.

    # NOTE Fxa5Bg103293 has crazy expression and has a SINE, Mutator and LTR in
    # front of it. Only one other homelog has 1/4th of the expression and 0.25 TE
    # Density in front, the other 2 basically have 0 expression and no TE.

    # Fxa6Bg101326 has high expression and compares well to its homeologs

    # NOTE now working on generating a simple bar plot of the expression and TE
    # density, each gene will have a bar for expression and a bar for TE density,
    # 4 pairs of bars will be plotted, one pair for each gene.
    # I will make a general bar plot of all the genes, and then a plot for my
    # biased set, and then maybe a single plot for a few genes of interest,
    # such as the ones described above.

    # Could do a groupby on this and make a general graph for everything
    # I want it all in one plot

    print(desired_output)

    # NOTE generate a bar graph for a single gene of interest, previously
    # manually identified.
    # Ok so I got a sample bar graph for ONE gene, how can I make this for all?
    gene_of_interest = "Fxa4Bg103345"
    Fxa5Bg103293 = desired_output.loc[desired_output["Gene_B_Name"] == gene_of_interest]
    groups = ["Homeolog_A", "Homeolog_B", "Homeolog_C", "Homeolog_D"]
    # BIG UGLY
    all_densities = [
        Fxa5Bg103293["Gene_A_Name_Density"].squeeze(),
        Fxa5Bg103293["Gene_B_Name_Density"].squeeze(),
        Fxa5Bg103293["Gene_C_Name_Density"].squeeze(),
        Fxa5Bg103293["Gene_D_Name_Density"].squeeze(),
    ]
    # BIG UGLY
    all_expression = [
        Fxa5Bg103293["Gene_A_Name_Avg_Expression"].squeeze(),
        Fxa5Bg103293["Gene_B_Name_Avg_Expression"].squeeze(),
        Fxa5Bg103293["Gene_C_Name_Avg_Expression"].squeeze(),
        Fxa5Bg103293["Gene_D_Name_Avg_Expression"].squeeze(),
    ]
    densities_and_expression = {"Density": all_densities, "Expression": all_expression}

    fig, ax1 = plt.subplots()
    x = np.arange(len(groups))  # the label locations
    width = 0.2  # MAGIC, to make the graph look nice, width of bars
    multiplier = 0  # MAGIC, to make the graph look nice
    for attribute_name, measurement in densities_and_expression.items():
        offset = width * multiplier
        if attribute_name == "Density":
            color = "tab:blue"
            rects = ax1.bar(
                x + offset, measurement, width, label=attribute_name, color=color
            )
            # ax1.bar_label(rects, padding=3)
            ax1.set_ylabel("Density", color=color)
            ax1.set_xticks(x + width, groups)
            ax1.set_ylim(0, 1.0)
        if attribute_name == "Expression":
            color = "tab:red"
            ax2 = ax1.twinx()  # share the x axis
            rects = ax2.bar(
                x + offset, measurement, width, label=attribute_name, color=color
            )
            # ax2.set_yscale("log")
            # ax2.bar_label(rects, padding=3)
            ax2.set_ylabel("Expression", color=color)
        multiplier += 1
        # ax1.tick_params(axis="y", labelcolor=color)
        # ax2.tick_params(axis="y", labelcolor=color)
    plt.title(gene_of_interest)
    plt.show()
