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


# NOTE hardcoded to only work with B genes at the moment
def plot_gene_of_interest_density_and_expression(dataframe, gene_of_interest):
    subsetted_data = dataframe.loc[dataframe["Gene_B_Name"] == gene_of_interest]
    groups = ["Homeolog_A", "Homeolog_B", "Homeolog_C", "Homeolog_D"]

    # BIG UGLY
    all_densities = [
        subsetted_data["Gene_A_Name_Density"].squeeze(),
        subsetted_data["Gene_B_Name_Density"].squeeze(),
        subsetted_data["Gene_C_Name_Density"].squeeze(),
        subsetted_data["Gene_D_Name_Density"].squeeze(),
    ]
    # BIG UGLY
    all_expression = [
        subsetted_data["Gene_A_Name_Avg_Expression"].squeeze(),
        subsetted_data["Gene_B_Name_Avg_Expression"].squeeze(),
        subsetted_data["Gene_C_Name_Avg_Expression"].squeeze(),
        subsetted_data["Gene_D_Name_Avg_Expression"].squeeze(),
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


def find_highest_expression_and_lowest_density(row):
    """
    Find the gene with the highest expression and lowest density in each row
    NOTE, semi duplicate code with find_highest_expression_and_highest_density
    """
    genes = ["Gene_A_Name", "Gene_B_Name", "Gene_C_Name", "Gene_D_Name"]
    highest_expression_gene = max(
        genes,
        key=lambda x: row[f"{x}_Avg_Expression"]
        if pd.notna(row[f"{x}_Avg_Expression"])
        else -float("inf"),
    )
    lowest_density_gene = min(
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
            "Lowest_Density_Gene": lowest_density_gene,
            "Lowest_Density_Value": row[f"{lowest_density_gene}_Density"],
        }
    )


def find_highest_expression_and_highest_density(row):
    """
    Find the gene with the highest expression and highest density in each row
    """
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
    # NOTE NOTE NOTE this is a toggle I am using
    # TODO TOGGLE THIS BACK for the bargraphs
    syntelog_data = syntelog_data.loc[syntelog_data["Set"] == "ABCD"]
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

    # ---------------------------------------------------
    # Now I need to make a bar plot of binned TE density vs expression for all
    # genes, and incorporate a variance to the expression. This figure is aimed
    # at emulating the green diamon figure from the Mimulus paper
    # Generate a bin plot of the density and expression data

    # Remove chloroplast genes
    density_and_expression = density_and_expression.loc[
        ~density_and_expression["Gene_Name"].str.contains("FxaCT")
    ]

    density_and_expression.set_index("Gene_Name", inplace=True)

    bins = np.arange(0, 1.05, 0.05)
    density_and_expression["Density_Bin"] = pd.cut(
        density_and_expression["Density"], bins=bins, include_lowest=True
    )

    # NOTE Jordan outlier rough estimates
    density_and_expression = density_and_expression.loc[
        (density_and_expression["Avg_Expression"] > 0.1)
        & (density_and_expression["Avg_Expression"] < 1000.0)
    ]

    print(density_and_expression)
    density_and_expression["Avg_Expression"] = np.log10(
        density_and_expression["Avg_Expression"] + 1
    )
    print(density_and_expression)

    grouped = (
        density_and_expression.groupby("Density_Bin")
        .agg(
            Gene_Count=("Density", "size"),
            Avg_Expression=("Avg_Expression", "mean"),
            Exp_Std_Dev=("Avg_Expression", "std"),
        )
        .reset_index()
    )

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.bar(
        grouped["Density_Bin"].astype(str),
        grouped["Gene_Count"],
        color="lightblue",
        label="Gene Count",
    )
    ax1.tick_params(axis="x", labelrotation=-45)
    ax2 = ax1.twinx()
    # ax2.plot(grouped["Density_Bin"].astype(str), grouped["Avg_Expression"], color="red")
    grouped.plot(
        x="Density_Bin",
        y="Avg_Expression",
        ax=ax2,
        color="black",
        yerr="Exp_Std_Dev",
        capsize=4,
        rot=0,
        logy=False,
    )
    # ax2.plot(grouped["Density_Bin"].astype(str), grouped["Avg_Expression"], color="red")
    # ax2.set_yscale("log")
    plt.show()
    raise ValueError("STOP")
    # ---------------------------------------------------

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

    # NOTE toggle
    result = merged_D.apply(find_highest_expression_and_lowest_density, axis=1)
    final_result = pd.concat([merged_D, result], axis=1)

    # subset the data so that I can find a good example manually, easily
    desired_output = final_result.loc[
        (final_result["Highest_Expression_Gene"] == "Gene_B_Name")
        & (final_result["Lowest_Density_Gene"] == "Gene_B_Name")
        & (final_result["Lowest_Density_Value"] < 0.20)
        & (final_result["Gene_A_Name_Density"] > 0.40)
        & (final_result["Gene_C_Name_Density"] > 0.40)
        & (final_result["Gene_D_Name_Density"] > 0.40)
    ]

    # NOTE toggle
    # gene_of_interest = "Fxa4Bg103345"
    # plot_gene_of_interest_density_and_expression(desired_output, gene_of_interest)
    # for gene_of_interest in desired_output["Gene_B_Name"].to_list():
    # plot_gene_of_interest_density_and_expression(desired_output, gene_of_interest)

    # NOTE toggle
    # result = merged_D.apply(find_highest_expression_and_highest_density, axis=1)
    # final_result = pd.concat([merged_D, result], axis=1)
    # desired_output = final_result.loc[
    #     (final_result["Highest_Expression_Gene"] == "Gene_B_Name")
    #     & (final_result["Highest_Density_Gene"] == "Gene_B_Name")
    #     & (final_result["Highest_Density_Value"] > 0.50)
    # ]

    # -----------------------------------------------------
    # NOTE
    # Example A: A gene that has the highest expression, and the highest TE
    # density
    # Fxa4Bg103345 is probably the best example
    # It has much higher expression than all other homeologs, and
    # has a CACTA with high identity right in front of it. This is a good example.
    # Other potential candidates:

    # Fxa5Bg100284 is a good example as it has crazy expression and has a big
    # piece of LTR, LINE, and some other upstream. Other homeologs have good
    # expression, 2 of the homeologs have 0 TE, and one has 60\% TE but has lowest
    # expression.

    # Fxa5Bg103293 has crazy expression and has a SINE, Mutator and LTR in
    # front of it. Only one other homelog has 1/4th of the expression and 0.25 TE
    # Density in front, the other 2 basically have 0 expression and no TE.

    # Fxa6Bg101326 has high expression and compares well to its homeologs

    # NOTE
    # Example B: A gene that has the highest expression, and the lowest TE density
    # Fxa5Bg101314 is the best example
    # OR Fxa3Bg203434
    # OR Fxa7Bg201270
    # OR Fxa7Bg201800
    # -----------------------------------------------------
