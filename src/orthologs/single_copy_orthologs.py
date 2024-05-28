#!/usr/bin/env python3

"""
- SCO stands for Single Copy Orthologs
- Calculate the TE density average of SCO genes...
- TODO
"""

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np
import os
import argparse
import logging
import coloredlogs
from scipy.stats import mannwhitneyu

from src.syntelog_differences.bargraphs import decode_te_window_direction_str
from src.extract_AED_score import read_aed_output_table
from src.go_analysis.find_abnormal_genes import subset_by_aed_score


def import_SCO_table(filepath):
    data = pd.read_csv(
        filepath,
        header=0,
        names=["Arabidopsis_Gene", "Function"],
        usecols=["Arabidopsis_Gene"],
        dtype={"Arabidopsis_Gene": str},
        sep="\t",
    )
    return data


def tag_sco_genes(strawberry_ortholog_table, SCO_table):
    """
    Add a column to the strawberry ortholog table that identifies if an AT gene
    is a single copy ortholog
    """
    # Value will be 'Y' if a gene in the column 'Arabidopsis_Gene' is in the
    # SCO table under the column 'Arabidopsis_Gene'
    sco_filter = strawberry_ortholog_table["Arabidopsis_Gene"].isin(
        SCO_table["Arabidopsis_Gene"]
    )
    strawberry_ortholog_table.loc[sco_filter, ["SCO_Status"]] = "Y"
    strawberry_ortholog_table.loc[~sco_filter, ["SCO_Status"]] = "N"
    return strawberry_ortholog_table


def subset_strawberry_orthologs_by_SCO_status(strawberry_ortholog_table):
    """
    TODO
    """
    return strawberry_ortholog_table.loc[strawberry_ortholog_table["SCO_Status"] == "Y"]


def merge_merged_SCO_table_w_TE_table(merged_orthologs, TE_table, common_gene_col):
    """
    Merge the merged strawberry|SCO table with the TE table, this will allow us
    to evaluate TE density in the context of SCOs
    """
    data = merged_orthologs.merge(TE_table, on=common_gene_col, how="inner")
    return data


def filter_merged_SCO_table():
    """
    Sanity check the merger of the merged strawberry|SCO table, there might be
    some genes were we have duplicates

    TODO checkpoint this step with Pat.
    """
    raise NotImplementedError("Function not implemented yet")


def calc_genome_wide_average_density():
    raise NotImplementedError("Function not implemented yet")


def identify_genes_not_in_SCO():
    """
    Return a table of form merged strawberry|SCO but only genes that are NOT in
    the single copy ortholog state

    TODO edit this docstring, I expect it to change
    """
    raise NotImplementedError("Function not implemented yet")


def identify_genes_in_SCO():
    """
    Return a table of form merged strawberry|SCO but only genes that are only
    in the single copy ortholog state
    """
    raise NotImplementedError("Function not implemented yet")


def calc_test_statistic(sco_array, not_sco_array):
    """
    Test that the means of the SCO and non-SCO genes are different

    TODO consider using the Welch's approximate t-test if the variances are
    unequal

    TODO we may also want to consider comparing the variances with an F-test
    but that is really sensitive to departures from normality

    TODO consider using a non-parametric test like the Mann-Whitney U test
    """
    # One-sided test, we want to know if the SCO genes have a lower TE density
    statistic = mannwhitneyu(sco_array, not_sco_array, alternative="less")
    return statistic


if __name__ == "__main__":

    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "preprocessed_density_table",
        type=str,
        help="TODO",
    )
    parser.add_argument("strawberry_ortholog_table", type=str, help="TODO")
    parser.add_argument("sco_table", type=str, help="TODO")
    parser.add_argument("aed_score_table", type=str, help="TODO")
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
    args.sco_table = os.path.abspath(args.sco_table)
    args.aed_score_table = os.path.abspath(args.aed_score_table)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -----------------------------------------------------------
    print()  # TODO remove
    # Load in the pre-filtered TE Density data
    te_table = pd.read_csv(args.preprocessed_density_table, header="infer", sep="\t")

    # Load in the AED scores
    aed_scores = read_aed_output_table(args.aed_score_table)

    # Load in the pre-filtered strawberry_ortholog information
    strawberry_ortholog_table = pd.read_csv(
        args.strawberry_ortholog_table, header="infer", sep="\t", low_memory=False
    )
    # Reorder columns for easier reading
    strawberry_ortholog_table.insert(
        0, "Arabidopsis_Gene", strawberry_ortholog_table.pop("Arabidopsis_Gene")
    )
    strawberry_ortholog_table.insert(
        1, "H4_Gene", strawberry_ortholog_table.pop("H4_Gene")
    )
    strawberry_ortholog_table.insert(
        2, "RR_Gene", strawberry_ortholog_table.pop("RR_Gene")
    )

    # Load in the SCO table
    # These are all unique AT genes
    SCO_table = import_SCO_table(args.sco_table)
    logger.info(f"Loaded {len(SCO_table)} unique SCO (AT) genes")

    # NOTE, duplicate code with find_abnormal_genes.py
    # TODO clean up
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

    te_col = f"{te_type}_{window}_{direction}"

    # -----------------------------------------------------------
    strawberry_ortholog_table = tag_sco_genes(strawberry_ortholog_table, SCO_table)
    strawberry_ortholog_table.sort_values(by="Arabidopsis_Gene", inplace=True)

    # Drop columns we don't need
    strawberry_ortholog_table.drop(
        columns=[
            "RR_Gene",
            "DN_Gene",
            "DN_Chromosome",
            "RR_Chromosome",
            "H4_Chromosome",
            "DN_RR_BLAST_E_Value",
            "DN_RR_Syntelog_E_Value",
            "DN_RR_Point_of_Origin",
            "H4_RR_BLAST_E_Value",
            "H4_RR_Syntelog_E_Value",
            "H4_RR_Point_of_Origin",
            "GO_ID",
            "GO_Term_Description",
        ],
        inplace=True,
    )

    # NOTE TODO REVIEW THE FOLLOWING CODE is heavily dependent on the order of
    # things being executed, clarify and do a code review with Pat.

    # H4 analysis
    # First drop duplicate gene pairs that are from when we merged with the
    # other strawberry genomes
    strawberry_ortholog_table.drop_duplicates(
        subset=["Arabidopsis_Gene", "H4_Gene"], inplace=True
    )

    # Merge the TE and merged ortholog table
    h4_merged = merge_merged_SCO_table_w_TE_table(
        strawberry_ortholog_table, te_table, "H4_Gene"
    )

    # TODO check this step with Pat
    scos = h4_merged.loc[h4_merged["SCO_Status"] == "Y"]
    scos = scos.copy(deep=True)
    logger.info(f"Identified {len(scos)} SCO AT genes in the H4 genome")
    # Remove any AT genes in my SCO-Strawberry table if those AT genes repeat
    scos.drop_duplicates(subset=["Arabidopsis_Gene"], keep=False, inplace=True)
    logger.info(
        f"""
        Removing duplicate AT genes to comply with the SCO
        definition, {len(scos)} remaining unique SCO AT genes
        """
    )

    # Create a table of AT-Strawberry genes where the AT gene is not a SCO
    # AT genes may repeat in this table
    not_scos = h4_merged.loc[h4_merged["SCO_Status"] == "N"]
    good_aed_score = 0.75  # MAGIC TODO make this a parameter

    # NOTE hard-coded value, will break when we get to RR
    aed_scores.rename(columns={"Gene_Name": "H4_Gene"}, inplace=True)
    not_scos = subset_by_aed_score(not_scos, aed_scores, "H4", good_aed_score)

    # TODO redefine this with a genome name, parametrize
    # NOTE the effect size ....
    te_col = f"H4_{te_col}"
    avg_sco = scos[te_col].mean()
    std_sco = scos[te_col].std()
    avg_not_sco = not_scos[te_col].mean()
    std_not_sco = not_scos[te_col].std()

    # Print out the effect size and standard deviation
    logger.info(f"SCO average: {avg_sco}, transformed: {avg_sco*5000:0.2f} BP")
    logger.info(
        f"Non-SCO average: {avg_not_sco}, transformed: {avg_not_sco*5000:0.2f} BP"
    )
    logger.info(f"BP difference: {avg_sco*5000 - avg_not_sco*5000:0.2f}")
    logger.info(f"SCO Standard Deviation: {std_sco}")
    logger.info(f"Non-SCO Standard Deviation: {std_not_sco}")

    statistic = calc_test_statistic(
        scos[te_col].to_numpy(), not_scos[te_col].to_numpy()
    )
    logger.info(f"My test statistic is {statistic}")
