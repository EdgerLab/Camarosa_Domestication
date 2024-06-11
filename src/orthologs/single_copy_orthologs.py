#!/usr/bin/env python3

"""
- SCO stands for Single Copy Orthologs
- Calculate the TE density average of SCO genes, and compare to non-SCO genes
- Perform a Mann-Whitney U test to determine a significant difference in TE density
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
from src.parse_ka_ks import read_filtered_ka_ks_table


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
    return strawberry_ortholog_table.loc[strawberry_ortholog_table["SCO_Status"] == "Y"]


def merge_merged_SCO_table_w_TE_table(merged_orthologs, TE_table, common_gene_col):
    """
    Merge the merged strawberry|SCO table with the TE table, this will allow us
    to evaluate TE density in the context of SCOs
    """
    data = merged_orthologs.merge(TE_table, on=common_gene_col, how="inner")
    return data


def calc_test_statistic(sco_array, not_sco_array):
    """
    Test that the means of the SCO and non-SCO genes are different
    """
    # One-sided test, we want to know if the SCO genes have a lower TE density
    statistic = mannwhitneyu(sco_array, not_sco_array, alternative="less")
    return statistic


def identify_h4_scos(strawberry_ortholog_table, SCO_table, logger):
    """
    Return the set of "TRUE" SCO genes. Merge the vanilla SCO table with my
    H4-AT gene table and remove any duplicates.
    """
    data = strawberry_ortholog_table.copy(deep=True)
    sco_filter = strawberry_ortholog_table["Arabidopsis_Gene"].isin(
        SCO_table["Arabidopsis_Gene"]
    )
    data.loc[sco_filter, ["SCO_Status"]] = "Y"
    data.loc[~sco_filter, ["SCO_Status"]] = "N"

    data.drop(columns=["RR_Gene"], inplace=True)
    data.drop_duplicates(subset=["Arabidopsis_Gene", "H4_Gene"], inplace=True)
    data = data.loc[data["SCO_Status"] == "Y"]
    logger.info(f"Identified {len(data)} SCO AT genes in the H4 genome")

    # Remove any AT genes in my SCO-Strawberry table if those AT genes repeat,
    # we want to be consistent with the SCO definition
    # TODO have Pat check this one more time
    data.drop_duplicates(subset=["Arabidopsis_Gene"], keep=False, inplace=True)
    logger.info(
        f"""
        Removing duplicate AT genes (rows) to comply with the SCO
        definition, {len(data)} remaining unique SCO AT genes
        """
    )
    return data


def calc_and_compare_sco_TE_density(
    scos,
    not_scos,
    te_table,
    te_col,
    genome_name,
    logger,
):

    scos = merge_merged_SCO_table_w_TE_table(scos, te_table, f"{genome_name}_Gene")
    not_scos = merge_merged_SCO_table_w_TE_table(
        not_scos, te_table, f"{genome_name}_Gene"
    )

    te_col = f"{genome_name}_{te_col}"
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


if __name__ == "__main__":

    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "h4_preprocessed_density_table",
        type=str,
    )
    parser.add_argument(
        "rr_preprocessed_density_table",
        type=str,
    )
    parser.add_argument("strawberry_ortholog_table", type=str)
    parser.add_argument("sco_table", type=str)
    parser.add_argument(
        "output_dir",
        type=str,
        help="parent directory to output results",
    )
    parser.add_argument("ka_ks_table", type=str)

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.h4_preprocessed_density_table = os.path.abspath(
        args.h4_preprocessed_density_table
    )
    args.rr_preprocessed_density_table = os.path.abspath(
        args.rr_preprocessed_density_table
    )
    args.strawberry_ortholog_table = os.path.abspath(args.strawberry_ortholog_table)
    args.sco_table = os.path.abspath(args.sco_table)
    args.output_dir = os.path.abspath(args.output_dir)
    args.ka_ks_table = os.path.abspath(args.ka_ks_table)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -----------------------------------------------------------
    print()
    # Load in the pre-filtered TE Density data
    h4_te_table = pd.read_csv(
        args.h4_preprocessed_density_table, header="infer", sep="\t"
    )
    rr_te_table = pd.read_csv(
        args.rr_preprocessed_density_table, header="infer", sep="\t"
    )

    # Load in the pre-filtered ka_ks data
    ka_ks_table = read_filtered_ka_ks_table(args.ka_ks_table, "RR")

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
    # NOTE THIS MUST BE THE SAME FOR H4 AND RR
    filename = os.path.splitext(os.path.basename(args.h4_preprocessed_density_table))[0]
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
    strawberry_ortholog_table.sort_values(by="Arabidopsis_Gene", inplace=True)

    # Drop columns we don't need
    strawberry_ortholog_table.drop(
        columns=[
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

    genome_name = "H4"  # MAGIC
    H4_AT_scos = identify_h4_scos(strawberry_ortholog_table, SCO_table, logger)
    true_AT_scos = H4_AT_scos["Arabidopsis_Gene"].tolist()

    # Identify genes that are not in the SCO table
    H4_not_sco = strawberry_ortholog_table.copy(deep=True)
    H4_not_sco = H4_not_sco.drop(columns=["RR_Gene"])
    H4_not_sco.drop_duplicates(subset=["Arabidopsis_Gene", "H4_Gene"], inplace=True)
    H4_not_sco = H4_not_sco.loc[~H4_not_sco["Arabidopsis_Gene"].isin(true_AT_scos)]
    H4_not_sco.loc[:, "SCO_Status"] = "N"

    logger.info("H4 SCO Analysis")
    calc_and_compare_sco_TE_density(
        H4_AT_scos,
        H4_not_sco,
        h4_te_table,
        te_col,
        "H4",
        logger,
    )
    print()

    # Now do the RR genes that Pat requested
    RR_not_sco = strawberry_ortholog_table.copy(deep=True)
    RR_not_sco = RR_not_sco.drop(columns=["H4_Gene"])

    # TODO check this usage of drop duplicates, do we want to drop all duplicates?
    RR_not_sco.drop_duplicates(subset=["Arabidopsis_Gene", "RR_Gene"], inplace=True)
    RR_not_sco = RR_not_sco.loc[~RR_not_sco["Arabidopsis_Gene"].isin(true_AT_scos)]
    RR_not_sco.loc[:, "SCO_Status"] = "N"
    print(RR_not_sco)
    raise ValueError("STOP")

    # Find the RR SCO genes, don't remove duplicates
    RR_sco = strawberry_ortholog_table.copy(deep=True)
    RR_sco = RR_sco.drop(columns=["H4_Gene"])
    RR_sco.drop_duplicates(subset=["Arabidopsis_Gene", "RR_Gene"], inplace=True)
    RR_sco = RR_sco.loc[RR_sco["Arabidopsis_Gene"].isin(true_AT_scos)]

    # Count the number of times an "Arabidopsis_Gene" only appears once
    gene_counts = RR_sco["Arabidopsis_Gene"].value_counts()
    only_once = gene_counts[gene_counts == 1].index
    only_four_times = gene_counts[gene_counts == 4].index
    rr_4 = RR_sco.loc[RR_sco["Arabidopsis_Gene"].isin(only_four_times)]
    rr_1 = RR_sco.loc[RR_sco["Arabidopsis_Gene"].isin(only_once)]

    logger.info("RR SCO Analysis")
    calc_and_compare_sco_TE_density(
        RR_sco,
        RR_not_sco,
        rr_te_table,
        te_col,
        "RR",
        logger,
    )

    merged_RR_sco = RR_sco.merge(ka_ks_table, on="RR_Gene", how="inner")

    # Merging ~10K gene list from KA_KS analysis with 64K gene list from
    # ortholog table
    merged_RR_nonsco = RR_not_sco.merge(ka_ks_table, on="RR_Gene", how="inner")

    # NOTE well this is good, KA_KS is basically 1 for the SCO genes
    print(merged_RR_sco)
    print(merged_RR_sco["KA_KS"].mean())

    # NOTE KA_KS for the non-SCO genes
    print(merged_RR_nonsco)
    print(merged_RR_nonsco["KA_KS"].mean())
