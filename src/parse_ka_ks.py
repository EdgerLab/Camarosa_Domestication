#!/usr/bin/env python3

"""
- Reads Pat's KA_KS tables from disk and filters them, saving a new table to disk.
"""

__author__ = "Scott Teresi"


import pandas as pd
import numpy as np
import argparse
import os
import logging
import coloredlogs

import matplotlib.pyplot as plt

from src.orthologs.utils import remove_str_from_val, map_names
from src.orthologs.replace_and_reformat_DN_RR_BLAST_results import (
    import_decoder_ring,
    blacklist_if_no_new_name,
)


def read_raw_ka_ks_table(filepath, other_species=None):
    if other_species is None:
        raise ValueError("Set a species name, trying to reduce duplicate code")
    table = pd.read_csv(
        filepath,
        delimiter="\t",
        header=0,
        names=[
            "KS",
            "KA",
            "H4_Chromosome",
            "H4_Gene",
            f"{other_species}_Chromosome",
            f"{other_species}_Gene",
        ],
        usecols=["KS", "KA", "H4_Gene", f"{other_species}_Gene"],
        dtype={
            "KS": np.float64,
            "KA": np.float64,
            "H4_Gene": str,
            f"{other_species}_Gene": str,
        },
    )
    table.dropna(subset=["KS", "KA"], inplace=True)
    table["KA_KS"] = table["KA"] / table["KS"]
    table.replace([np.inf, -np.inf], np.nan, inplace=True)
    table.dropna(subset=["KA_KS"], inplace=True)

    # Reorder columns for easier reading
    table = table.reindex(
        columns=[
            "KA_KS",
            "KA",
            "KS",
            "H4_Gene",
            f"{other_species}_Gene",
        ],
        copy=True,
    )

    table = filter_table(table)
    return table


def filter_h4_rr_more(h4_rr):
    # Remove non primary transcripts
    h4_rr = h4_rr.loc[h4_rr["RR_Gene"].str.endswith(".1")]
    h4_rr = h4_rr.copy(deep=True)
    h4_rr.loc[:, "RR_Gene"] = h4_rr.loc[:, "RR_Gene"].str[:-2]
    return h4_rr


def filter_h4_dn_more(h4_dn):
    # Remove non primary transcripts
    h4_dn = h4_dn.loc[h4_dn["DN_Gene"].str.endswith("-mRNA-1")]
    h4_dn = h4_dn.copy(deep=True)
    h4_dn.loc[:, "DN_Gene"] = h4_dn.loc[:, "DN_Gene"].str.removesuffix("-mRNA-1")

    # Replace the old names because Pat used an annotation that contains some
    # of the old names
    h4_dn = map_names(h4_dn, "DN_Gene", decoder_ring)
    h4_dn = blacklist_if_no_new_name(h4_dn, decoder_ring, column="DN_Gene")

    return h4_dn


def filter_table(table):
    # Drop rows with KS greater than 5
    table = table.loc[table["KS"] <= 5]

    # Remove non primary H4 transcripts
    table = table.loc[table["H4_Gene"].str.endswith(".1")]
    table.loc[:, "H4_Gene"] = table.loc[:, "H4_Gene"].str[:-2]

    # Remove anything that is NaN (0/0) or other causes
    table = table.loc[~table["KA_KS"].isna()]
    return table


def read_filtered_ka_ks_table(filepath, other_species=None):
    if other_species is None:
        raise ValueError("Set a species name, trying to reduce duplicate code")
    table = pd.read_csv(
        filepath,
        delimiter="\t",
        header=0,
        names=[
            "KA_KS",
            "KS",
            "KA",
            "H4_Gene",
            f"{other_species}_Gene",
        ],
        dtype={
            "KA_KS": np.float64,
            "KS": np.float64,
            "KA": np.float64,
            "H4_Gene": str,
            f"{other_species}_Gene": str,
        },
    )
    return table


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument("h4_rr_raw_ka_ks", type=str, help="TODO")
    parser.add_argument("h4_dn_raw_ka_ks", type=str, help="TODO")
    parser.add_argument(
        "decoder_ring_input_file",
        type=str,
        help="parent path of decoder ring file",
    )

    parser.add_argument("h4_rr_filtered_ka_ks", type=str, help="TODO")
    parser.add_argument("h4_dn_filtered_ka_ks", type=str, help="TODO")

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.h4_rr_raw_ka_ks = os.path.abspath(args.h4_rr_raw_ka_ks)
    args.h4_dn_raw_ka_ks = os.path.abspath(args.h4_dn_raw_ka_ks)
    args.h4_rr_filtered_ka_ks = os.path.abspath(args.h4_rr_filtered_ka_ks)
    args.h4_dn_filtered_ka_ks = os.path.abspath(args.h4_dn_filtered_ka_ks)
    args.decoder_ring_input_file = os.path.abspath(args.decoder_ring_input_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # ------------------------------------------------------------------
    decoder_ring = import_decoder_ring(args.decoder_ring_input_file)

    h4_rr = read_raw_ka_ks_table(args.h4_rr_raw_ka_ks, other_species="RR")
    h4_rr = filter_h4_rr_more(h4_rr)

    h4_dn = read_raw_ka_ks_table(args.h4_dn_raw_ka_ks, other_species="DN")
    h4_dn = filter_h4_dn_more(h4_dn)

    print(h4_dn["KA_KS"].describe())
    print(h4_rr["KA_KS"].describe())

    # print(h4_dn.loc[h4_dn["DN_Gene"].str.contains("RagTag")])
    # print()
    # print(h4_dn.loc[h4_dn["DN_Gene"].str.contains("pilon")])
    # print()
    # print(h4_dn.loc[h4_dn["DN_Gene"].str.contains("genemark")])

    # TODO resolve this issue later
    h4_dn = h4_dn.loc[~h4_dn["DN_Gene"].str.contains("RagTag")]
    h4_dn = h4_dn.loc[~h4_dn["DN_Gene"].str.contains("pilon")]

    h4_dn.to_csv(args.h4_dn_filtered_ka_ks, sep="\t", index=False)
    h4_rr.to_csv(args.h4_rr_filtered_ka_ks, sep="\t", index=False)
