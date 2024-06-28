__author__ = "Scott Teresi"

"""
TODO
"""

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd

from src.go_analysis.upset_plot import read_go_enrichment_table


def read_domestication_sweep_table(filename):
    data = pd.read_csv(
        filename,
        sep=",",
        names=[
            "Sweep_ID",
            "Chromsome_ID",
            "Number_of_Sig_Windows_in_Sweep",
            "Start",
            "Stop",
            "Max_XPLCR_Val",
            "Window_Start_For_Max_XPLCR_Val",
            "Window_Stop_For_Max_XPLCR_Val",
            "Size_of_Sweep",
            "Selection_Coefficient",
            "Category",
        ],
        header=0,
        dtype={
            "Sweep_ID": str,
            "Chromsome_ID": str,
            "Number_of_Sig_Windows_in_Sweep": int,
            "Start": int,
            "Stop": int,
            "Max_XPLCR_Val": float,
            "Window_Start_For_Max_XPLCR_Val": int,
            "Window_Stop_For_Max_XPLCR_Val": int,
            "Size_of_Sweep": int,
            "Selection_Coefficient": float,
            "Category": str,
        },
    )
    return data


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "preprocessed_go_enrichment_table",
        type=str,
        help="TODO",
    )
    parser.add_argument("genome_name", type=str, help="Shorthand name of the genome")

    parser.add_argument("domestication_sweep_table", type=str, help="TODO")
    args = parser.parse_args()
    args.preprocessed_go_enrichment_table = os.path.abspath(
        args.preprocessed_go_enrichment_table
    )
    args.domestication_sweep_table = os.path.abspath(args.domestication_sweep_table)

    go_enrichment_table = read_go_enrichment_table(
        args.preprocessed_go_enrichment_table
    )

    sweep_table = read_domestication_sweep_table(args.domestication_sweep_table)
    print(go_enrichment_table)
    print(go_enrichment_table.columns)
    print(sweep_table)

    print(sweep_table["Sweep_ID"].max())
