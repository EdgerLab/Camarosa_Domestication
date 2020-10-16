#!/usr/bin/env python3

"""
Import expression matrix
"""

__author__ = "Scott Teresi"

import pandas as pd


def import_exp_dataframe(filename):
    """

    """

    data = pd.read_csv(filename, sep="\t+", header="infer", engine="python",)
    data.set_index("Gene_Name", inplace=True)

    return data
