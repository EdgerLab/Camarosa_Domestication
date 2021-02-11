#!/usr/bin/env python3

"""
Import the homolog data and provide a class for its access
"""

__author__ = "Scott Teresi"

import logging
import pandas as pd


def import_homologs(homolog_input_file, genome_name):
    """
    Import the homologs from the raw file and manage data filtration
    """
    col_names = [
        "Query",
        "Subject",
        "Percent ID",
        "Alignment Length",
        "Mismatches",
        "Gap Openings",
        "Q_Start",
        "Q_Stop",
        "Subject_Start",
        "Subject_Stop",
        "E_Value",
        "Bit_Score",
    ]

    col_to_use = [
        "Query",
        "Subject",
        "E_Value",
    ]

    homolog_dataframe = pd.read_csv(
        homolog_input_file,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
    )
    homolog_dataframe.E_Value = homolog_dataframe.E_Value.astype("float64")

    # Get the correct name for the arabidopsis genes
    homolog_dataframe["Subject"] = homolog_dataframe["Subject"].str.split("\|\|").str[3]

    # Get the correct name for the strawberry genes
    homolog_dataframe["Query"] = homolog_dataframe["Query"].str.split("-mRNA-1").str[0]

    # Trim E-values less than 0.05
    homolog_dataframe = homolog_dataframe.loc[homolog_dataframe["E_Value"] < 0.05]

    # Need to take first occurrence of the gene, the one with the smallest
    # E-Value
    homolog_dataframe = homolog_dataframe.loc[
        homolog_dataframe.groupby("Query")["E_Value"].idxmin()
    ]

    # Rename columns
    homolog_dataframe.rename(
        columns={"Query": genome_name, "Subject": "Arabidopsis"}, inplace=True
    )

    # Add column with identifier so we can later see what source we derived the
    # gene pair from
    homolog_dataframe["Point_of_Origin"] = "BLAST"

    return homolog_dataframe


class Homolog_Data(object):
    """
    Wrappers for input data, multiple homolog pairs.

    Used to provide a common interface and fast calculations with numpy.

    """

    def __init__(self, homolog_dataframe, logger=None):
        """Initialize.

        Args:
            homolog_dataframe (DataFrame): homolog data frame.
        """
        self._logger = logger or logging.getLogger(__name__)
        self.dataframe = homolog_dataframe

    def save_to_disk(self, filename):
        """
        Save the syntelogs to disk in a 2-column format.
        Arabidopsis in left-hand column, strawberry in right-hand column.

        Args:
            filename (str): path for the file
        """
        self.dataframe.to_csv(filename, sep="\t", header=True, index=False)

    def __repr__(self):
        """
        String representation for developer.
        """
        return "Homolog_Data{}".format(self.dataframe)
