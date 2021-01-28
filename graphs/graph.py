#!/usr/bin/env/python

"""
Execute graphing commands
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd
import h5py

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

from scripts.density_data import DensityData
from graphs.dotplots import *
from transposon.gene_data import GeneData


def test_data():
    """remove at later point"""
    # TEST DATA
    table = np.array(range(30))
    table = table.reshape(3, 2, 5)
    # Functionally equivalent
    # table[1][1][:]
    # table[1, 1, :]


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = "/home/scott/Documents/Uni/Research/Projects/Domestication_Data/"
    parser = argparse.ArgumentParser(description="generate graphs")
    test_file = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/ShortTest_Fvb1-1.h5"

    genes = GeneData(
        pd.read_csv(
            "/home/scott/Documents/Uni/Research/Projects/TE_Data/filtered_input_data/Cleaned_ShortFvb1-1_Genes.tsv",
            sep="\t",
            header="infer",
            dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
            index_col="Gene_Name",  # this is crucial
        ),
        "Camarosa",
    )
    strands_as_numpy = genes.strands.to_numpy(copy=False)
    # one_gene_list = np.nonzero(strands_as_numpy)[0]  # gets indices of 1s
    zero_gene_list = np.where(strands_as_numpy == 0)[0]  # gets indices of 0s
    # NB the 0s are the antisense
    my_names = genes.data_frame.iloc[zero_gene_list, :].index.tolist()

    density_data = DensityData(test_file, my_names)

    plot_intra_density(density_data, "Superfamily", output_default)
    plot_intra_density(density_data, "Order", output_default)
    plot_density_all(density_data, "Superfamily", output_default)
    plot_density_all(density_data, "Order", output_default)
