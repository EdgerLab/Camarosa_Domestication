import numpy as np
import pandas as pd


def calc_fpkm(count_matrix, lengths):
    """Calculates fragments per kilobase transcript per million reads.

    Args:
        count_matrix (pandas.DataFrame): A pandas dataframe of counts. Shape is
        (N_Genes, N_Samples).

        lengths (pandas.Data.Frame): A pandas dataframe of gene lengths. Shape
        is (N_Genes,).

    Returns:
        normalized (numpy array): A numpy array of FPKM values. Shape is
        (N_Genes, N_Samples).
    """

    lengths = lengths.to_numpy()
    counts = count_matrix.to_numpy()
    N = counts.sum(axis=0)
    L = lengths
    C = counts
    normalized = 1e9 * C / (N[np.newaxis, :] * L[:, np.newaxis])

    return normalized
