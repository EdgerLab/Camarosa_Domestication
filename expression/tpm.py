import numpy as np
import pandas as pd


def calc_tpm(count_matrix, lengths):
    """Calculate TPM matrix given a matrix of counts and a matrix of summed
    exon lengths.

    Args:
        count_matrix (pandas.DataFrame): A pandas dataframe of counts. Shape is
        (N_Genes, N_Samples).

        lengths (pandas.Data.Frame): A pandas dataframe of gene lengths. Shape
        is (N_Genes,).

    Returns:
        normalized (numpy array): A numpy array of TPM values. Shape is
        (N_Genes, N_Samples).
    """

    lengths = lengths.to_numpy()
    counts = count_matrix.to_numpy()
    N = counts.sum(axis=0)
    L = lengths
    C = counts

    A = (C * 1e3) / (L[:, np.newaxis])
    TPM = A * 1e6 * (1 / A.sum(axis=0))

    return TPM
