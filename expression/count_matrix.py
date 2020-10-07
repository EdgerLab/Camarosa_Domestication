import pandas as pd


def import_count_matrix(matrix_input_file):
    """
    matrix_input_file (str): String path to count matrix file.
    """
    counts = pd.read_csv(matrix_input_file, sep="\t", header="infer")
    # Remove all the extraneous text from the HTSeq and R transformations
    counts.columns = counts.columns.str.rstrip("_001_trimmed_P1_Results.count")
    counts.rename(columns={"G": "Gene_Name"}, inplace=True)  # rename column
    counts.drop(counts.tail(5).index, inplace=True)  # drop last n rows,
    # inherent to our count file, extraneous output

    return counts
