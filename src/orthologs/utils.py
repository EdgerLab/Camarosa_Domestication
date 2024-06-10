def dn_chromosome_dict():
    """
    Dictionary to replace the chromosome names in the DN dataset to match the
    RR chromosome names, this dictionary was taken from the Hardigan paper
    """
    mapper = {
        "Fvb1-1": "1D",
        "Fvb1-2": "1B",
        "Fvb1-3": "1C",
        "Fvb1-4": "1A",
        "Fvb2-1": "2C",
        "Fvb2-2": "2A",
        "Fvb2-3": "2D",
        "Fvb2-4": "2B",
        "Fvb3-1": "3D",
        "Fvb3-2": "3B",
        "Fvb3-3": "3C",
        "Fvb3-4": "3A",
        "Fvb4-1": "4D",
        "Fvb4-2": "4C",
        "Fvb4-3": "4A",
        "Fvb4-4": "4B",
        "Fvb5-1": "5A",
        "Fvb5-2": "5D",
        "Fvb5-3": "5B",
        "Fvb5-4": "5C",
        "Fvb6-1": "6A",
        "Fvb6-2": "6C",
        "Fvb6-3": "6B",
        "Fvb6-4": "6D",
        "Fvb7-1": "7C",
        "Fvb7-2": "7A",
        "Fvb7-3": "7B",
        "Fvb7-4": "7D",
    }
    return mapper


def map_names(dataframe, chromosome_col, mapper=dn_chromosome_dict()):
    """
    Using a dictionary replace the values in a dataframe column with the key:values in
    the dictionary

    Args:
        dataframe (pd.DataFrame)
        chromosome_col (str): The name of the column to map, doesn't
            necessarily have to be a chromosome column
        mapper (dict): The dictionary to use to map the values

    Returns:
        dataframe (pd.DataFrame): With updated values
    """
    dataframe[chromosome_col] = dataframe[chromosome_col].replace(mapper)
    return dataframe


def reformat_chromosomes_from_SynMap(dataframe, chromosome_col):
    """
    Chromosomes that are of the form a6734_1A, we want to remove the a6734_ and
    take the second (MAGIC first index) part of the string.
    """
    dataframe[chromosome_col] = dataframe[chromosome_col].str.split("_", n=1).str[1]
    return dataframe


def reformat_gene_names_from_SynMap(dataframe, gene_col):
    """
    Take the SynMap output table in Pandas DataFrame format and reformat the
    column that pertains to the gene names. It remove the extraneous '|'
    characters by performing a string split with a MAGIC number
    """
    dataframe[gene_col] = dataframe[gene_col].str.split("\|\|").str[3]
    return dataframe


def reformat_gene_names_with_period(dataframe, gene_col):
    """
    Ex: Fxa1Ag100003.1
    Return Fxa1Ag100003
    """
    dataframe[gene_col] = dataframe[gene_col].str.split(".").str[0]
    return dataframe


def drop_rows_with_bad_val_in_col(dataframe, bad_val, col):
    """
    Helper function to drop rows in a dataframe where a column contains the
    'bad_val'.
    """

    dataframe = dataframe.loc[~dataframe[col].str.contains(bad_val)]
    return dataframe


def remove_str_from_val(dataframe, str_val, col):
    """
    DO NOT USE THIS TO FILTER NUMERIC "STRINGS" FROM A COLUMN

    Helper function to remove a string from each value in a column in a dataframe
    This is done by replacing the string with an empty string "".

    Args:
        dataframe (pd.DataFrame)
        str_val (str): The value to remove
        col (str): The name of the column
    """
    dataframe[col] = dataframe[col].str.replace(str_val, "")
    return dataframe
