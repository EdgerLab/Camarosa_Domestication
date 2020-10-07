import pandas as pd


def import_genes(genes_input_path):
    """Import gene annotation file as a Pandas DataFrame.

    Args:
        genes_input_path (str, argparse): Specify the input path of the gene
        annotation data. GFF format.

    Returns:
        exon_sums (Pandas DataFrame): Pandas dataframe of shape (N_Genes,).
        Second column is the length of the gene, which is calculated by summing
        the exon lengths.

    """

    col_names = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Score",
        "Strand",
        "Frame",
        "FullName",
    ]

    col_to_use = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Strand",
        "FullName",
    ]

    Gene_Data = pd.read_csv(
        genes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
    )

    exons = Gene_Data[Gene_Data.Feature == "exon"].copy(deep=True)

    # These commands will need to be changed if you are using a different
    # annotation as the naming scheme may be a little different.
    exons[["Name1", "Gene_Portion"]] = exons.FullName.str.split(";Parent=", expand=True)
    exons[["Gene_Name", "mRNA"]] = exons.Gene_Portion.str.split("-mRNA-", expand=True)
    exons.drop(
        ["FullName", "Name1", "Gene_Portion", "mRNA", "Feature", "Software", "Strand"],
        axis=1,
        inplace=True,
    )
    # Calculate the correct 'gene' lengths for normalization in the fpkm
    # calculation
    exons["Total_Exon_Length"] = exons["Stop"] - exons["Start"] + 1
    # group by gene name to get the correct summation
    exon_sums = exons.groupby(["Gene_Name"]).Total_Exon_Length.sum()

    return exon_sums
