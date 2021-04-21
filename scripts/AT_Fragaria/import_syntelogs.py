#!/usr/bin/env python3

__author__ = "Scott Teresi"

import logging
import pandas as pd


def import_syntelogs(syntelog_input_file, genome_name):
    """
    Import the syntelogs from the raw file and manage data filtration
    """

    col_names = [
        "OrgA_Chromosome",
        "OrgA_Gene_Region",
        "OrgA_Start",
        "OrgA_Stop",
        "OrgB_Chromosome",
        "OrgB_Gene_Region",
        "OrgB_Start",
        "OrgB_Stop",
        "E_Value",
        "Diagonal_Score",
        "Web_Link",
    ]

    col_to_use = [
        "OrgA_Chromosome",
        "OrgA_Gene_Region",
        "OrgB_Chromosome",
        "OrgB_Gene_Region",
        "E_Value",
        "Diagonal_Score",
    ]

    gene_data = pd.read_csv(
        syntelog_input_file,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
    )

    # Set the correct data types
    gene_data.OrgA_Chromosome = gene_data.OrgA_Chromosome.astype(str)
    gene_data.OrgB_Chromosome = gene_data.OrgB_Chromosome.astype(str)
    gene_data.OrgA_Gene_Region = gene_data.OrgA_Gene_Region.astype(str)
    gene_data.OrgB_Gene_Region = gene_data.OrgB_Gene_Region.astype(str)
    gene_data.E_Value = gene_data.E_Value.astype("float64")
    gene_data.Diagonal_Score = gene_data.Diagonal_Score.astype("int32")

    # Get the correct name for the genes
    # MAGIC
    gene_data["OrgA_Gene_Region"] = (
        gene_data["OrgA_Gene_Region"].str.split("\|\|").str[3]
    )
    # MAGIC
    gene_data["OrgB_Gene_Region"] = (
        gene_data["OrgB_Gene_Region"].str.split("\|\|").str[3]
    )

    if genome_name == "H4":
        # Fix Cam genes
        gene_data["OrgA_Gene_Region"] = (
            gene_data["OrgA_Gene_Region"].str.split("-mRNA-1").str[0]
        )
        # Fix H4 genes
        # MAGIC
        gene_data["OrgB_Gene_Region"] = (
            gene_data["OrgB_Gene_Region"].str.split("\.1").str[0]
        )
        # This step is important, it could differ if your data input is different.
        gene_data.rename(
            columns={"OrgA_Gene_Region": "Camarosa", "OrgB_Gene_Region": genome_name,},
            inplace=True,
        )
    elif genome_name == "Camarosa":
        gene_data["OrgB_Gene_Region"] = (
            gene_data["OrgB_Gene_Region"].str.split("-mRNA-1").str[0]
        )
        # This step is important, it could differ if your data input is different.
        gene_data.rename(
            columns={
                "OrgA_Gene_Region": "Arabidopsis",
                "OrgB_Gene_Region": genome_name,
            },
            inplace=True,
        )
    elif genome_name == "Del_Norte":
        # Fix Cam genes
        gene_data["OrgA_Gene_Region"] = (
            gene_data["OrgA_Gene_Region"].str.split("-mRNA-1").str[0]
        )
        # Fix Del Norte genes
        gene_data["OrgB_Gene_Region"] = (
            gene_data["OrgB_Gene_Region"].str.split("-mRNA-1").str[0]
        )
        # This step is important, it could differ if your data input is different.
        gene_data.rename(
            columns={"OrgA_Gene_Region": "Camarosa", "OrgB_Gene_Region": genome_name,},
            inplace=True,
        )

    elif genome_name == "Iinumae":
        pass  # no renaming needs to be done for Iinumae
    else:
        raise ValueError("Genome name not Camarosa, H4, or Iinumae")

    # Get the correct name for the chromosome
    # gene_data["OrgA_Chromosome"] = gene_data["OrgA_Chromosome"].str.split("_").str[1]
    # gene_data["OrgB_Chromosome"] = gene_data["OrgB_Chromosome"].str.split("_").str[1]

    # Trim E-values less than 0.05
    # MAGIC
    # NOTE could make input argument
    gene_data = gene_data.loc[gene_data["E_Value"] < 0.05]

    gene_data.drop(
        columns=["OrgA_Chromosome", "OrgB_Chromosome", "Diagonal_Score"], inplace=True,
    )

    # Add column with identifier so we can later see what source we derived the
    # gene pair from
    gene_data["Point_of_Origin"] = "Synteny"

    # Need to take first occurrence of the gene, the one with the smallest
    # E-Value

    # TODO Debug
    # Pat thinks I am discarding a huge chunk of Camarosa
    # The master file is just missing some H4 entries for some Camarosa genes,
    # which would also be reflected in this proto-table here.
    # Basically I just have a shit ton of Null H4 values in the dataframe as a
    # result of this step.
    # It is throwing away whatever is not the top hit in the group by
    # Throwing awya 25% of H4 and 75% of Cam.

    # For every gene in Cam get the best hit in the H4 no matter
    # if it repeats elsewhere

    # For every gene in Del Norte get the best hit in Cam (the cam gene won't
    # repeat).

    # gene_data = gene_data.loc[gene_data.groupby(genome_name)["E_Value"].idxmin()]

    return gene_data


class Syntelog_Data(object):
    """
    Wrappers for input data, multiple syntelog pairs.

    Used to provide a common interface and fast calculations with numpy.
    """

    def __init__(self, syntelog_dataframe, logger=None):
        """Initialize.

        Args:
            syntelog_dataframe (DataFrame): syntelog data frame.
        """
        self._logger = logger or logging.getLogger(__name__)
        self.dataframe = syntelog_dataframe

    def save_to_disk(self, filename):
        """
        Save the syntelogs to disk in a 2-column format.
        Arabidopsis in left-hand column, strawberry in right-hand column.

        Args:
            filename (str): path for the file
        """
        # self.to_write = self.dataframe.copy(deep=True)
        # self.to_write.drop(
        # columns=["E_Value", "OrgA_Chromosome", "OrgB_Chromosome", "Diagonal_Score"],
        # inplace=True,
        # )
        self.dataframe.to_csv(filename, sep="\t", header=True, index=False)

    def __repr__(self):
        """
        String representation for developer.
        """
        return "Syntelog_Data{}".format(self.dataframe)
