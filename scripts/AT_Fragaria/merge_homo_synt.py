#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import logging


def merge_homo_synt(wrapped_syntelogs, wrapped_homologs, genome_name):
    """
    Only add a strawberry gene from the BLAST set (wrapped homologs) if that
    gene is not currently in the synteny set.

    Args:
        wrapped_syntelogs (Syntelog_Data object): wrapped data

        wrapped_homologs (Homolog_Data object ): wrapped data

    Returns:
        merged_all (pandas.core.Data.Frame): Pandas dataframe that is the
        amalgamation of the syntelog and homolog data, with duplicate genes
        coming from the homolog data being dropped.
    """
    missing_links = _identify_missing_links(
        wrapped_syntelogs, wrapped_homologs, genome_name
    )
    merged_all = pd.concat(
        [wrapped_syntelogs.dataframe, missing_links], axis=0, join="outer"
    )
    # merged_all.rename(
    # columns={"Arabidopsis": "Arabidopsis_Gene", genome_name: genome_name,},
    # inplace=True,
    # )
    return merged_all


def _identify_missing_links(wrapped_syntelogs, wrapped_homologs, genome_name):
    """
    Create and use a boolean mask to identify the strawberry genes in the
    homolog data that are missing from the syntelog data.

    Args:
        wrapped_syntelogs (Syntelog_Data object): wrapped data

        wrapped_homologs (Homolog_Data object ): wrapped data

    Returns:
        missing_links (pandas.core.Data.Frame): Pandas dataframe that is the
        amalgamation of the syntelog and homolog data, with duplicate genes
        coming from the homolog data being dropped.
    """
    syntelog_strawberry_list = wrapped_syntelogs.dataframe[genome_name].to_list()
    boolean_mask = wrapped_homologs.dataframe[genome_name].isin(
        syntelog_strawberry_list
    )
    missing_links = wrapped_homologs.dataframe[~boolean_mask]
    return missing_links


class Merged_Data(object):
    """
    Wrappers for ouput data, multiple syntelog and ortholog pairs.
    """

    def __init__(self, merged_dataframe, logger=None):
        """Initialize.

        Args:
            merged_dataframe (DataFrame): data frame.
        """
        self._logger = logger or logging.getLogger(__name__)
        self.dataframe = merged_dataframe

    def save_to_disk(self, filename):
        """
        Save the syntelogs to disk in a 2-column format.
        Arabidopsis in left-hand column, strawberry in right-hand column.

        Args:
            filename (str): path for file
        """
        self.to_write = self.dataframe.copy(deep=True)
        self.to_write.to_csv(filename, sep="\t", header=True, index=False)

    def __repr__(self):
        """
        String representation for developer.
        """
        return "Merged_Data{}".format(self.dataframe)

    @property
    def total_count_summary(self):
        """
        Summary of source of genes

        Returns:
            synteny_counts, blast_counts (tuple)
        """
        # total_pairs = self.dataframe.shape[0]
        synteny_counts = self.dataframe["Point_of_Origin"].value_counts().Synteny
        blast_counts = self.dataframe["Point_of_Origin"].value_counts().BLAST
        return synteny_counts, blast_counts
