#!/usr/bin/env/python

"""
Read TE Density H5 files
"""

__author__ = "Scott Teresi"

import h5py
import numpy as np


class DensityData:
    def __init__(self, input_h5, genes_to_swap):
        """
        input_h5 (str): Path to h5 file of TE Density output.
        """

        self.data_frame = h5py.File(input_h5, "r+")
        self.key_list = list(self.data_frame.keys())
        self.gene_list = self.data_frame["GENE_NAMES"][:]
        self.num_genes = len(self.gene_list)
        self.chromosomes = self.data_frame["CHROMOSOME_ID"]
        self.windows = self.data_frame["WINDOWS"]  # not int
        self.window_list = [int(i) for i in self.windows[:]]  # list of ints

        self.order_list = sorted(list(self.data_frame["ORDER_NAMES"][:]))
        self.super_list = sorted(list(self.data_frame["SUPERFAMILY_NAMES"][:]))

        # NB. Shape for these is (type of TE, window, gene)
        self.left_orders = self.data_frame["RHO_ORDERS_LEFT"]
        self.intra_orders = self.data_frame["RHO_ORDERS_INTRA"]
        self.right_orders = self.data_frame["RHO_ORDERS_RIGHT"]

        self.left_supers = self.data_frame["RHO_SUPERFAMILIES_LEFT"]
        self.intra_supers = self.data_frame["RHO_SUPERFAMILIES_INTRA"]
        self.right_supers = self.data_frame["RHO_SUPERFAMILIES_RIGHT"]

        self._swap_strand_vals(genes_to_swap)

    def index_of_gene(self, gene_string):
        """Return the index of a gene given the name of the gene
        Args:
            gene_string (str): string representing the name of the gene
        Returns:
            Returns an index of the gene in the H5 dataset
        """
        return np.where(self.gene_list == gene_string)[0][0]  # MAGIC

    def nonzero_indices(self):
        pass
        # print(np.nonzero(data.left_orders[:]))

    @property
    def order_index_dict(self):
        """Returns a dictionary of TE order names as keys and indices as values"""
        order_dict = {}
        for i in range(len(self.order_list)):
            order_dict[self.order_list[i]] = i
        return order_dict

    @property
    def super_index_dict(self):
        """Returns a dictionary of TE superfamily names as keys and indices as
        values"""
        super_dict = {}
        for i in range(len(self.super_list)):
            super_dict[self.super_list[i]] = i
        return super_dict

    def _swap_strand_vals(self, gene_names):
        """Switch density values for the genes in which it is antisense due to
        the fact that antisense genes point in the opposite direction to sense
        genes

        Args:
            gene_names(list of str):
        """
        for name in gene_names:
            index_to_switch = self.index_of_gene(name)

            # SUPERFAMILIES
            left_val_super = self.data_frame["RHO_SUPERFAMILIES_LEFT"][
                :, :, index_to_switch
            ]
            right_val_super = self.data_frame["RHO_SUPERFAMILIES_RIGHT"][
                :, :, index_to_switch
            ]
            # ORDERS
            left_val_order = self.data_frame["RHO_ORDERS_LEFT"][:, :, index_to_switch]

            right_val_order = self.data_frame["RHO_ORDERS_RIGHT"][:, :, index_to_switch]

            # Reassign
            self.data_frame["RHO_SUPERFAMILIES_RIGHT"][
                :, :, index_to_switch
            ] = left_val_super

            self.data_frame["RHO_SUPERFAMILIES_LEFT"][
                :, :, index_to_switch
            ] = right_val_super
            self.data_frame["RHO_ORDERS_RIGHT"][:, :, index_to_switch] = left_val_order

            self.data_frame["RHO_ORDERS_LEFT"][:, :, index_to_switch] = right_val_order

    def __repr__(self):
        """String representation for developer."""
        info = """
                DensityData shape key: (type of TE, windows, genes)
                DensityData left_order shape: {self.left_orders.shape}
                DensityData intra_order shape: {self.intra_orders.shape}
                DensityData right_order shape: {self.right_orders.shape}
                DensityData left_supers shape: {self.left_supers.shape}
                DensityData intra_supers shape: {self.intra_supers.shape}
                DensityData right_supers shape: {self.right_supers.shape}
            """
        return info.format(self=self)
