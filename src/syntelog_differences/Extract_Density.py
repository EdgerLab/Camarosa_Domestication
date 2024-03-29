"""Parse the density files for Strawberry purposes"""


from transposon.import_filtered_genes import import_filtered_genes
from transposon.gene_data import GeneData
from transposon.density_data import DensityData
from transposon.density_utils import (
    add_hdf5_indices_to_gene_data_from_list_hdf5,
    add_te_vals_to_gene_info_pandas_from_list_hdf5,
    add_te_vals_to_gene_info_pandas,
    get_specific_slice,
    add_hdf5_indices_to_gene_data,
    info_of_gene,
)

import os


class Strawberry_Specific_Density:
    def __init__(
        self,
        gene_frame_with_indices,
        processed_density_data,
        genome_name,
        major_group,
        te_type,
        direction,
        window,
    ):
        self.gene_frame_with_indices = gene_frame_with_indices
        self.processed_density_data = processed_density_data
        self.genome_name = genome_name
        self.major_group = major_group
        self.te_type = te_type
        self.direction = direction
        self.window = window
        self.te_window_direction_str = f"{self.te_type}_{self.window}_{self.direction}"
        self.table = self.make_table_for_te_type_and_direction()

    def make_table_for_te_type_and_direction(self):
        """Make a table for a specific TE type and direction"""
        table = add_te_vals_to_gene_info_pandas_from_list_hdf5(
            self.gene_frame_with_indices,
            self.processed_density_data,
            self.major_group,
            self.te_type,
            self.direction,
            self.window,
        )

        table.rename(
            columns={
                self.te_window_direction_str: f"{self.genome_name}_{self.te_window_direction_str}",
                "Gene_Name": f"{self.genome_name}_Gene",
            },
            inplace=True,
        )

        # Drop uninformative columns
        table.drop(
            columns=[
                "Feature",
                "Index_Val",
            ],
            inplace=True,
        )
        return table

    def save_table_to_disk(self, output_dir):
        """Save the table to disk"""
        outfile = os.path.join(
            output_dir, f"{self.genome_name}_{self.te_window_direction_str}.tsv"
        )
        self.table.to_csv(outfile, header=True, index=False, sep="\t")
