#!/usr/bin/env/python

"""
Compare TE Density values between orthologs
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd
import h5py
from collections import defaultdict
import time
import matplotlib.pyplot as plt

from graphs.graph import _split_wrt_chromosome
from transposon.gene_data import GeneData
from scripts.density_data import DensityData


def read_ortholog_table(strawberry_at_ortholog_table):
    dataframe = pd.read_csv(strawberry_at_ortholog_table, sep="\t", header="infer")
    return dataframe


def compare_two_chrom_test(orthologs, dd_obj1, dd_obj2, output_dir):
    """

    """
    cam_orth_genes = orthologs["Camarosa_Gene"].tolist()
    h4_orth_genes = orthologs["H4"].tolist()
    dd1_h4_indices_and_genes = {
        dd_gene: dd_index for dd_index, dd_gene in enumerate(dd_obj1.gene_list)
    }
    orthologs["H4_Indices"] = [
        dd1_h4_indices_and_genes.get(gene_name, None) for gene_name in h4_orth_genes
    ]

    dd2_cam_indices_and_genes = {
        dd_gene: dd_index for dd_index, dd_gene in enumerate(dd_obj2.gene_list)
    }

    orthologs["Camarosa_Indices"] = [
        dd2_cam_indices_and_genes.get(gene_name, None) for gene_name in cam_orth_genes
    ]
    # SAVE
    orthologs.to_csv(
        os.path.join(output_dir, "Test_Indices.tsv"), header=True, index=False, sep="\t"
    )

    # Now we want to compare values
    # So now just get the union for the two index columns that both have data
    # h4_to_cam_matches = orthologs.loc[orthologs_new["H4"].notnull()]
    h4_to_cam_matches = orthologs[
        orthologs[["Camarosa_Indices", "H4_Indices"]].notnull().all(1)
    ]  # the all command wants them both to be true
    # print(h4_to_cam_matches)
    h4_indices = h4_to_cam_matches["H4_Indices"].tolist()
    cam_indices = h4_to_cam_matches["Camarosa_Indices"].tolist()

    # But note that if I feed the indices to index the h5, they aren't in
    # order, but I lose the ability to compare values if I sort it, because I
    # lose the pair information

    differences = []
    for h4_val, cam_val in zip(h4_indices, cam_indices):
        # TODO verify that the transposon names are in the same order
        diff_val = (
            dd_obj1.data_frame["RHO_ORDERS_LEFT"][1, 1, h4_val]
            - dd_obj2.data_frame["RHO_ORDERS_LEFT"][1, 1, cam_val]
        )  # remember shape is type, window, gene
        differences.append(diff_val)
    print(dd_obj1.order_list)
    # print(dd_obj2.order_list)

    test_diff = differences

    print(max(test_diff))
    print(min(test_diff))
    x = test_diff
    num_bins = 15
    n, bins, patches = plt.hist(x, num_bins, facecolor="blue", alpha=0.5)
    plt.show()


def compare_two_chrom(orthologs, dd_obj1, dd_obj2, output_dir):
    # All entries of a camarosa and H4 gene pair, N.B some pairs will not be
    # completed with an Arabidopsis match
    orthologs_new = orthologs.drop(columns=["Del_Norte", "Arabidopsis_Gene"])
    h4_to_cam_matches = orthologs_new.loc[orthologs_new["H4"].notnull()]
    h4_cam_no_dupes = h4_to_cam_matches.drop_duplicates()
    camarosa_genes = h4_to_cam_matches["Camarosa_Gene"].tolist()
    # Genes repeat in the camarosa list, but they won't repeat in DD object

    # y = []
    # for i, item in enumerate(dd_obj2.gene_list):
    # if item in h4_cam_no_dupes["Camarosa_Gene"].tolist():
    # y.append(i)
    # else:
    # y.append(None)

    # print(len(y))
    # print(y)

    # NOTE it seems y is what I want
    y = []
    for gene in camarosa_genes:  # loop over camarosa genes in orthology table
        for i, d_gene in enumerate(dd_obj2.gene_list):  # look over dd obj
            if gene == d_gene:  # if gene in O table is equiv to dd obj gene
                y.append(i)  # store the gene
        else:
            y.append(None)
    print(len(y))
    print(y)
    print(len(h4_cam_no_dupes["Camarosa_Gene"].tolist()))

    # h4_cam_no_dupes["Camarosa_Gene_DD_Idx"] = y
    # print(h4_cam_no_dupes.head())
    raise ValueError

    # gene_dict = h4_cam_no_dupes.groupby("Camarosa_Gene")["H4"].apply(list).to_dict()
    # gene_frame = h4_cam_no_dupes.groupby("Camarosa_Gene")["H4"].apply(list).to_frame()
    # A camarosa gene can only be represented once and it can have multiple
    # H4 join pointing to it as a list in the divtionary
    # gene_frame.reset_index(level=0, inplace=True)

    # gene_frame.to_csv(
    # os.path.join(output_dir, "test.tsv"), sep="\t", header=True, index=False,
    # )
    # cam_genes = gene_frame["Camarosa_Gene"].tolist()
    # h4_genes = gene_frame["H4"].tolist()
    # h4_genes = list(
    # set([item for sublist in h4_genes for item in sublist])
    # )  # flatten a list of lists

    # I need to find the index of the Camarosa genes of my gene_frame in the
    # density data

    # cam_list = []
    # for count, d_gene in enumerate(dd_obj2.gene_list):
    # for gene in cam_genes:
    # if gene == d_gene:
    # cam_list.append(count)
    # cam_list.append(None)
    # print(cam_list)
    # print(len(dd_obj2.gene_list))
    # print(len(cam_genes))
    # print(len(cam_list))

    raise ValueError
    idx_of_cam_genes_in_density_w_pairs = [
        i for i, item in enumerate(dd_obj2.gene_list) if item in cam_genes
    ]
    idx_of_h4_genes_in_density_w_pairs = [
        i for i, item in enumerate(dd_obj1.gene_list) if item in h4_genes
    ]
    print(len(idx_of_cam_genes_in_density_w_pairs))
    print(len(idx_of_h4_genes_in_density_w_pairs))

    print(dd_obj1.left_orders.shape)
    print(dd_obj2.left_orders.shape)

    # print(cam_indices)
    # print(np.array_equal(x, cam_indices))

    # my_bool = np.in1d(dd_obj2.gene_list, gene_frame["Camarosa_Gene"].tolist())
    # mask = np.isin(gene_frame["Camarosa_Gene"].tolist(), dd_obj2.gene_list)
    # print(np.where(mask, dd_obj2.gene_list))
    # print(np.where(gene_frame["Camarosa_Gene"].tolist() in dd_obj2.gene_list))
    # gene_frame["Cam_Gene_DD_Index"] = np.where(
    # np.in1d(dd_obj2.gene_list, gene_frame["Camarosa_Gene"])
    # )[
    # 0
    # ]  # MAGIC

    # print(gene_frame)
    raise ValueError

    # Make a copy of the ortholog set that is the union between the two boolean
    # frames
    # I have a list of genes that I have data for, for both density data
    # objects
    # union_orth = h4_to_cam_matches[(h4_bool) & (camarosa_bool)]

    # Write the subsetted ortholog set to disk for testing purposes
    # union_orth.to_csv(
    # os.path.join(output_dir, "Fvb1_andFvb1_4_matches_to_DD_data.tsv"),
    # sep="\t",
    # header=True,
    # index=False,
    # )

    list_of_d1_densities = []
    list_of_d2_densities = []

    # I now need to subset the DD objects, removing the data entries for genes
    # that are NOT in the list of present genes in the orth union
    # present_cam_genes_list = union_orth.Camarosa_Gene.to_list()
    # present_h4_genes_list = union_orth.H4.to_list()
    # print(len(present_h4_genes_list))
    # print(len(present_cam_genes_list))
    # print(len(set(present_h4_genes_list)))
    # print(len(set(present_cam_genes_list)))

    # SO it appears that I am having an issue because there is a case where a
    # Camarosa gene is being paired with multiple H4 genes

    # Show where A is found in B, and get the indices of those True values in
    # A
    # The lists, which were of same length, may not all be found in the density
    # data, because those lists were made from ALL genes. Here when looking at
    # DD data we are looking at 1 chromosome at a time. SO to enforce a 1-1
    # ratio for the differences in density values we need to make sure the....
    # Are the numbers not equal here because of translocations and such?
    # So we were looking at the en

    cam_indices = np.where(np.in1d(dd_obj2.gene_list, np.asarray(cam_genes)))[
        0
    ]  # MAGIC

    cam_indices = np.where(
        np.in1d(dd_obj2.gene_list, np.asarray(present_cam_genes_list))
    )[
        0
    ]  # MAGIC
    cam_indices = cam_indices.tolist()
    h4_indices = np.where(
        np.in1d(dd_obj1.gene_list, np.asarray(present_h4_genes_list))
    )[
        0
    ]  # MAGIC
    h4_indices = h4_indices.tolist()
    # Why are the lengths of the indice matrices not the same? They should have
    # the same hits
    print(len(set(present_h4_genes_list)))
    print(len(set(present_cam_genes_list)))
    print(len(set(dd_obj1.gene_list)))
    print(len(set(dd_obj2.gene_list)))

    for gene in present_cam_genes_list:
        if gene not in dd_obj2.gene_list:
            print("Camarosa")
            print(gene)
    for gene in present_h4_genes_list:
        if gene not in dd_obj1.gene_list:
            print("H4")
            print(gene)
    raise ValueError

    # NOTE
    # I need to get the indices of the PAIRS and calculate the difference
    # between the PAIRS. Maybe make an intermediate pandas set? Or will that
    # need to be done in H5? Can I just make a list of values from the
    # differences?

    # Now I should be able to index the DD object by the values and plot those,
    # subsetting the total number of data points
    # this will of course need to be mirrored in the other DD object.

    print(np.unique(dd_obj2.left_orders.shape))
    dd_obj2.left_orders = dd_obj2.data_frame["RHO_ORDERS_LEFT"][:, :, cam_indices]
    dd_obj1.left_orders = dd_obj1.data_frame["RHO_ORDERS_LEFT"][:, :, h4_indices]
    print(np.unique(dd_obj2.left_orders.shape))
    print(np.unique(dd_obj1.left_orders.shape))


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    data_dir = os.path.abspath(os.path.join(dir_main, "../../", "Domestication_Data"))
    parser = argparse.ArgumentParser(
        description="compare density values between orthologs"
    )

    parser.add_argument(
        "ortholog_input_file", type=str, help="parent path to ortholog file",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=data_dir,
        help="parent directory to output results",
    )
    args = parser.parse_args()
    args.ortholog_input_file = os.path.abspath(args.ortholog_input_file)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    test_file_h4 = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/H4/5KB_Up/H4_Fvb1.h5"
    test_file_cam = "/home/scott/Documents/Uni/Research/Projects/TE_Data/finalized_data/Fvb1_to_Fvb1_4_Compare_5KB_Up/Camarosa_Fvb1-4.h5"
    all_genes_h4 = pd.read_csv(
        "/home/scott/Documents/Uni/Research/Projects/TE_Data/filtered_input_data/Cleaned_H4_Genes.tsv",
        sep="\t",
        header="infer",
        dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
        index_col="Gene_Name",  # this is crucial
    )
    all_genes_cam = pd.read_csv(
        "/home/scott/Documents/Uni/Research/Projects/TE_Data/filtered_input_data/Cleaned_Camarosa_Genes.tsv",
        sep="\t",
        header="infer",
        dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
        index_col="Gene_Name",  # this is crucial
    )
    h5_file_list_h4 = [test_file_h4]
    h5_file_list_cam = [test_file_cam]

    genes_split = _split_wrt_chromosome(all_genes_h4)
    genes_split = [
        GeneData(geneframe, ("H4_" + geneframe.Chromosome.unique()[0]))
        for geneframe in genes_split
    ]
    chromosomes = [chromosome for chromosome in genes_split]  # gene data
    chromosomes = chromosomes[0]  # MAGIC select the right chromosome for our
    # test set here
    dd_objects_h4 = DensityData(h5_file_list_h4[0], chromosomes, logger)  # MAGIC

    # Load Camarosa
    genes_split = _split_wrt_chromosome(all_genes_cam)
    genes_split = [
        GeneData(geneframe, ("Cam_" + geneframe.Chromosome.unique()[0]))
        for geneframe in genes_split
    ]
    chromosomes = [chromosome for chromosome in genes_split]  # gene data
    chromosomes = chromosomes[3]  # MAGIC select the right chromosome for our
    # test set here

    # NOTE OLD method
    # dd_objects_cam = [
    # DensityData(h5_file, chromosome, logger)
    # for h5_file, chromosome, in zip(h5_file_list_cam, chromosomes)
    # ]
    dd_objects_cam = DensityData(h5_file_list_cam[0], chromosomes, logger)  # MAGIC

    # Execute
    logger.info("Reading ortholog input file: %s" % (args.ortholog_input_file))
    orthologs = read_ortholog_table(args.ortholog_input_file)

    compare_two_chrom_test(orthologs, dd_objects_h4, dd_objects_cam, args.output_dir)
