#!/usr/bin/env python3

"""
Master code file. Control filtration of syntelog data and homolog data. Manage
merging of the two dataframes and orchestrate all file commands.
"""

__author__ = "Scott Teresi"

import argparse
import os

import logging
import coloredlogs
import pandas as pd


from import_syntelogs import import_syntelogs
from import_syntelogs import Syntelog_Data
from import_homologs import import_homologs
from import_homologs import Homolog_Data
from merge_homo_synt import merge_homo_synt
from merge_homo_synt import Merged_Data
from verify_cache import verify_BLAST_cache


def process_CAM(
    syntelog_input_file, homolog_input_file, genome_name, data_output_path,
):
    # Import the synteny data from raw file
    logger.info("Importing raw syntelogs: %s" % syntelog_input_file)
    syntelogs = import_syntelogs(syntelog_input_file, genome_name)

    # Wrap the data
    logger.debug("Wrapping Syntelog_Data...")
    instance_Syntelog_Data = Syntelog_Data(syntelogs)
    file_to_save = os.path.join(
        data_output_path, ("Synteny_AT_" + genome_name + ".tsv")
    )
    logger.info("Writing syntelog data to disk: %s" % file_to_save)
    instance_Syntelog_Data.save_to_disk(file_to_save)

    # Import the homology (BLAST) data from raw file
    logger.info("Importing raw BLAST results: %s" % homolog_input_file)
    filtered_blast_data_file = os.path.join(
        data_output_path, str("Homologs_AT_" + genome_name + ".tsv")
    )

    homologs = verify_BLAST_cache(
        homolog_input_file, genome_name, filtered_blast_data_file, logger
    )

    # Wrap the data
    logger.debug("Wrapping Homolog_Data...")
    instance_Homolog_Data = Homolog_Data(homologs)
    file_to_save = os.path.join(
        data_output_path, ("Homologs_AT_" + genome_name + ".tsv")
    )
    logger.info("Writing homolog data to disk: %s" % file_to_save)
    instance_Homolog_Data.save_to_disk(file_to_save)

    # Merge the synteny and homology data
    logger.debug("Merging the data...")
    merged_all = merge_homo_synt(
        instance_Syntelog_Data, instance_Homolog_Data, genome_name
    )

    # Wrap the data
    logger.debug("Wrapping the merged data...")
    instance_Merged_Data = Merged_Data(merged_all)
    # Save to disk
    file_to_save = os.path.join(
        data_output_path, str(genome_name + "_merged_homo_syn.tsv")
    )
    logger.info("Writing merged data to disk: %s" % file_to_save)
    instance_Merged_Data.save_to_disk(file_to_save)
    synteny_count, blast_count = instance_Merged_Data.total_count_summary
    # NOTE look at the above if you want the numbers of BLAST vs Synteny
    return instance_Merged_Data


def process_H4(syntelog_input_file, genome_name, data_output_path):
    # Import the synteny data from raw file
    logger.info("Importing raw syntelogs: %s" % syntelog_input_file)
    syntelogs = import_syntelogs(syntelog_input_file, genome_name)

    # Wrap the data
    logger.debug("Wrapping Syntelog_Data...")
    instance_Syntelog_Data = Syntelog_Data(syntelogs)
    file_to_save = os.path.join(
        data_output_path, ("Synteny_CAM_" + genome_name + ".tsv")
    )
    logger.info("Writing syntelog data to disk: %s" % file_to_save)
    instance_Syntelog_Data.save_to_disk(file_to_save)
    # NOTE look at the above if you want the numbers of BLAST vs Synteny
    return instance_Syntelog_Data


def make_table(camarosa_merged, h4_synteny):
    # Remove the point of origin column for H4 because it is always Synteny
    h4_synteny.dataframe.drop("Point_of_Origin", inplace=True, axis=1)
    merged_all = camarosa_merged.dataframe.merge(
        h4_synteny.dataframe, on="Camarosa_Gene", how="outer"
    ).fillna("NA")
    merged_all.to_csv("test_total_merge.tsv", sep="\t", header=True, index=False)


if __name__ == "__main__":
    """Command line interface to link syntelogs together."""

    parser = argparse.ArgumentParser(description="Filter syntelogs")
    path_main = os.path.abspath(__file__)
    parser.add_argument(
        "cam_syntelog_input_file",
        type=str,
        help="parent path of camarosa syntelog file",
    )

    parser.add_argument(
        "cam_homolog_input_file", type=str, help="parent path of homolog file"
    )
    parser.add_argument(
        "HFour_syntelog_input_file",
        type=str,
        help="parent path of H4 (vesca) syntelog file",
    )

    parser.add_argument(
        "--output_directory",
        type=str,
        help="parent path of output directory",
        default=os.path.join(path_main, "../../../../Domestication_Data"),
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.cam_syntelog_input_file = os.path.abspath(args.cam_syntelog_input_file)
    args.cam_homolog_input_file = os.path.abspath(args.cam_homolog_input_file)
    args.HFour_syntelog_input_file = os.path.abspath(args.HFour_syntelog_input_file)
    args.output_directory = os.path.abspath(args.output_directory)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Process Camarosa
    logger.info("Starting filtration...")
    camarosa_merged = process_CAM(
        args.cam_syntelog_input_file,
        args.cam_homolog_input_file,
        "Camarosa",
        args.output_directory,
    )
    # Process H4
    h4_synteny = process_H4(args.HFour_syntelog_input_file, "H4", args.output_directory)
    make_table(camarosa_merged, h4_synteny)
