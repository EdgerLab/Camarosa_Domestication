#!/usr/bin/env python3

__author__ = "Scott Teresi"

import os
import pandas as pd
from import_homologs import *


def verify_BLAST_cache(
    raw_blast_data_file, genome_name, filtered_blast_data_file, logger
):
    """Determine whether or not previously saved BLAST data exists on disk"""

    if os.path.exists(filtered_blast_data_file):
        raw_file_creation_time = os.path.getmtime(raw_blast_data_file)
        filtered_file_creation_time = os.path.getmtime(filtered_blast_data_file)

        if raw_file_creation_time > filtered_file_creation_time:
            logger.info("Cache is too old for `%s`" % filtered_file_creation_time)
            logger.info("Create: %s" % filtered_blast_data_file)
            homolog_data = import_homologs(homolog_input_file, genome_name)

        else:
            logger.info(
                "BLAST cache exists, reading cache `%s`" % filtered_blast_data_file
            )
            homolog_data = pd.read_csv(
                filtered_blast_data_file, sep="\t", header="infer"
            )
    else:
        logger.info("BLAST cache DNE `%s`" % filtered_blast_data_file)
        logger.info("Create: %s" % filtered_blast_data_file)
        homolog_data = import_homologs(homolog_input_file, genome_name)

    return homolog_data
