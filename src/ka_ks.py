__author__ = "Scott Teresi"


"""
TODO
"""

import pandas as pd
import numpy as np
import argparse
import os
import logging
import coloredlogs


def read_raw_h4_dn_ka_ks_table():
    raise NotImplementedError
    table = filter_h4_dn(table)


def filter_h4_dn(table):
    raise NotImplementedError


def read_raw_h4_rr_ka_ks_table():
    raise NotImplementedError
    table = filter_h4_rr(table)


def filter_h4_rr(table):
    raise NotImplementedError


def read_filtered_h4_dn_ka_ks_table():
    raise NotImplementedError


def read_filtered_h4_rr_ka_ks_table():
    raise NotImplementedError


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TODO")
    parser.add_argument(
        "gene_input_file", type=str, help="Parent path of gene annotation file"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.gene_input_file = os.path.abspath(args.gene_input_file)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # ------------------------------------------------------------------
