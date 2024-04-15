#!/usr/bin/env/python
"""
TODO Description
"""
__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd
from collections import namedtuple

if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "todo_file",
        type=str,
        help="TODO",
    )

    parser.add_argument("strawberry_ortholog_table", type=str, help="TODO")

    parser.add_argument(
        "output_dir",
        type=str,
        help="parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.todo_file = os.path.abspath(args.todo_file)
    args.strawberry_ortholog_table = os.path.abspath(args.strawberry_ortholog_table)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -------------------------------------------------------------------------

    file = pd.read_csv(args.todo_file, sep="\t", header="infer")
    print(file)
