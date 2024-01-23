#!/usr/bin/env python3

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np

import os
import argparse
import logging
import coloredlogs

"""
- Provide helper reader function to read the cleaned data from disk
"""


def read_cleaned_homologs(path):
    """
    Read cleaned homolog data from disk
    """
    return pd.read_csv(path, sep="\t", header="infer")
