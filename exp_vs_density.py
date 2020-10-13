#!/user/bin/env/python3

"""

"""

__author__ = "Scott Teresi"

import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def exp_vs_density(window, te_type, direction, output_dir):
    """
    Upstream or Downstream:
        X-Axis: Density in a specific window from 0 to 1
        Y-Axis: Gene expression (TPM)

    Args:
        window (int): Integer representing the window for the analysis
        te_type (str): String of the TE name for the analysis
        direction (str): Upstream or downstream?
        output_dir (str): Path representing folder to save to

    """
    # TODO
    # Raise ValueError if direction does not fit within the categories.
    plt.plot(x, y)
    plt.title(os.path.join(output_dir(str(window)+ str(te_type) +'_Density_vs_Exp.png')))
    plt.xlabel('Transposon Density')

    # NOTE are we using log base 2?
    plt.ylabel('Gene Expression (Log(2) TPM')

if __name__ = '__main__':
    parser = argparse.ArgumentParser()
    path_main = os.path.abspath(__file__)
    raise NotImplementedError
    #parser.add_
