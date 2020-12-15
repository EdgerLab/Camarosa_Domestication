#!/usr/bin/env/python

"""
Unit test DensityData
"""

__author__ = "Scott Teresi"

import h5py
import numpy as np
import pytest

from scripts.density_data import DensityData


def write_vlen_str_h5py(h5file, strings, dataset_key):
    """Write to an H5 File an iterable of variable length unicode

    Args:
        h5file (h5py.File): opened destination file
        strings (iterable(str)): string list
        dataset_key(str): name of data set to write

    """
    vlen = h5py.special_dtype(vlen=str)
    n_strings = sum(1 for s in strings)
    dset = h5file.create_dataset(dataset_key, (n_strings,), dtype=vlen)
    dset[:] = strings


@pytest.fixture
def chromosomes():
    """Return a list of chromosomes for h5py file"""
    return ["TestChrom1", "TestChrom2", "TestChrom3"]


@pytest.fixture
def gene_names():
    """Yield a list of genes"""
    gene_names = [
        "Gene1",
        "Gene2",
        "Gene3",
        "Gene4",
        "Gene5",
        "Gene6",
        "Gene7",
        "Gene8",
    ]
    return gene_names


@pytest.fixture
def order_names():
    """Return a list of order types"""
    order_names = ["LTR", "DNA"]
    return order_names


@pytest.fixture
def super_names():
    """Return a list of super names"""
    super_names = ["Copia", "Gypsy", "HAT"]
    return super_names


@pytest.fixture
def windows():
    """Return a list of window values"""
    windows = ["500", "1000", "1500"]
    return windows


@pytest.fixture
def total_orders(order_names):
    """Return an integer which is the total number of orders"""
    return sum(1 for order in order_names)


@pytest.fixture
def total_windows(windows):
    """Return an integer which is the total number of windows"""
    return sum(1 for window in windows)


@pytest.fixture
def total_genes(gene_names):
    """Return an integer which is the total number of windows"""
    return sum(1 for gene in gene_names)


@pytest.fixture
def rho_o_left(total_orders, total_windows, total_genes):
    """Return an array of order left values"""
    # NOTE currently it will perform np.arange(48) and reshape to (2,3,8)
    matrix_num = total_orders * total_windows * total_genes
    arr = np.arange(matrix_num).reshape((total_orders, total_windows, total_genes))
    return arr


@pytest.fixture
def rho_o_intra(total_orders, total_genes):
    """Return an array of order left values"""
    # NOTE currently it will perform np.arange(48)
    matrix_num = total_orders * total_genes
    arr = np.arange(matrix_num).reshape((total_orders, 1, total_genes))
    return arr


@pytest.fixture
def rho_o_right(total_orders, total_windows, total_genes):
    """Return an array of order left values"""
    # NOTE currently it will perform np.arange(48) *2 starting from 48
    # so that we may differentiate from rho o left, and reshape to (2,3,8)
    matrix_num = total_orders * total_windows * total_genes * 2
    arr = np.arange(48, matrix_num).reshape((total_orders, total_windows, total_genes))
    return arr


@pytest.fixture
def density_data_test_obj(
    chromosomes,
    gene_names,
    order_names,
    super_names,
    windows,
    rho_o_left,
    rho_o_intra,
    rho_o_right,
):
    """Create a test object for DensityData, reads from file"""
    f = h5py.File(
        "/home/scott/Documents/Uni/Research/Projects/Camarosa_Domestication/tests/input_data/testfile.hdf5",
        "w",
    )
    write_vlen_str_h5py(f, chromosomes, "CHROMOSOME_ID")
    write_vlen_str_h5py(f, gene_names, "GENE_NAMES")
    write_vlen_str_h5py(f, order_names, "ORDER_NAMES")
    write_vlen_str_h5py(f, super_names, "SUPERFAMILY_NAMES")
    write_vlen_str_h5py(f, windows, "WINDOWS")

    dset = f.create_dataset("RHO_ORDERS_LEFT", data=rho_o_left)
    dset = f.create_dataset("RHO_ORDERS_INTRA", data=rho_o_intra)
    dset = f.create_dataset("RHO_ORDERS_RIGHT", data=rho_o_right)

    # NB just re-doing the values for the supers because not testing supers
    dset = f.create_dataset("RHO_SUPERFAMILIES_LEFT", data=rho_o_left)
    dset = f.create_dataset("RHO_SUPERFAMILIES_INTRA", data=rho_o_intra)
    dset = f.create_dataset("RHO_SUPERFAMILIES_RIGHT", data=rho_o_right)
    f.close()
    return DensityData(
        "/home/scott/Documents/Uni/Research/Projects/Camarosa_Domestication/tests/input_data/testfile.hdf5",
        ["Gene1"],  # dummy switch,
    )


def test_swap_density_vals(density_data_test_obj):
    """Test whether or not left and right density values are swapped
    correctly"""
    # Set values for easy testing and checking
    # Shape of left/right orders is (no. of TE types, window, no. genes).
    # Here we do index 1, which corresponds to Gene2
    density_data_test_obj.left_orders[:, :, 1] = 100
    density_data_test_obj.right_orders[:, :, 1] = 200

    # Call the value swapper
    density_data_test_obj._swap_strand_vals(["Gene2"])

    # Check the values
    assert np.array_equal(
        density_data_test_obj.left_orders[:, :, 1], np.full((2, 3), 200)
    )
    assert np.array_equal(
        density_data_test_obj.right_orders[:, :, 1], np.full((2, 3), 100)
    )


if __name__ == "__main__":
    pytest.main(["-s", __file__])  # for convenience
