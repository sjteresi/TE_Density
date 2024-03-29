#!/usr/bin/env/python

"""
Unit test DensityData
"""

__author__ = "Scott Teresi"

import h5py
import numpy as np
import pytest
import logging
import coloredlogs

from transposon import write_vlen_str_h5py
from transposon.density_data import DensityData
from transposon.import_filtered_genes import import_filtered_genes
from transposon.gene_data import GeneData

LOGGER = logging.getLogger(__name__)
coloredlogs.install(level=logging.DEBUG)
TEST_GENE_FILE = "tests/input_data/Test_Genes_DensityData.tsv"


@pytest.fixture
def chromosomes():
    """Return a list of chromosomes for h5py file"""
    return ["TestChrom1"]


@pytest.fixture
def gene_names():
    """Yield a list of genes"""
    gene_names = [
        "dummy1",
        "dummy2",
        "dummy3",
        "dummy4",
        "dummy5",
        "dummy6",
        "dummy7",
        "dummy8",
    ]
    return gene_names


@pytest.fixture
def dummy_gene_data():
    """Yield a GeneData for testing"""
    gene_pandas = import_filtered_genes(TEST_GENE_FILE, LOGGER)
    return GeneData(gene_pandas, "dummy_genome_name", LOGGER)


@pytest.fixture
def order_names():
    """Return a list of order types"""
    order_names = ["LTR", "TIR"]
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
    """Return an integer which is the total number of genes"""
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
    # All this does is make the values different, not the shape
    matrix_num = total_orders * total_windows * total_genes * 2
    arr = np.arange(48, matrix_num).reshape((total_orders, total_windows, total_genes))
    return arr


@pytest.fixture
def density_data_test_obj_swap_vals(
    chromosomes,
    dummy_gene_data,
    order_names,
    super_names,
    windows,
    rho_o_left,
    rho_o_intra,
    rho_o_right,
):
    """Create a test object for DensityData, reads from file"""
    # TODO set path to be a variable not hard-coded, plus it repeats further
    # down
    f = h5py.File(
        "tests/input_data/test_swap_file.h5",
        "w",
    )
    gene_names = list(dummy_gene_data.names)
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
    f.create_dataset("RHO_SUPERFAMILIES_RIGHT", data=rho_o_right)
    f.close()
    return DensityData(
        "tests/input_data/test_swap_file.h5",
        dummy_gene_data,
        LOGGER,
        sense_swap=False,
    )


def test_swap_density_vals(density_data_test_obj_swap_vals):
    """Test whether or not left and right density values are swapped
    correctly"""
    # Set values for easy testing and checking
    # Shape of left/right orders is (no. of TE types, window, no. genes).
    # Here we do index 1, which corresponds to Gene2

    density_data_test_obj_swap_vals.left_orders[:, :, 1] = 100
    density_data_test_obj_swap_vals.right_orders[:, :, 1] = 200

    # Call the value swapper, manually, because it was set to NOT do it upon
    # initialization
    density_data_test_obj_swap_vals._swap_strand_vals(["dummy2"])

    # Check the values
    assert np.array_equal(
        density_data_test_obj_swap_vals.left_orders[:, :, 1], np.full((2, 3), 200)
    )
    assert np.array_equal(
        density_data_test_obj_swap_vals.right_orders[:, :, 1], np.full((2, 3), 100)
    )


if __name__ == "__main__":
    pytest.main(["-s", __file__])  # for convenience
