#!/usr/bin/env python3

"""
Unit test NormMatrix
"""

__author__ = "Scott Teresi"

import os
import pytest
import numpy as np

from transposon.normalization_matrix import NormMatrix
from transposon.gene_data import GeneData
from transposon.import_genes import import_genes


@pytest.fixture
def GeneData_test_obj():
    path_main = os.path.abspath(__file__)
    default_path = os.path.join(path_main, '../../', 'test_geneframe.tsv')
    # default_path = os.path.join(path_main, '../../', 'Cam_Genes.gtf')  # use
    # this if you want to use the full real dataset
    default_path = os.path.abspath(default_path)
    genes_input_dataframe = import_genes(default_path)
    sample_genome = GeneData(genes_input_dataframe, 'Mock_Camarosa')
    return sample_genome


@pytest.fixture
def Windows():
    return [x for x in range(500, 10001, 500)]


# NOTE expected test val below with a window of 500
# NOTE window value is a fixture, TODO candidate for future simplification
LEFT_DIVISORS = np.array([41, 501, 501, 501, 501, 501, 501, 501, 501, 501])
INTRA_DIVISORS = np.array([2356, 2423, 311, 298, 285, 2106, 3259, 1778, 1559, 1495])
RIGHT_DIVISORS = np.array([501, 501, 501, 501, 501, 501, 501, 501, 501, 501])


@pytest.mark.parametrize(
                        'true_divisors',
                        [LEFT_DIVISORS])
def test_divisor_left(GeneData_test_obj, Windows, true_divisors):
    divisor_frame = NormMatrix.divisor_left(GeneData_test_obj, Windows)
    assert np.array_equal(true_divisors, divisor_frame[0])  # MAGIC NUMBER 0 is
    # the index of the resultant dataframe that we need to compare against


@pytest.mark.parametrize(
                        'true_divisors',
                        [INTRA_DIVISORS])
def test_divisor_intra(GeneData_test_obj, Windows, true_divisors):
    divisor_frame = NormMatrix.divisor_intra(GeneData_test_obj, Windows)
    assert np.array_equal(true_divisors, divisor_frame[0])  # MAGIC NUMBER 0 is
    # the index of the resultant dataframe that we need to compare against


@pytest.mark.parametrize(
                        'true_divisors',
                        [RIGHT_DIVISORS])
def test_divisor_right(GeneData_test_obj, Windows, true_divisors):
    divisor_frame = NormMatrix.divisor_right(GeneData_test_obj, Windows)
    assert np.array_equal(true_divisors, divisor_frame[0])  # MAGIC NUMBER 0 is
    # the index of the resultant dataframe that we need to compare against


if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
