#!/usr/bin/env python3

"""
Unit test NormMatrix
"""

__author__ = "Scott Teresi"

import os
import pytest
import numpy as np
import pandas as pd
from io import StringIO

from transposon.normalization_matrix import NormMatrix
from transposon.gene_data import GeneData
from transposon.import_genes import import_genes


@pytest.fixture
def GeneData_test_obj():
    gene_file = 'tests/input_data/Test_Genes_NormMatrix.tsv'
    gene_pandas = pd.read_csv(gene_file, header='infer', sep='\t',
                              dtype={'Start': 'float32', 'Stop': 'float32',
                              'Length': 'float32'}, index_col='Gene_Name')
    sample_genome = GeneData(gene_pandas, 'Mock_Camarosa')
    return sample_genome


LEFT_DIVISORS = {500: np.array([41, 501, 501, 501, 501, 501, 501, 501, 501,
                                501, 501, 501, 501, 501, 501, 501, 501, 501,
                                501, 501, 501, 501, 501, 501, 501, 501, 501,
                                501, 501, 501, 501, 501, 501, 501, 501, 501,
                                501, 501, 501, 501]),
                 3500: np.array([41, 3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501,
                                3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501,
                                3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501,
                                3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501,
                                3501, 3501, 3501, 3501]),
                 7000: np.array([41, 5556, 7001, 7001, 7001, 7001, 7001, 7001, 7001,
                                7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001,
                                7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001,
                                7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001,
                                7001, 7001, 7001, 7001])}

# Should always just be the gene length, so even though we test over multiple
# windows the normalization value should stay the same
INTRA_DIVISORS = {500: np.array([2356, 2423, 311, 298, 285, 2106, 3259, 1778,
                                 1559, 1495, 363, 7434, 1550, 1900, 12502, 9274, 2802,
                                 189, 522, 3312, 4469, 3301, 4558, 2193, 5681, 204,
                                 1466, 1883, 4476, 7998, 1179, 4160, 12075, 312, 1073,
                                 1751, 1950, 312, 4224, 264], dtype=float),
                 3500: np.array([2356, 2423, 311, 298, 285, 2106, 3259, 1778,
                                 1559, 1495, 363, 7434, 1550, 1900, 12502, 9274, 2802,
                                 189, 522, 3312, 4469, 3301, 4558, 2193, 5681, 204,
                                 1466, 1883, 4476, 7998, 1179, 4160, 12075, 312, 1073,
                                 1751, 1950, 312, 4224, 264], dtype=float),
                 7000: np.array([2356, 2423, 311, 298, 285, 2106, 3259, 1778,
                                 1559, 1495, 363, 7434, 1550, 1900, 12502, 9274, 2802,
                                 189, 522, 3312, 4469, 3301, 4558, 2193, 5681, 204,
                                 1466, 1883, 4476, 7998, 1179, 4160, 12075, 312, 1073,
                                 1751, 1950, 312, 4224, 264], dtype=float)}

RIGHT_DIVISORS = {500: np.array([501, 501, 501, 501, 501, 501, 501, 501, 501,
                                501, 501, 501, 501, 501, 501, 501, 501, 501,
                                501, 501, 501, 501, 501, 501, 501, 501, 501,
                                501, 501, 501, 501, 501, 501, 501, 501, 501,
                                501, 501, 501, 501]),
                 3500: np.array([3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501,
                                3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501,
                                3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501,
                                3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501, 3501,
                                3501, 3501, 3501, 3501]),
                 7000: np.array([7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001,
                                7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001,
                                7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001,
                                7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001, 7001,
                                7001, 7001, 7001, 7001])}


@pytest.mark.parametrize(
                        'true_divisors',
                        [LEFT_DIVISORS])
def test_divisor_left(GeneData_test_obj, true_divisors):
    for window, true_val_array in true_divisors.items():
        window = [window]
        divisor_frame = NormMatrix.divisor_left(GeneData_test_obj, window)
        assert np.array_equal(true_val_array, divisor_frame[0])  # MAGIC number
        # Have to use magic number because divisor frame will essentially be a
        # list of lists during production. But for the sake of tests, which we
        # are doing one at a time, we need to index it for comparisons.


@pytest.mark.parametrize(
                        'true_divisors',
                        [INTRA_DIVISORS])
def test_divisor_intra(GeneData_test_obj, true_divisors):
    for window, true_val_array in true_divisors.items():
        window = [window]
        divisor_frame = NormMatrix.divisor_intra(GeneData_test_obj, window)
        assert np.array_equal(true_val_array, divisor_frame[0])  # MAGIC number


@pytest.mark.parametrize(
                        'true_divisors',
                        [RIGHT_DIVISORS])
def test_divisor_right(GeneData_test_obj, true_divisors):
    for window, true_val_array in true_divisors.items():
        window = [window]
        divisor_frame = NormMatrix.divisor_right(GeneData_test_obj, window)
        assert np.array_equal(true_val_array, divisor_frame[0])  # MAGIC number
        # Have to use magic number because divisor frame will essentially be a
        # list of lists during production. But for the sake of tests, which we
        # are doing one at a time, we need to index it for comparisons.


if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
