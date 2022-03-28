#!/usr/bin/env python3

"""
Unit test GeneDatum
"""

__author__ = "Scott Teresi"

import os
import pytest
import numpy as np
import pandas as pd
from io import StringIO

from transposon.gene_data import GeneData
from transposon.gene_datum import GeneDatum

# pytestmark = pytest.mark.skip  # skip all tests in this file

LEFT_DIVISORS = {
    500: np.array(
        [
            41,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
        ]
    ),
    3500: np.array(
        [
            41,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
        ]
    ),
    7000: np.array(
        [
            41,
            5556,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
        ]
    ),
}


INTRA_DIVISORS = {
    None: np.array(
        [
            2356,
            2423,
            311,
            298,
            285,
            2106,
            3259,
            1778,
            1559,
            1495,
            363,
            7434,
            1550,
            1900,
            12502,
            9274,
            2802,
            189,
            522,
            3312,
            4469,
            3301,
            4558,
            2193,
            5681,
            204,
            1466,
            1883,
            4476,
            7998,
            1179,
            4160,
            12075,
            312,
            1073,
            1751,
            1950,
            312,
            4224,
            264,
        ],
        dtype=float,
    ),
    None: np.array(
        [
            2356,
            2423,
            311,
            298,
            285,
            2106,
            3259,
            1778,
            1559,
            1495,
            363,
            7434,
            1550,
            1900,
            12502,
            9274,
            2802,
            189,
            522,
            3312,
            4469,
            3301,
            4558,
            2193,
            5681,
            204,
            1466,
            1883,
            4476,
            7998,
            1179,
            4160,
            12075,
            312,
            1073,
            1751,
            1950,
            312,
            4224,
            264,
        ],
        dtype=float,
    ),
    None: np.array(
        [
            2356,
            2423,
            311,
            298,
            285,
            2106,
            3259,
            1778,
            1559,
            1495,
            363,
            7434,
            1550,
            1900,
            12502,
            9274,
            2802,
            189,
            522,
            3312,
            4469,
            3301,
            4558,
            2193,
            5681,
            204,
            1466,
            1883,
            4476,
            7998,
            1179,
            4160,
            12075,
            312,
            1073,
            1751,
            1950,
            312,
            4224,
            264,
        ],
        dtype=float,
    ),
}

INTRA_INVALID_WINDOW_DIVISORS = {
    500: np.array(
        [
            2356,
            2423,
            311,
            298,
            285,
            2106,
            3259,
            1778,
            1559,
            1495,
            363,
            7434,
            1550,
            1900,
            12502,
            9274,
            2802,
            189,
            522,
            3312,
            4469,
            3301,
            4558,
            2193,
            5681,
            204,
            1466,
            1883,
            4476,
            7998,
            1179,
            4160,
            12075,
            312,
            1073,
            1751,
            1950,
            312,
            4224,
            264,
        ],
        dtype=float,
    ),
    1000: np.array(
        [
            2356,
            2423,
            311,
            298,
            285,
            2106,
            3259,
            1778,
            1559,
            1495,
            363,
            7434,
            1550,
            1900,
            12502,
            9274,
            2802,
            189,
            522,
            3312,
            4469,
            3301,
            4558,
            2193,
            5681,
            204,
            1466,
            1883,
            4476,
            7998,
            1179,
            4160,
            12075,
            312,
            1073,
            1751,
            1950,
            312,
            4224,
            264,
        ],
        dtype=float,
    ),
    8000: np.array(
        [
            2356,
            2423,
            311,
            298,
            285,
            2106,
            3259,
            1778,
            1559,
            1495,
            363,
            7434,
            1550,
            1900,
            12502,
            9274,
            2802,
            189,
            522,
            3312,
            4469,
            3301,
            4558,
            2193,
            5681,
            204,
            1466,
            1883,
            4476,
            7998,
            1179,
            4160,
            12075,
            312,
            1073,
            1751,
            1950,
            312,
            4224,
            264,
        ],
        dtype=float,
    ),
}


RIGHT_DIVISORS = {
    500: np.array(
        [
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
            501,
        ]
    ),
    3500: np.array(
        [
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
            3501,
        ]
    ),
    7000: np.array(
        [
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
            7001,
        ]
    ),
}


@pytest.fixture
def GeneDatum_dict_test_obj():
    gene_file = "tests/input_data/Test_Genes_NormMatrix.tsv"
    gene_pandas = pd.read_csv(
        gene_file,
        header="infer",
        sep="\t",
        dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
        index_col="Gene_Name",
    )
    sample_genome = GeneData(gene_pandas, "Mock_Camarosa")

    gene_storage_dict = {}
    for name in sample_genome.names:
        gene_storage_dict[name] = sample_genome.get_gene(name)
    return gene_storage_dict


@pytest.mark.parametrize("true_divisors", [LEFT_DIVISORS])
def test_divisor_left_normal(GeneDatum_dict_test_obj, true_divisors):
    """Does the divisor_left function produce the correct values?"""

    for window, true_val_array in true_divisors.items():
        i = 0
        for g_name, g_datum in GeneDatum_dict_test_obj.items():
            expected = g_datum.divisor_left(window)
            assert np.array_equal(true_val_array[i], expected)
            i += 1


@pytest.mark.parametrize("true_divisors", [INTRA_DIVISORS])
def test_divisor_intra_normal(GeneDatum_dict_test_obj, true_divisors):
    """Does the divisor_intra function produce the correct values?"""
    for window, true_val_array in true_divisors.items():
        i = 0
        for g_name, g_datum in GeneDatum_dict_test_obj.items():
            expected = g_datum.divisor_intra(window)
            assert np.array_equal(true_val_array[i], expected)
            i += 1


@pytest.mark.parametrize("true_divisors", [RIGHT_DIVISORS])
def test_divisor_right_normal(GeneDatum_dict_test_obj, true_divisors):
    """Does the divisor_right function produce the correct values?"""
    for window, true_val_array in true_divisors.items():
        i = 0
        for g_name, g_datum in GeneDatum_dict_test_obj.items():
            expected = g_datum.divisor_right(window)
            assert np.array_equal(true_val_array[i], expected)
            i += 1


@pytest.mark.parametrize("true_divisors", [INTRA_INVALID_WINDOW_DIVISORS])
def test_divisor_intra_invalid_windows(GeneDatum_dict_test_obj, true_divisors):
    """Does the divisor_intra function raise the correct error when the window
    is not None?
    """
    with pytest.raises(ValueError):
        for window, true_val_array in true_divisors.items():
            i = 0
            for g_name, g_datum in GeneDatum_dict_test_obj.items():
                expected = g_datum.divisor_intra(window)
                assert np.array_equal(true_val_array[i], expected)
                i += 1


if __name__ == "__main__":
    pytest.main(["-s", __file__])  # for convenience
