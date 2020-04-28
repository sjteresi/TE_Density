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


ground_truth = \
"""
Fvb1-1	maker	gene	41	2396	.	+	.	ID=maker-Fvb1-1-snap-gene-0.15;Name=maker-Fvb1-1-snap-gene-0.15
Fvb1-1	maker	gene	5556	7978	.	-	.	ID=maker-Fvb1-1-augustus-gene-0.13;Name=maker-Fvb1-1-augustus-gene-0.13
Fvb1-1	maker	gene	8487	8797	.	-	.	ID=maker-Fvb1-1-snap-gene-0.18;Name=maker-Fvb1-1-snap-gene-0.18
Fvb1-1	maker	gene	9361	9658	.	+	.	ID=snap_masked-Fvb1-1-processed-gene-0.6;Name=snap_masked-Fvb1-1-processed-gene-0.6
Fvb1-1	maker	gene	11127	11411	.	-	.	ID=augustus_masked-Fvb1-1-processed-gene-0.4;Name=augustus_masked-Fvb1-1-processed-gene-0.4
Fvb1-1	maker	gene	84598	86703	.	+	.	ID=maker-Fvb1-1-snap-gene-0.16;Name=maker-Fvb1-1-snap-gene-0.16
Fvb1-1	maker	gene	314397	317655	.	-	.	ID=maker-Fvb1-1-augustus-gene-3.19;Name=maker-Fvb1-1-augustus-gene-3.19
Fvb1-1	maker	gene	315831	317608	.	+	.	ID=maker-Fvb1-1-snap-gene-3.20;Name=maker-Fvb1-1-snap-gene-3.20
Fvb1-1	maker	gene	319026	320584	.	+	.	ID=augustus_masked-Fvb1-1-processed-gene-3.0;Name=augustus_masked-Fvb1-1-processed-gene-3.0
Fvb1-1	maker	gene	356220	357714	.	+	.	ID=maker-Fvb1-1-augustus-gene-3.17;Name=maker-Fvb1-1-augustus-gene-3.17
"""


@pytest.fixture
def GeneData_test_obj():
    ground_truth_io = StringIO(ground_truth)
    genes_input_dataframe = import_genes(ground_truth_io)
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
