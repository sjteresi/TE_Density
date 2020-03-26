#!/usr/bin/env python3

"""
Unit test transposon_data.py
"""

__author__ = "Michael Teresi"

import pytest

import numpy as np
import pandas as pd

from transposon.gene_data import GeneData
from transposon.gene_datum import GeneDatum
from transposon.transposon_data import TransposonData

@pytest.fixture
def te_data():
    """Default TransposonData instance."""

    return TransposonData.mock()

def test_init(te_data):
    """Does the initializer fail?"""

    pass

def test_subset_id_unique():
    """Does the chromosome identifier work if the IDs are the same?"""

    transposons = TransposonData.mock(np.array([[0, 9], [10, 19]]))
    same_chromosome_name = "this_name_is_consistent"
    transposons.chromosomes[0] = same_chromosome_name
    transposons.chromosomes[1] = same_chromosome_name
    assert transposons.chromosome_unique_id == same_chromosome_name

def test_subset_id_missing():
    """Does the chromosome identifier raise if the IDs are missing?"""

    transposons = TransposonData.mock(np.array([]))
    with pytest.raises(RuntimeError) as excinfo:
        transposons.chromosome_unique_id

def test_subset_id_not_unique():
    """Does the property raise if the chromosome IDs aren't unique?"""

    transposons = TransposonData.mock(np.array([[3, 9], [12, 25]]))
    transposons.chromosomes[0] = "first_not_unique_chromosome_name"
    transposons.chromosomes[1] = 'Different_chromosome_name_from_default'
    with pytest.raises(RuntimeError) as excinfo:
        transposons.chromosome_unique_id



if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
