#!/usr/bin/env python3

"""
Unit test data.py
"""

__author__ = "Michael Teresi"

import pytest

import numpy as np
import pandas as pd

from transposon.data import GeneData, TransposonData

@pytest.fixture
def gene_data():
    """Default GeneData instance."""

    return GeneData.mock()

def test_init(gene_data):
    """Does the initializer fail?"""

    pass

def test_subset_id_unique():
    """Does the chromosome identifier work if the IDs are the same?"""

    genes = GeneData.mock(np.array([[0, 9], [10, 19]]))
    same_chromosome_name = "this_name_is_consistent"
    genes.chromosomes[0] = same_chromosome_name
    genes.chromosomes[1] = same_chromosome_name
    assert genes.chromosome_unique_id == same_chromosome_name

def test_subset_id_missing():
    """Does the chromosome identifier raise if the IDs are missing?"""

    genes = GeneData.mock(np.array([]))
    with pytest.raises(RuntimeError) as excinfo:
        genes.chromosome_unique_id

def test_subset_id_not_unique():
    """Does the property raise if the chromosome IDs aren't unique?"""

    genes = GeneData.mock(np.array([[0, 9], [10, 19]]))
    genes.chromosomes[0] = "first_not_unique_chromosome_name"
    genes.chromosomes[1] = "this_one_is_different"
    with pytest.raises(RuntimeError) as excinfo:
        genes.chromosome_unique_id


if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
