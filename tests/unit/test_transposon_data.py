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

    raise NotImplementedError()

def test_subset_id_missing():
    """Does the chromosome identifier raise if the IDs are missing?"""

    raise NotImplementedError()

def test_subset_id_not_unique():
    """Does the property raise if the chromosome IDs aren't unique?"""

    raise NotImplementedError()



if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
