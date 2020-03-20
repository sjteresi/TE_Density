#!/usr/bin/env python3

"""
Unit test OverlapData.
"""

__author__ = "Michael Teresi"

import logging
import pytest
import tempfile

import coloredlogs
import numpy as np
import pandas as pd

from transposon.gene_data import GeneData
from transposon.overlap import OverlapData

N_TRANSPOSONS = 4
WINDOWS = [10, 20]
LOGGER = logging.getLogger(__name__)
coloredlogs.install(level=logging.DEBUG)


@pytest.fixture
def gene_data():
    """Default GeneData instance."""

    return GeneData.mock()

@pytest.yield_fixture
def temp_dir():
    """Temporary directory."""

    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir

@pytest.yield_fixture
def overlap_data(gene_data):
    """Default OverlapData instance."""

    with tempfile.TemporaryDirectory() as temp_dir:
        yield OverlapData(gene_data, 2, WINDOWS, temp_dir, LOGGER)

def test_init_raise(overlap_data):
    """Does the initializer raise?"""

    pass

def test_from_param_raise(gene_data, temp_dir):
    """Does the from param factory raise?"""

    overlap_data = OverlapData.from_param(gene_data, N_TRANSPOSONS, temp_dir, LOGGER)

def test_from_file_raise(temp_dir):
    """Does the from file factory raise?"""

    overlap_data = OverlapData.from_file(temp_dir, LOGGER)

if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
