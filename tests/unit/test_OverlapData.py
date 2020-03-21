#!/usr/bin/env python3

"""
Unit test OverlapData.

OverlapData can be used in two modes: outut or input file, SEE `from_param`, `from_file`.
Tests are grouped accordingly.
"""

__author__ = "Michael Teresi"

import logging
import os
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

def default_data_out(temp_dir):
    """Return default output OverlapData instance."""

    return OverlapData.from_param(
        GeneData.mock(), N_TRANSPOSONS, WINDOWS, temp_dir, logger=LOGGER
    )

@pytest.yield_fixture
def overlap_output(gene_data, temp_dir):
    """Default OverlapData instance for writing."""

    with tempfile.TemporaryDirectory() as temp_dir:
        with default_data_out(temp_dir) as active_output:
            yield active_output

# TODO create mocked config source
# TODO test check ram
# TODO test init throw, source

def test_from_param_raise(gene_data, temp_dir):
    """Does the from param factory raise?"""

    OverlapData.from_param(gene_data, N_TRANSPOSONS, WINDOWS, temp_dir, logger=LOGGER)

def test_from_param_raise_enter_exit(overlap_output):
    """Does the context manager raise for an output file?"""

    pass

def test_open_dispatch_bad(temp_dir):
    """Does the open dispatch raise on an invalid config?"""

    class DummyClass():
        pass
    od = OverlapData(DummyClass())
    with pytest.raises(TypeError) as excinfo:
        od._open_dispatcher()

def test_open_dispatch_sink(temp_dir):
    """Does the open dispatch create a file?"""

    with default_data_out(temp_dir) as io:
        assert os.path.isfile(io.h5_file.filename)

def _test_open_dispatch_source(temp_dir):
    """Does the open dispatch call the source initializer?"""

    raise NotImplementedError()

def test_from_file_raise(temp_dir):
    """Does the from file factory raise?"""

    overlap_data = OverlapData.from_file(temp_dir, LOGGER)

if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
