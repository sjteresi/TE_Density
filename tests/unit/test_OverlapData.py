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

@pytest.fixture
def temp_dir():
    """Temporary directory."""

    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


@pytest.fixture
def temp_file():
    """Temporary directory."""

    with tempfile.NamedTemporaryFile(suffix="."+OverlapData.EXT) as temp_file:
        yield temp_file.name

@pytest.fixture
def default_data_out(temp_file):
    """Return default output OverlapData instance."""

    return OverlapData.from_param(
        GeneData.mock(), N_TRANSPOSONS, WINDOWS, temp_file, logger=LOGGER
    )

@pytest.fixture
def active_output(default_data_out):
    """Default OverlapData instance for writing."""

    with default_data_out as active_output:
        yield active_output

@pytest.fixture
def active_input(default_data_out):
    """Default OverlapData instance for reading."""

    filepath = None
    with default_data_out as io:
        filepath = io.filepath
    with OverlapData.from_file(filepath) as io:
        yield io

@pytest.fixture
def serialized_deserialized(default_data_out):
    """Yield an overlap written to disk and one read from the first."""

    filepath = None
    with default_data_out as output:
        # MAGIC NUMBER dummy data
        output.left[:] = np.ones(output.left.shape) * 2
        output.right[:] = np.ones(output.right.shape) * 3
        output.intra[:] = np.ones(output.intra.shape) * 4
        output._h5_file.flush()
        filepath = output.filepath
        with OverlapData.from_file(filepath) as input:
            yield (input, output)

def test_from_param_raise(gene_data, temp_file):
    """Does the from param factory raise?"""

    OverlapData.from_param(gene_data, N_TRANSPOSONS, WINDOWS, temp_file, logger=LOGGER)

def test_from_param_raise_enter_exit(active_output):
    """Does the context manager raise for an output file?"""

    pass

def test_open_dispatch_bad(temp_dir):
    """Does the open dispatch raise on an invalid config?"""

    class DummyClass():
        pass
    od = OverlapData(DummyClass())
    with pytest.raises(TypeError) as excinfo:
        od._open_dispatcher()

def test_open_dispatch_sink(active_output):
    """Does the open dispatch create a file?"""

    assert os.path.isfile(active_output._h5_file.filename)

def _test_open_dispatch_source(temp_dir):
    """Does the open dispatch call the source initializer?"""

    raise NotImplementedError()

def test_from_file_raise_valid(default_data_out):
    """Does the from file factory raise for valid data?"""

    filepath = None
    with default_data_out as io:
        filepath = io.filepath
    overlap_data = OverlapData.from_file(filepath, LOGGER)

def test_from_file_raise(temp_dir):
    """Does the from file factory raise for invalid data?"""

    with pytest.raises(ValueError) as excinfo:
        overlap_data = OverlapData.from_file('not a file')

def test_from_file_raise_enter_exit(default_data_out):
    """Does the context manager raise?"""

    filepath = None
    with default_data_out as io:
        filepath = io.filepath
    with OverlapData.from_file(filepath) as io:
        assert io.filepath == filepath

def test_left_right_shape(active_output):
    """Do the shapes match for left / right overlap?"""

    assert np.all(active_output.left.shape == active_output.right.shape)

def test_io_windows(serialized_deserialized):
    """Can it serialize / deserialize the windows?"""

    input, output = serialized_deserialized
    assert input.windows == output.windows

def test_io_genes(serialized_deserialized):
    """Can it serialize / deserialize the gene names?"""

    input, output = serialized_deserialized
    assert input.gene_names == output.gene_names

def test_io_chromosome_id(serialized_deserialized):
    """Can it serialize / deserialize the chromosome id?"""

    input, output = serialized_deserialized
    assert input.chromosome_id == output.chromosome_id

def test_io_chromosome_id(serialized_deserialized):
    """Can it serialize / deserialize the chromosome id?"""

    input, output = serialized_deserialized
    assert input.genome_id == output.genome_id

def test_io_left(serialized_deserialized):
    """Can it serialize / deserialize the left overlap"""

    input, output = serialized_deserialized
    assert np.all(input.left == output.left)

def test_io_intra(serialized_deserialized):
    """Can it serialize / deserialize the intra overlap?"""

    input, output = serialized_deserialized
    assert np.all(input.intra == output.intra)

def test_io_right(serialized_deserialized):
    """Can it serialize / deserialize the right overlap?"""

    input, output = serialized_deserialized
    assert np.all(input.right == output.right)

# TODO test slicing

if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
