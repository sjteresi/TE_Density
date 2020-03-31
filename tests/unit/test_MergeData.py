#!/usr/bin/env python3

"""
Unit test MergeData.
"""

__author__ = "Michael Teresi"

import logging
import os
import pytest
import tempfile

import coloredlogs
import numpy as np
import pandas as pd

from transposon.merge import MergeData
from transposon.merge import _MergeConfigSink, _MergeConfigSource
from transposon.transposon_data import TransposonData
from transposon.gene_data import GeneData
from transposon.overlap import OverlapData

MAX_RAM = 1024*4
N_TRANSPOSONS = 4
WINDOWS = [10, 20]
CHROME_ID = "mocked_chromosome_id"
LOGGER = logging.getLogger(__name__)
coloredlogs.install(level=logging.DEBUG)
SETS = [
    MergeData._S_LEFT,
    MergeData._S_RIGHT,
    MergeData._S_INTRA,
    MergeData._O_LEFT,
    MergeData._O_RIGHT,
    MergeData._O_INTRA
]


@pytest.fixture
def gene_data():
    """Default GeneData instance."""

    return GeneData.mock()


@pytest.fixture
def te_data():
    """Default TransposonData instance."""

    return TransposonData.mock(chromosome=CHROME_ID)


@pytest.yield_fixture
def temp_dir():
    """Temporary directory."""

    with tempfile.TemporaryDirectory() as dir:
        yield dir


@pytest.yield_fixture
def temp_file(temp_dir):
    """Temporary h5 file."""

    yield os.path.join(temp_dir, 'temp_merge.h5')


@pytest.yield_fixture
def config_sink(te_data, gene_data, temp_file):

    sink = _MergeConfigSink(
        transposons = te_data,
        gene_names = gene_data.names,
        windows = WINDOWS,
        filepath = temp_file,
        ram_bytes = MAX_RAM
    )
    yield sink


@pytest.yield_fixture
def merge_sink(te_data, gene_data, temp_dir):

    yield MergeData.from_param(te_data, gene_data.names, WINDOWS, temp_dir)


@pytest.yield_fixture
def active_merge_sink(merge_sink):

    with merge_sink as active:
        yield active

def test_init_sink(config_sink):
    """Does the initializer raise for a valid sink?"""

    MergeData(config_sink)


def test_from_param(merge_sink):
    """Does the factory from_param raise if valid?"""

    assert isinstance(merge_sink, MergeData)


def test_context_mgr_sink(merge_sink):
    """Does the context manager return the merge data?"""

    with merge_sink as active_sink:
        assert isinstance(merge_sink, MergeData)


@pytest.mark.parametrize("set_name", SETS)
def test_create_set_(active_merge_sink, set_name):
    """Does create_set add the right data sets?"""

    active_merge_sink._h5_file[set_name]


def test_create_set_superfamily(active_merge_sink):
    """Does create_set add the 'public' superfamily densitites?"""

    assert active_merge_sink.superfamily is not None


def test_create_set_super(active_merge_sink):
    """Does create_set add the 'public' order densitites?"""

    assert active_merge_sink.order is not None


# TODO check open_new_file, exists, win /gene names / chr ID /
# TODO test check compression, ok
# TODO test check compression, bad
# TODO test init w/ config source
# TODO test context manager source





if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
