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
from io import StringIO

from transposon.merge_data import MergeData
from transposon.merge_data import _MergeConfigSink, _MergeConfigSource, _SummationArgs
from transposon.transposon_data import TransposonData
from transposon.gene_data import GeneData
from transposon.overlap import OverlapData, OverlapWorker


MAX_RAM = 1024 * 4
N_TRANSPOSONS = 4
WINDOWS = [10, 20, 30, 40]
CHROME_ID = "fake_chromosome_id"
LOGGER = logging.getLogger(__name__)
coloredlogs.install(level=logging.DEBUG)
SETS = [
    MergeData._S_LEFT,
    MergeData._S_RIGHT,
    MergeData._S_INTRA,
    MergeData._O_LEFT,
    MergeData._O_RIGHT,
    MergeData._O_INTRA,
]


def _gene_data():

    return GeneData.mock(genome_id=CHROME_ID)


def _te_data():

    return TransposonData.mock(chromosome=CHROME_ID)


@pytest.fixture
def gene_data():
    """Default GeneData instance."""

    return _gene_data()


@pytest.fixture
def te_data():
    """Default TransposonData instance."""

    return _te_data()


@pytest.yield_fixture()
def temp_dir():
    """Temporary directory."""

    with tempfile.TemporaryDirectory() as dir:
        yield dir


@pytest.yield_fixture
def temp_file(temp_dir):
    """Temporary h5 file."""

    yield os.path.join(temp_dir, "temp_merge.h5")


@pytest.yield_fixture
def config_sink(te_data, gene_data, temp_file):

    sink = _MergeConfigSink(
        transposons=te_data,
        gene_names=gene_data.names,
        windows=WINDOWS,
        filepath=temp_file,
        ram_bytes=MAX_RAM,
    )
    yield sink


@pytest.yield_fixture
def merge_sink(te_data, gene_data, temp_dir):
    """Yield a MergeData instance.
    The instance is yielded b/c it has a file resource.
    The resource is managed by the yield.
    """

    yield MergeData.from_param(te_data, gene_data.names, WINDOWS, temp_dir)


# NOTE Scott edit here
# Use the 'real' data sets imported at the top of the file
# Edit the parameters to point to the new data

# -----------------------------------------------------
# SCOTT
@pytest.fixture
def GeneData_test_obj():
    gene_file = "tests/input_data/Test_Genes_MergeData.tsv"
    gene_pandas = pd.read_csv(
        gene_file,
        header="infer",
        sep="\t",
        dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
        index_col="Gene_Name",
    )
    sample_genome = GeneData(gene_pandas, "Mock_Camarosa")
    return sample_genome


@pytest.fixture
def TransposonData_test_obj():
    te_file = "tests/input_data/Test_TEs_MergeData.tsv"
    te_pandas = pd.read_csv(
        te_file,
        header="infer",
        sep="\t",
        dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
    )
    sample_genome = TransposonData(te_pandas)
    return sample_genome


# Then use merge_sink_real as an input to another functions to check the math
@pytest.yield_fixture
def merge_sink_real(TransposonData_test_obj, GeneData_test_obj, temp_dir):
    """Yield a MergeData instance."""
    # TODO add fixture for the inputs (or just put it here) to the real data
    # basically get the real data into this instead of the dumb data
    return MergeData.from_param(
        TransposonData_test_obj, GeneData_test_obj.names, windows_real, temp_dir
    )


TRUE_SUPER_VALS = 1

windows_real = [500, 1000, 1500, 2000]


@pytest.mark.parametrize("true_summed_overlaps", [TRUE_SUPER_VALS])
def test_supers_real(merge_sink_real, true_summed_overlaps):
    print()
    # Superfamily|Order, window, gene
    with merge_sink_real as active_merge_sink_real:
        my_slice = MergeData.left_right_slice(1, 0, 1)
        print(active_merge_sink_real.order_names)
        print(active_merge_sink_real.superfamily_names)
        print(my_slice)
        print()
        # Matrix, Row, Column
        print(active_merge_sink_real.superfamily[my_slice[2]])


# SCOTT
# -----------------------------------------------------
# ACCESSORY


# SCOTT
# -----------------------------------------------------
@pytest.yield_fixture
def active_merge_sink(merge_sink):

    with merge_sink as active:
        yield active


@pytest.fixture(scope="module")
def overlap_source():

    with tempfile.TemporaryDirectory() as temp_dir:
        worker = OverlapWorker(temp_dir)
        gene_data = _gene_data()
        te_data = _te_data()
        overlap_file = worker.calculate(gene_data, te_data, WINDOWS, gene_data.names)
        with OverlapData.from_file(overlap_file) as active_source:
            yield active_source


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


def test_list_sum_input_outputs_super_no_throw(active_merge_sink, overlap_source):
    """Does it produce the right tuple for superfamily inputs?"""

    out = active_merge_sink._list_sum_input_outputs(
        overlap_source,
        active_merge_sink.superfamily,
        active_merge_sink._config.transposons.superfamilies,
        active_merge_sink._config.transposons.superfamily_name_set,
        active_merge_sink._superfam_2_idx,
        overlap_source.windows,
    )
    for o in out:
        assert isinstance(o, _SummationArgs)


def test_list_sum_input_outputs_order_no_throw(active_merge_sink, overlap_source):
    """Does it produce the right tuple for order inputs?"""

    out = active_merge_sink._list_sum_input_outputs(
        overlap_source,
        active_merge_sink.order,
        active_merge_sink._config.transposons.orders,
        active_merge_sink._config.transposons.order_name_set,
        active_merge_sink._order_2_idx,
        overlap_source.windows,
    )
    for o in out:
        assert isinstance(o, _SummationArgs)


def test_list_sum_args_no_throw(active_merge_sink, overlap_source):
    """Does it produce the right tuples?"""

    out = active_merge_sink._list_sum_args(overlap_source)
    for o in out:
        assert isinstance(o, _SummationArgs)


def test_sum_no_throw(active_merge_sink, overlap_source):

    active_merge_sink.sum(overlap_source, None)


def test_process_sum():
    pass


def test_real_init(merge_sink_real):
    """Can we acquire the real data resource w/o it raising?"""

    with merge_sink_real as active:
        pass


if __name__ == "__main__":
    pytest.main(["-s", __file__])  # for convenience
