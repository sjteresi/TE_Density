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
from collections import namedtuple

from transposon.merge_data import MergeData
from transposon.merge_data import _MergeConfigSink, _MergeConfigSource, _SummationArgs
from transposon.transposon_data import TransposonData
from transposon.gene_data import GeneData
from transposon.overlap import OverlapData, OverlapWorker

# pytestmark = pytest.mark.skip


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

    file = tempfile.NamedTemporaryFile(dir=temp_dir, suffix=".h5")
    with file as temp:
        yield temp.name


@pytest.yield_fixture
def config_sink(te_data, gene_data, temp_file):
    sink = _MergeConfigSink(
        transposons=te_data,
        gene_data=gene_data,
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

    yield MergeData.from_param(te_data, gene_data, WINDOWS, temp_dir)


# -----------------------------------------------------
# SCOTT


windows_real = [500, 1000, 1500, 2000]  # NOTE used for the
# MergeData.from_param parts


@pytest.fixture
def genedata_test_obj():
    """Create test object for GeneData information, reads from file"""
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
def transposondata_test_obj():
    """Create test object for TransposonData information, reads from file"""
    te_file = "tests/input_data/Test_TEs_MergeData.tsv"
    te_pandas = pd.read_csv(
        te_file,
        header="infer",
        sep="\t",
        dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
    )
    sample_genome = TransposonData(te_pandas, "Mock_Camarosa")
    return sample_genome


@pytest.yield_fixture
def active_merge_sink_real(transposondata_test_obj, genedata_test_obj, temp_dir):
    """Yield an active MergeData instance from real data"""
    my_merge = MergeData.from_param(
        transposondata_test_obj, genedata_test_obj, windows_real, temp_dir
    )
    with my_merge as active_merge:
        yield active_merge


@pytest.yield_fixture
def active_real_overlap_data(genedata_test_obj, transposondata_test_obj, temp_file):
    """Yield an active OverlapData instance from real data"""
    my_overlap_data = OverlapData.from_param(
        genedata_test_obj,
        transposondata_test_obj.number_elements,
        windows_real,
        temp_file,
    )
    with my_overlap_data as active_overlap:
        yield active_overlap


@pytest.mark.skip(reason="TODO")
def test_scott_process_sum(
    active_merge_sink_real,
    active_real_overlap_data,
    genedata_test_obj,
    list_sum_args_real,
):

    # This was my attempt at reproducing the sum issue.
    # But I wasn't able to get the stuff preceeding it to work.

    for sum_args in list_sum_args_real:

        for gene_name in active_real_overlap_data.gene_names:
            gene_datum = genedata_test_obj.get_gene(gene_name)
            g_idx = active_merge_sink_real._gene_2_idx[gene_name]

            for te_ in zip(*sum_args.te_idx_name):  # for every superfam | order
                te_idx, te_name = te_  # the supefamily / order index and string

                for window in sum_args.windows:  # for every window value
                    w_idx = self._window_2_idx.get(window, None)

                    slice_out = sum_args.slice_out(
                        window_idx=w_idx, gene_idx=g_idx, group_idx=te_idx
                    )
                    slice_in = sum_args.slice_in(w_idx, g_idx)
                    # find which genes match the superfam | order, and sum those
                    superfam_or_order_match = sum_args.where(te_name)
                    slice_in_only_te = (superfam_or_order_match, *slice_in[1:])

                    # BUG below this is what fails
                    # slice_sum = sum_args.input[slice_in_only_te]


@pytest.fixture
def list_sum_args_real(active_merge_sink_real, active_real_overlap_data):
    sums = active_merge_sink_real._list_sum_args(active_real_overlap_data)
    return sums


@pytest.mark.skip(reason="TODO")
def test_merge_real_summed(
    active_merge_sink_real, active_real_overlap_data, genedata_test_obj
):
    """Can MergeData calculate the sum of overlaps"""
    merge_real_summed = active_merge_sink_real.sum(
        active_real_overlap_data, genedata_test_obj
    )


def test_slice_left_right(active_merge_sink_real):
    """Can we get a slice from on active merge sink of real data?"""

    active_merge_sink_real.left_right_slice(1, 0, 1)


def test_slice_intra(active_merge_sink_real):
    """Can we get a slice from on active merge sink of real data?"""

    active_merge_sink_real.intra_slice(1, 0, None)


# -----------------------------------------------------
@pytest.yield_fixture
def active_merge_sink(merge_sink):

    with merge_sink as active:
        yield active


@pytest.fixture(scope="module")
def overlap_source():

    with tempfile.NamedTemporaryFile(suffix=".h5") as temp_file:
        worker = OverlapWorker(temp_file.name)
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


@pytest.mark.skip(reason="TODO")
def test_sum_no_throw(active_merge_sink, overlap_source):
    # NOTE FAILS
    pass
    # active_merge_sink.sum(overlap_source, None)


@pytest.mark.skip(reason="TODO")
def test_process_sum():
    pass


if __name__ == "__main__":
    pytest.main(["-s", __file__])  # for convenience
