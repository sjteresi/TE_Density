#!/usr/bin/env python3

"""
Unit test MergeData.
"""

__author__ = "Michael Teresi, Scott Teresi"

import coloredlogs
import numpy as np
import pandas as pd
import logging
import pytest
import tempfile


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
# SCOTT tests with real data below
# NOTE window values used in tests
windows_real = [500, 1000, 1500, 2000]


@pytest.fixture
def genedata_test_obj():
    """Create test object for GeneData information, reads from file"""
    # NOTE this path is relative to the main directory
    gene_file = "tests/input_data/Test_Genes_MergeData.tsv"
    gene_pandas = pd.read_csv(
        gene_file,
        header="infer",
        sep="\t",
        dtype={"Start": "float64", "Stop": "float64", "Length": "float64"},
        index_col="Gene_Name",
    )
    sample_genome = GeneData(gene_pandas, "Mock_Camarosa")
    return sample_genome


@pytest.fixture
def transposondata_test_obj():
    """Create test object for TransposonData information, reads from file"""
    # NOTE this path is relative to the main directory
    te_file = "tests/input_data/Test_TEs_MergeData.tsv"
    te_pandas = pd.read_csv(
        te_file,
        header="infer",
        sep="\t",
        dtype={"Start": "float64", "Stop": "float64", "Length": "float64"},
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
def active_overlap_data_real(genedata_test_obj, transposondata_test_obj, temp_file):
    """Yield an active OverlapData instance from real data"""

    overlap_worker = OverlapWorker(temp_file)
    overlap_worker.calculate(
        genedata_test_obj,
        transposondata_test_obj,
        windows_real,
        genedata_test_obj.names,
    )

    with OverlapData.from_file(temp_file) as active_overlap:
        yield active_overlap


@pytest.fixture
def list_sum_args_real(active_merge_sink_real, active_overlap_data_real):
    sums = active_merge_sink_real._list_density_args(active_overlap_data_real)
    return sums


def test_merge_summed_real(
    active_merge_sink_real, active_overlap_data_real, genedata_test_obj
):
    """
    Can MergeData calculate the sum of overlaps?
    NB, MergeData.sum also does the division step, this will later be
    refactored but important to have a system test for this.

    System test for the final output of MergeData because difficult to unit
    test the many pieces of MergeData at the moment.
    """
    logging.info("test_merge_summed_real...")
    active_merge_sink_real.sum(active_overlap_data_real, genedata_test_obj)

    # print(active_merge_sink_real.filepath)
    # print(active_merge_sink_real.gene_names)
    # print(active_merge_sink_real.windows)
    # print(active_merge_sink_real.order_names)
    # print(active_merge_sink_real.superfamily_names)
    # print(active_merge_sink_real._order_2_idx)

    order_idx = active_merge_sink_real._order_2_idx["Total_TE_Density"]

    # TEST ONE
    gene_name = "dummy_gene_1"
    gene_idx = active_merge_sink_real._gene_2_idx[gene_name]
    lr_slice = MergeData.left_right_slice(order_idx, gene_idx, 3)
    gene_datum = genedata_test_obj.get_gene(gene_name)
    divisor = gene_datum.divisor_right(2000)
    numpy_val_to_compare = active_merge_sink_real.order.right[lr_slice].ravel().item()
    expected = np.divide((4100 - 3350 + 1), divisor)
    assert numpy_val_to_compare == pytest.approx(expected, rel=1e-6)

    # TEST TWO
    gene_name = "dummy_gene_3"
    gene_idx = active_merge_sink_real._gene_2_idx[gene_name]
    lr_slice = MergeData.left_right_slice(order_idx, gene_idx, 3)
    gene_datum = genedata_test_obj.get_gene(gene_name)
    divisor = gene_datum.divisor_right(2000)
    numpy_val_to_compare = active_merge_sink_real.order.right[lr_slice].ravel().item()
    expected = np.divide(
        ((9370 - 8500 + 1) + (9677 - 9556 + 1) + (10500 - 9678 + 1)), divisor
    )
    assert numpy_val_to_compare == pytest.approx(expected, rel=1e-6)

    # TEST THREE
    gene_name = "dummy_gene_2"
    gene_idx = active_merge_sink_real._gene_2_idx[gene_name]
    lr_slice = MergeData.left_right_slice(order_idx, gene_idx, 1)
    gene_datum = genedata_test_obj.get_gene(gene_name)
    divisor = gene_datum.divisor_left(1000)
    numpy_val_to_compare = active_merge_sink_real.order.left[lr_slice].ravel().item()
    # NOTE confusing because window actually goes from (gene_start - 1 -
    # window) because the gene start is part of the gene, not the window, so
    # there is a 1 int offset. Since this TE starts outside of the window, and
    # ends inside, we have to make sure we subtract the TE stop by the correct
    # window start, which in this unique case is 4999, not 5000 which would be
    # (gene start - window value). Here 4999 represents (gene start - window -
    # 1)
    expected = np.divide((5229 - 4999 + 1), divisor)
    assert numpy_val_to_compare == pytest.approx(expected, rel=1e-6)

    # TEST FOUR (INTRA)
    gene_name = "dummy_gene_3"
    gene_idx = active_merge_sink_real._gene_2_idx[gene_name]
    lr_slice = MergeData.intra_slice(order_idx, gene_idx)
    gene_datum = genedata_test_obj.get_gene(gene_name)
    divisor = gene_datum.divisor_intra(None)
    numpy_val_to_compare = active_merge_sink_real.order.intra[lr_slice].ravel().item()
    expected = np.divide((8500 - 8459 + 1), divisor)
    assert numpy_val_to_compare == pytest.approx(expected, rel=1e-6)


def test_slice_left_right_real(active_merge_sink_real):
    """Can we get a slice from on active merge sink of real data?"""

    active_merge_sink_real.left_right_slice(1, 0, 1)


def test_slice_intra_real(active_merge_sink_real):
    """Can we get a slice from on active merge sink of real data?"""

    active_merge_sink_real.intra_slice(1, 0, None)


def test_list_sum_args_no_throw_real(active_merge_sink_real, active_overlap_data_real):
    """Show useful information for each calculation, tests to see if can
    construct the arguments correctly"""

    out = active_merge_sink_real._list_density_args(active_overlap_data_real)
    for o in out:
        # print(o)
        # print()
        assert isinstance(o, _SummationArgs)


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

    # NOTE we don't use active_sink here?
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

    out = active_merge_sink._list_density_args(overlap_source)
    for o in out:
        assert isinstance(o, _SummationArgs)


@pytest.mark.skip(reason="TODO")
def test_sum_no_throw(active_merge_sink, overlap_source):
    # NOTE FAILS
    pass
    # active_merge_sink.sum(overlap_source, None)
    # Not the right place to start


@pytest.mark.skip(reason="TODO")
def test_process_sum():
    # Not the right place to start
    pass


if __name__ == "__main__":
    pytest.main(["-s", __file__])  # for convenience
