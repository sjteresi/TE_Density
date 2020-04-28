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


ground_truth_te = \
"""
Fvb1-1  RepeatMasker    DNA/Mutator     3350    4133    5.0     +       4992    Fxa_V1_4_26930
Fvb1-1  RepeatMasker    LTR/Gypsy       4209    4873    11.7    -       2969    Fvb4-4:6485593..6486917_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       4871    5232    16.0    -       1185    Fvb4-4:6485593..6486917_LTR
Fvb1-1  RepeatMasker    LTR/unknown     15404   16346   19.1    -       3262    Fvb1-1:9463282..9465926_LTR
Fvb1-1  RepeatMasker    LTR/unknown     16347   20224   3.8     -       30004   Fvb3-2:23517681..23521581_INT-int
Fvb1-1  RepeatMasker    LTR/unknown     20225   22848   15.4    -       12678   Fvb1-1:9463282..9465926_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       23450   24094   20.3    -       2636    Fvb1-1:5519794..5522056_LTR
Fvb1-1  RepeatMasker    LTR/unknown     24133   24415   19.0    -       1082    Fvb1-3:25198782..25200749_LTR
Fvb1-1  RepeatMasker    DNA/CMC-EnSpm   24412   24910   8.5     +       2257    family-72
Fvb1-1  RepeatMasker    DNA/CMC-EnSpm   25024   25752   18.8    +       1671    family-72
Fvb1-1  RepeatMasker    DNA/CMC-EnSpm   25751   25840   21.4    +       351     family-72
Fvb1-1  RepeatMasker    LTR/unknown     25842   27393   5.9     +       10968   Fvb6-4:17811941..17814568_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       27398   27999   19.6    +       1680    Fvb4-2:14377692..14383914_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       28097   28708   21.1    +       1651    Fvb4-2:14377692..14383914_LTR
Fvb1-1  RepeatMasker    LTR/unknown     28709   29048   4.4     -       2681    Fvb6-3:14035592..14036822_LTR
Fvb1-1  RepeatMasker    LTR/unknown     29051   30034   5.9     -       7242    Fvb7-3:21810512..21811992_LTR
Fvb1-1  RepeatMasker    LTR/unknown     30038   30447   20.0    +       1455    Fvb1-2:27872851..27875385_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       30440   31593   15.4    -       5517    Fvb1-1:1901529..1903797_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       30986   31928   12.8    +       3981    Fvb1-2:23736536..23738711_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       32011   33581   7.7     +       11550   Fvb6-4:12909946..12918059_INT-int
Fvb1-1  RepeatMasker    LTR/unknown     33582   33661   11.2    -       524     Fvb4-1:5166560..5169132_LTR
Fvb1-1  RepeatMasker    LTR/unknown     33646   33720   17.3    +       409     Fvb1-3:21882267..21885047_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       33721   34100   11.8    +       2503    Fvb4-2:19625958..19630035_INT-int
Fvb1-1  RepeatMasker    LTR/Gypsy       34101   36264   20.8    +       7112    Fvb3-1:4132469..4134698_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       36265   38219   20.9    +       6733    Fvb5-1:12085737..12087980_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       38514   40380   19.2    +       6488    Fvb6-2:9966870..9969479_LTR
Fvb1-1  RepeatMasker    LTR/unknown     40383   41342   5.3     -       7215    Fvb1-2:20827613..20830221_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       41343   41699   12.5    -       7196    Fvb1-1:1901529..1903797_LTR
Fvb1-1  RepeatMasker    LTR/Gypsy       41469   42958   8.3     +       9764    Fvb1-4:19462751..19465717_LTR
Fvb1-1  RepeatMasker    LTR/unknown     42959   45118   8.4     +       14420   Fvb6-4:10839685..10850190_INT-int
"""

@pytest.fixture
def TransposonData_test_obj():
    ground_truth_io = StringIO(ground_truth_te)
    transposon_input_dataframe = import_transposons(ground_truth_io)
    sample_genome = TransposonData(transposon_input_dataframe, 'Mock_Camarosa')
    return sample_tes

# TODO add this to __init__ or to the class along w/ the class method as discussed
ground_truth_genes = \
"""
Fvb1-1	maker	gene	41	2396	.	+	.	ID=maker-Fvb1-1-snap-gene-0.15;Name=maker-Fvb1-1-snap-gene-0.15
Fvb1-1	maker	gene	5556	7978	.	-	.	ID=maker-Fvb1-1-augustus-gene-0.13;Name=maker-Fvb1-1-augustus-gene-0.13
Fvb1-1	maker	gene	8487	8797	.	-	.	ID=maker-Fvb1-1-snap-gene-0.18;Name=maker-Fvb1-1-snap-gene-0.18
Fvb1-1	maker	gene	9361	9658	.	+	.	ID=snap_masked-Fvb1-1-processed-gene-0.6;Name=snap_masked-Fvb1-1-processed-gene-0.6
Fvb1-1	maker	gene	11127	11411	.	-	.	ID=augustus_masked-Fvb1-1-processed-gene-0.4;Name=augustus_masked-Fvb1-1-processed-gene-0.4
Fvb1-1	maker	gene	84598	86703	.	+	.	ID=maker-Fvb1-1-snap-gene-0.16;Name=maker-Fvb1-1-snap-gene-0.16
Fvb1-1	maker	gene	314397	317655	.	-	.	ID=maker-Fvb1-1-augustus-gene-3.19;Name=maker-Fvb1-1-augustus-gene-3.19
Fvb1-1	maker	gene	315831	317608	.	+	.	ID=maker-Fvb1-1-snap-gene-3.20;Name=maker-Fvb1-1-snap-gene-3.20
Fvb1-1	maker	gene	319026	320584	.	+	.	ID=augustus_masked-Fvb1-1-processed-gene-3.0;Name=augustus_masked-Fvb1-1-processed-gene-3.0
Fvb1-1	maker	gene	356220	357714	.	+	.	ID=maker-Fvb1-1-augustus-gene-3.17;Name=maker-Fvb1-1-augustus-gene-3.17
"""


@pytest.fixture
def GeneData_test_obj():
    ground_truth_io = StringIO(ground_truth_genes)
    genes_input_dataframe = import_genes(ground_truth_io)
    sample_genome = GeneData(genes_input_dataframe, 'Mock_Camarosa')
    return sample_genes


MAX_RAM = 1024*4
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
    MergeData._O_INTRA
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
    """Yield a MergeData instance.
    The instance is yielded b/c it has a file resource.
    The resource is managed by the yield.
    """

    yield MergeData.from_param(te_data, gene_data.names, WINDOWS, temp_dir)


# NOTE Scott edit here
# Use the 'real' data sets imported at the top of the file
# Edit the parameters to point to the new data

# Then use merge_sink_real as an input to another functions to check the math
@pytest.yield_fixture
def merge_sink_real(te_data, gene_data, temp_dir):
    """Yield a MergeData instance."""
    # TODO add fixture for the inputs (or just put it here) to the real data
    # basically get the real data into this instead of the dumb data

    return MergeData.from_param(te_data, gene_data.names, WINDOWS, temp_dir)

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
        overlap_source.windows
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
        overlap_source.windows
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
    pytest.main(['-s', __file__])  # for convenience
