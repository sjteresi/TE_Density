#!/usr/bin/env python3

"""
Unit test ReviseAnno.
"""

__author__ = "Scott Teresi"

import logging
import os
import pytest
import coloredlogs
import numpy as np
import pandas as pd
import logging

from transposon.transposon_data import TransposonData
from transposon.revise_annotation import Revise_Anno

# -------------------------------------------------------------
# TEST VALUES
TRUE_SingleC_SingleE_SUPERFAM = [
    784,
    1021,
    912,
    1694,
    634,
    212,
    150,
    94,
    1432,
    947,
    623,
    177,
    385,
    229,
    2382,
    131,
    189,
    1170,
    501,
    351,
]

TRUE_SingleC_MultiE_SUPERFAM = [
    1412,
    560,
    1219,
    2209,
    212,
    150,
    94,
    2452,
    623,
    177,
    385,
    2611,
    131,
    189,
    1170,
    501,
    351,
]


TRUE_SingleC_ConcOverlap_SUPERFAM = [608, 792, 201]


TRUE_SingleC_SingleE_ORDER = [
    784,
    1021,
    912,
    1694,
    1170,
    501,
    351,
]

TRUE_SingleC_MultiE_ORDER = [1880, 1219, 2209, 301, 212, 150, 94, 2452, 2611]
TRUE_SingleC_ConcOverlap_ORDER = [608, 792, 201]


TRUE_SingleC_SingleE_NAMELESS = [
    784,
    1021,
    912,
    1694,
    1170,
    551,
]
TRUE_SingleC_MultiE_NAMELESS = [1880, 1219, 2209, 212, 150, 2545, 2611]
TRUE_SingleC_ConcOverlap_NAMELESS = [608, 792, 201]


# -------------------------------------------------------------


@pytest.fixture
def logger_obj():
    """
    Dummy logging object for Revise_Anno constructor
    """
    logger = logging.getLogger(__name__)
    return logger


@pytest.fixture
def h5_cache_loc():
    """
    Location for outputting h5 revised annotation files.
    """
    path_main = os.path.abspath(__file__)
    h5_cache_loc = "tests/test_h5_cache_loc/"
    if not os.path.isdir(os.path.abspath(os.path.join(path_main, h5_cache_loc))):
        os.makedirs(h5_cache_loc, exist_ok=True)
    return h5_cache_loc


@pytest.fixture
def revised_te_annotation_loc():
    """
    Location for outputing revised annotation files for visual inspection.
    """
    path_main = os.path.abspath(__file__)
    data_output = "tests/output_data/"
    if not os.path.isdir(os.path.abspath(os.path.join(path_main, data_output))):
        os.makedirs(data_output, exist_ok=True)
    return data_output


@pytest.fixture
def superfam_name():
    """
    Dummy name used for Revise_Anno constructor
    """
    return "test_superfam_set"


@pytest.fixture
def order_name():
    """
    Dummy name used for Revise_Anno constructor
    """
    return "test_order_set"


@pytest.fixture
def nameless_name():
    """
    Dummy name used for Revise_Anno constructor
    """
    return "test_nameless_set"


# Input Data
path_main = os.path.abspath(__file__)
data_input = "tests/input_data/"
if not os.path.isdir(os.path.abspath(os.path.join(path_main, data_input))):
    os.makedirs(data_input, exist_ok=True)
# Supers
SingleC_SingleElongate_Super = (
    "tests/input_data/Test_SingleC_SingleElongate_Superfam_Revision.tsv"
)
SingleC_MultiElongate_Super = (
    "tests/input_data/Test_SingleC_MultiElongate_Superfam_Revision.tsv"
)
SingleC_Conc_Super = (
    "tests/input_data/Test_SingleC_ConcurrentOverlap_Superfam_Revision.tsv"
)
# Orders
SingleC_SingleElongate_Order = (
    "tests/input_data/Test_SingleC_SingleElongate_Order_Revision.tsv"
)
SingleC_MultiElongate_Order = (
    "tests/input_data/Test_SingleC_MultiElongate_Order_Revision.tsv"
)
SingleC_Conc_Order = (
    "tests/input_data/Test_SingleC_ConcurrentOverlap_Order_Revision.tsv"
)
# Nameless
SingleC_SingleElongate_Nameless = (
    "tests/input_data/Test_SingleC_SingleElongate_Nameless_Revision.tsv"
)
SingleC_MultiElongate_Nameless = (
    "tests/input_data/Test_SingleC_MultiElongate_Nameless_Revision.tsv"
)
SingleC_Conc_Nameless = (
    "tests/input_data/Test_SingleC_ConcurrentOverlap_Nameless_Revision.tsv"
)


@pytest.fixture
def TEData_TestObj(request):
    """
    Give a TransposonData object based on an input dataframe
    """
    # NOTE since we are parametrizing this fixture, the argument has to be
    # named request.... I do not know why, can't figure out
    transposon_input_dataframe = pd.read_csv(
        request.param,
        header="infer",
        dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
        sep="\t",
    )
    te_data = TransposonData(transposon_input_dataframe, "Mock_Camarosa")
    return te_data


@pytest.mark.parametrize(
    "TEData_TestObj, true_values, output_filenames",
    [
        (
            SingleC_SingleElongate_Super,
            TRUE_SingleC_SingleE_SUPERFAM,
            "SingleC_SingleE_Super.tsv",
        ),
        (
            SingleC_MultiElongate_Super,
            TRUE_SingleC_MultiE_SUPERFAM,
            "SingleC_MultiE_Super.tsv",
        ),
        (
            SingleC_Conc_Super,
            TRUE_SingleC_ConcOverlap_SUPERFAM,
            "SingleC_Conc_Super.tsv",
        ),
    ],
    indirect=["TEData_TestObj"],
)
def test_superfam(
    TEData_TestObj,
    true_values,
    output_filenames,
    h5_cache_loc,
    superfam_name,
    revised_te_annotation_loc,
):
    """Create superfamily revisions"""
    revise_anno_obj = Revise_Anno(
        TEData_TestObj.data_frame, "test.tsv", h5_cache_loc, superfam_name
    )
    revise_anno_obj.create_superfam()
    observed = revise_anno_obj.updated_te_annotation.Length.to_numpy(copy=False)

    # NB this is just for manual inspection
    revise_anno_obj._write(
        revise_anno_obj.updated_te_annotation,
        os.path.join(revised_te_annotation_loc, output_filenames),
    )
    assert np.array_equal(observed, true_values)


@pytest.mark.parametrize(
    "TEData_TestObj, true_values, output_filenames",
    [
        (
            SingleC_SingleElongate_Order,
            TRUE_SingleC_SingleE_ORDER,
            "SingleC_SingleE_Order.tsv",
        ),
        (
            SingleC_MultiElongate_Order,
            TRUE_SingleC_MultiE_ORDER,
            "SingleC_MultiE_Order.tsv",
        ),
        (
            SingleC_Conc_Order,
            TRUE_SingleC_ConcOverlap_ORDER,
            "SingleC_Conc_Order.tsv",
        ),
    ],
    indirect=["TEData_TestObj"],
)
def test_order(
    TEData_TestObj,
    true_values,
    output_filenames,
    h5_cache_loc,
    order_name,
    revised_te_annotation_loc,
):
    """Create order revisions"""
    revise_anno_obj = Revise_Anno(
        TEData_TestObj.data_frame, "test.tsv", h5_cache_loc, order_name
    )
    revise_anno_obj.create_order()
    observed = revise_anno_obj.updated_te_annotation.Length.to_numpy(copy=False)

    # NB this is just for manual inspection
    revise_anno_obj._write(
        revise_anno_obj.updated_te_annotation,
        os.path.join(revised_te_annotation_loc, output_filenames),
    )
    assert np.array_equal(observed, true_values)


@pytest.mark.parametrize(
    "TEData_TestObj, true_values, output_filenames",
    [
        (
            SingleC_SingleElongate_Nameless,
            TRUE_SingleC_SingleE_NAMELESS,
            "SingleC_SingleE_Nameless.tsv",
        ),
        (
            SingleC_MultiElongate_Nameless,
            TRUE_SingleC_MultiE_NAMELESS,
            "SingleC_MultiE_Nameless.tsv",
        ),
        (
            SingleC_Conc_Nameless,
            TRUE_SingleC_ConcOverlap_NAMELESS,
            "SingleC_Conc_Nameless.tsv",
        ),
    ],
    indirect=["TEData_TestObj"],
)
def test_nameless(
    TEData_TestObj,
    true_values,
    output_filenames,
    h5_cache_loc,
    nameless_name,
    revised_te_annotation_loc,
):
    """Create nameless revisions"""
    revise_anno_obj = Revise_Anno(
        TEData_TestObj.data_frame, "test.tsv", h5_cache_loc, nameless_name
    )
    revise_anno_obj.create_nameless()
    observed = revise_anno_obj.updated_te_annotation.Length.to_numpy(copy=False)

    # NB this is just for manual inspection
    revise_anno_obj._write(
        revise_anno_obj.updated_te_annotation,
        os.path.join(revised_te_annotation_loc, output_filenames),
    )
    assert np.array_equal(observed, true_values)


if __name__ == "__main__":
    pytest.main(["-s", __file__])  # for convenience
