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
    h5_cache_loc = 'tests/test_h5_cache_loc/'
    return h5_cache_loc


@pytest.fixture
def revised_te_annotation_loc():
    """
    Location for outputing revised annotation files for visual inspection.
    """
    data_output = 'tests/output_data/'
    return data_output


@pytest.fixture
def superfam_name():
    """
    Dummy name used for Revise_Anno constructor
    """
    return 'test_superfam_set'


@pytest.fixture
def TEData_TestObj_SingleC_SingleElongate_SupF():
    """
    Create the TransposonData test object for a single chromosome of data. This
    data has been specially modified for superfamily focused tests.

    Use case (single chromosome):
        Merge two TEs that overlap
    """
    te_file = 'tests/input_data/Test_SingleC_SingleElongate_Superfam_Revision.tsv'
    transposon_input_dataframe = pd.read_csv(te_file,
                                             header='infer',
                                             dtype={'Start': 'float32', 'Stop':
                                                    'float32', 'Length':
                                                    'float32'}, sep='\t')
    sample_tes = TransposonData(transposon_input_dataframe, 'Mock_Camarosa')
    return sample_tes


@pytest.fixture
def TEData_TestObj_SingleC_MultiElongate_SupF():
    """
    Create the TransposonData test object for a single chromosome of data. This
    data has been specially modified for superfamily focused tests.

    Use case (single chromosome):
        Situation where multiple TEs overlap with one another. Possible to
        string together multiple TEs.
    """
    te_file = 'tests/input_data/Test_SingleC_MultiElongate_Superfam_Revision.tsv'
    transposon_input_dataframe = pd.read_csv(te_file,
                                             header='infer',
                                             dtype={'Start': 'float32', 'Stop':
                                                    'float32', 'Length':
                                                    'float32'}, sep='\t')
    sample_tes = TransposonData(transposon_input_dataframe, 'Mock_Camarosa')
    return sample_tes


@pytest.fixture
def TEData_TestObj_SingleC_ConcOverlap_SupF():
    """
    Create the TransposonData test object for a single chromosome of data. This
    data has been specially modified for superfamily focused tests.

    Use case (single chromosome):
        A TE that is completely inside another TE.
        A TE that has the exact same Start and Stop value as another
        A TE that Starts on another's Stop.
    """
    te_file = 'tests/input_data/Test_SingleC_ConcurrentOverlap_Superfam_Revision.tsv'
    transposon_input_dataframe = pd.read_csv(te_file,
                                             header='infer',
                                             dtype={'Start': 'float32', 'Stop':
                                                    'float32', 'Length':
                                                    'float32'}, sep='\t')
    sample_tes = TransposonData(transposon_input_dataframe, 'Mock_Camarosa')
    return sample_tes


@pytest.fixture
def ReviseAnno_TestObj_SingleC_SingleElongate_SupF(
                         TEData_TestObj_SingleC_SingleElongate_SupF,
                         h5_cache_loc, logger_obj, superfam_name,
                         revised_te_annotation_loc):
    """
    Test whether or not Revise_Anno will produce the correct output when given
    data of TEs that should be merged once.

    This is a fixture which uses the pipeline. While pipeline may fail, and
    this fixture would not work, we need this fixture to create the data in
    order to compare the expected numbers vs the observed numbers.
    """
    revise_anno_obj = Revise_Anno(TEData_TestObj_SingleC_SingleElongate_SupF.data_frame,
                                  h5_cache_loc,
                                  superfam_name)
    revise_anno_obj.create_superfam()
    revise_anno_obj.save_for_dev(revise_anno_obj.updated_te_annotation,
                                 os.path.join(revised_te_annotation_loc,
                                              'SingleC_SingleElongate_Super_Revision.tsv')
                                 )
    return Revise_Anno._read(revise_anno_obj.superfam_cache_loc)


@pytest.fixture
def ReviseAnno_TestObj_SingleC_MultiElongate_SupF(
                         TEData_TestObj_SingleC_MultiElongate_SupF,
                         h5_cache_loc, logger_obj, superfam_name,
                         revised_te_annotation_loc):
    """
    Test whether or not Revise_Anno will produce the correct output when given
    data of TEs that should be merged multiple times.

    This is a fixture which uses the pipeline. While pipeline may fail, and
    this fixture would not work, we need this fixture to create the data in
    order to compare the expected numbers vs the observed numbers.
    """
    revise_anno_obj = Revise_Anno(TEData_TestObj_SingleC_MultiElongate_SupF.data_frame,
                                  h5_cache_loc,
                                  superfam_name)
    revise_anno_obj.create_superfam()
    revise_anno_obj.save_for_dev(revise_anno_obj.updated_te_annotation,
                                 os.path.join(revised_te_annotation_loc,
                                              'SingleC_MultiElongate_Super_Revision.tsv')
                                 )
    return Revise_Anno._read(revise_anno_obj.superfam_cache_loc)


@pytest.fixture
def ReviseAnno_TestObj_SingleC_ConcOverlap_SupF(
                         TEData_TestObj_SingleC_ConcOverlap_SupF,
                         h5_cache_loc, logger_obj, superfam_name,
                         revised_te_annotation_loc):
    """
    Test whether or not Revise_Anno will produce the correct output when given
    data of TEs that overlap completely with one another and or start/stop on
    the edges of one another.

    This is a fixture which uses the pipeline. While pipeline may fail, and
    this fixture would not work, we need this fixture to create the data in
    order to compare the expected numbers vs the observed numbers.
    """
    revise_anno_obj = Revise_Anno(TEData_TestObj_SingleC_ConcOverlap_SupF.data_frame,
                                  h5_cache_loc,
                                  superfam_name)
    revise_anno_obj.create_superfam()
    revise_anno_obj.save_for_dev(revise_anno_obj.updated_te_annotation,
                                 os.path.join(revised_te_annotation_loc,
                                              'SingleC_ConcurrentOverlap_Super_Revision.tsv')
                                 )
    return Revise_Anno._read(revise_anno_obj.superfam_cache_loc)


# -------------------------------------------------------------
# TEST VALUES
TRUE_SingleC_SingleE_SUPERFAM = [784, 1021, 912, 1694, 634,
                                 212, 150, 94, 1432, 947, 623,
                                 177, 385, 229, 2382, 131,
                                 189, 1170, 501, 351]

TRUE_SingleC_MultiE_SUPERFAM = [1412, 560, 1219, 2209, 212, 150, 94, 2452, 623,
                                177, 385, 2611, 131, 189, 1170, 501, 351]

TRUE_SingleC_ConcOverlap_SUPERFAM = [608, 792, 201]
# -------------------------------------------------------------
# TESTS


@pytest.mark.parametrize('true_values',
                         [TRUE_SingleC_SingleE_SUPERFAM])
def test_superfam_SingleC_SingleE(ReviseAnno_TestObj_SingleC_SingleElongate_SupF,
                                  true_values):
    """
    Test whether or not Revise_Anno will produce the correct output when given
    data of TEs that should be merged once.
    """
    observed = ReviseAnno_TestObj_SingleC_SingleElongate_SupF.Length.to_numpy()
    assert np.array_equal(observed, true_values)


@pytest.mark.parametrize('true_values',
                         [TRUE_SingleC_MultiE_SUPERFAM])
def test_superfam_SingleC_MultiE(ReviseAnno_TestObj_SingleC_MultiElongate_SupF,
                                 true_values):
    """
    Test whether or not Revise_Anno will produce the correct output when given
    data of TEs that should be merged multiple times.
    """
    observed = ReviseAnno_TestObj_SingleC_MultiElongate_SupF.Length.to_numpy()
    assert np.array_equal(observed, true_values)


@pytest.mark.parametrize('true_values',
                         [TRUE_SingleC_ConcOverlap_SUPERFAM])
def test_superfam_SingleC_ConcOverlap(ReviseAnno_TestObj_SingleC_ConcOverlap_SupF,
                                      true_values):
    """
    Test whether or not Revise_Anno will produce the correct output when given
    data of TEs that overlap completely with one another and or start/stop on
    the edges of one another.
    """
    observed = ReviseAnno_TestObj_SingleC_ConcOverlap_SupF.Length.to_numpy()
    assert np.array_equal(observed, true_values)


if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
