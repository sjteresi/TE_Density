# IDEAS
# Write a barebones tests of fake data to get an idea as to what kind of
# interface we need.

import pytest
import pandas as pd

# Use case:
# 1) Give me TE Density (2D) for all genes, for a specific TE type, for a
# specific window, for a specific direction
# It could be annoying to type in the combination of all 3 things all the
# time. Perhaps


class DensityData:
    def __init__(self):
        pass

    # fake_data = pd.DataFrame.from_dict(fake_data)

    def windows():
        return [500, 1000, 1500]  # TODO implement
        # actual: return

    def get_data(self, te_type_str, window_int, direction_str) -> pd.DataFrame:
        return pd.DataFrame()


@pytest.fixture
def fake_data():
    """
    Sample DensityData
    """
    return DensityData()


def test_get_windows(fake_data):
    """
    Get windows. Can I get all the windows given the data table (h5, pandas?)
    """
    pass


def test_get_data_for_specific_combo(fake_data):
    # Can I get an array of TE Density Data for
    # data = x.get_data(te_type, window, direction)
    # This would be for all genes
    # This would be 2D
    # This would a numpy array (list of floats?)

    # Would these be parameters?
    te_type = "LTR"  # Maybe need to also say this is an Order NOT SuperFamily?
    window = 500
    direction = "Upstream"

    answer = fake_data.get_data(te_type, window, direction)
    # fake_data_frame would be 2D, and have shape (number of genes,1)
    # type (string,float)
    # The way the data would look for a specific combo:
    # Gene_ID   | Float
    # Gene_1 | 0.6
    # Gene_2 | 0.9

    truth = {"Gene_ID": ["gene_1", "gene_2"], "Density": [0.6, 0.9]}
    truth = pd.DataFrame.from_dict(truth)
    assert answer.equals(truth)
