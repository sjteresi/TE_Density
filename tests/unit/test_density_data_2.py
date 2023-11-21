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

from transposon.transposon_data import TransposonData
from transposon.gene_data import GeneData
from transposon.test_utils import temp_dir, temp_h5_file

from transposon.density_data_2 import DensityData


N_TRANSPOSONS = 4
WINDOWS = [10, 20, 30, 40]
CHROME_ID = "fake_chromosome_id"
LOGGER = logging.getLogger(__name__)
coloredlogs.install(level=logging.DEBUG)


@pytest.fixture
def gene_data():
    """Default GeneData instance."""

    # TODO is a fake, not mock
    return GeneData.mock(genome_id=CHROME_ID)


@pytest.fixture
def te_data():
    """Default TransposonData instance."""

    # TODO is a fake, not mock
    return TransposonData.mock(chromosome=CHROME_ID)


@pytest.fixture
def density(gene_data, te_data):

    # MAGIC gene_0 is the format in GeneData.mock,
    # TODO add an iterator to get all the datums?
    rho = DensityData(gene_data, te_data, WINDOWS)
    with rho as handle:
        yield handle


def test_density_init(density):
    """Does the DensityData instance enter / exit?"""

    pass
