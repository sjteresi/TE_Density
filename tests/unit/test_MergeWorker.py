#!/usr/bin/env python3

"""
Unit test MergeWorker.
"""

__author__ = "Michael Teresi"

import logging
import os
import pytest
import tempfile

import coloredlogs
import numpy as np
import pandas as pd

from transposon.merge_data import MergeData
from transposon.merge_data import _MergeConfigSink, _MergeConfigSource, _SummationArgs
from transposon.transposon_data import TransposonData
from transposon.gene_data import GeneData
from transposon.overlap import OverlapData, OverlapWorker

@pytest.fixture()
def temp_dir():
    """Temporary directory."""

    with tempfile.TemporaryDirectory() as dir:
        yield dir

@pytest.fixture(scope="module")
def overlap_data():
    pass
    # scope=module b/c we can reuse the files b/c they are read only
