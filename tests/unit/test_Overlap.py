#!/usr/bin/env python3

"""
Unit test Overlap.
"""

__author__ = "Scott Teresi, Michael Teresi"

import pytest

import numpy as np
import pandas as pd

from transposon.gene_data import GeneData
from transposon.overlap import Overlap

# NOTE these should be similar to the density tests, you are refactoring them sto that
# the scope is smaller so that we can do the pseudo split merge pattern more easily
# TODO add mock data
# TODO add left tests
# TODO add intra tests
# TODO add right tests

if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
