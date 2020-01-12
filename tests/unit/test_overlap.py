#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unit test overlap calculations.
"""

__author__ = "Michael Teresi"

import pytest
import logging, coloredlogs

import numpy as np
import pandas as pd

from transposon.overlap import Overlap, OverlapData, OverlapReducer
from transposon.data import GeneData, TransposonData



if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience, run only these tests
