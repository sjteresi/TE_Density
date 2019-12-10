#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Unit test overlap calculations.
"""

__author__ = "Michael Teresi"

import pytest

import logging, coloredlogs

from transposon.overlap import Overlap, OverlapData, OverlapMerger

if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience, run only these tests
