#!/usr/bin/env python

"""
Orchestrates the summation of many OverlapData into many MergeData.
"""

__author__ = "Michael Teresi"

from collections import namedtuple
import logging
import os
from functools import partial
import pprint

import h5py
import numpy as np

import transposon


class MergeWorker():
    """Orchestrates merging multiple OverlapData."""

    def __init__(self, overlaps, transposons, genome, windows, output_dir, logger):
        """Initialzer.

        Args:
            overlaps (iterable(str)): filepaths to OverlapData.
            transposons (str): path to TransposonData file.
            genes (str): path to GeneData file.
            windows (iterable(int)): input window values.
            output_dir (str): path to output results.
        """

        pass


    def scrape_chromosome_ids(self):
        pass


    # scrape the windows from the overlap files, pass to MergeData
    # scrape genes from the overlap files, pass to MergeData
