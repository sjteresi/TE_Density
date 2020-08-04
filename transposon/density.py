#!/usr/bin/env python3

"""
Calculate transposable element density.
"""

__author__ = "Scott Teresi, Michael Teresi"

import argparse
import os

import logging
import coloredlogs
import numpy as np
import pandas as pd
from tqdm import tqdm
from configparser import ConfigParser
import sys
import time

from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.overlap import OverlapWorker
from transposon.replace_names import te_annot_renamer
from transposon.verify_cache import (verify_chromosome_h5_cache,
    verify_TE_cache, verify_gene_cache, revise_annotation)
from transposon.revise_annotation import Revise_Anno



# TODO remove
def rho_left_window(gene_data, gene_name, transposon_data, window):
    """Density to the left (downstream) of a gene.
    When TE is between gene and window

    Relevant things: No. TE bases / window
    GeneStart so that we may calculate the WindowStart
    WindowStart is unique to each gene. Calculated via Gstart - Window
    WindowStart is the leftmost window

    Args:
        window (int): integer value of the current window
        gene_data (transponson.data.GeneData): gene container
        gene_name (hashable): name of gene to use
        transposon_data (transponson.data.TransposonData): transposon container
    """
    transposon_data.check_shape()
    gene_datum = gene_data.get_gene(gene_name)
    g_start, g_stop, g_length = gene_datum.start_stop_len

    # Define windows
    win_length = gene_datum.win_length(window)
    win_start = gene_datum.left_win_start(window)
    win_stop = gene_datum.left_win_stop
    win_length = validate_window(win_start, g_start, win_length)

    # Set bounds and perform density calculation
    lower_bound = np.maximum(win_start, transposon_data.starts)
    upper_bound = np.minimum(win_stop, transposon_data.stops)
    te_overlaps = np.maximum(0, (upper_bound - lower_bound + 1))
    densities = np.divide(
        te_overlaps,
        win_length,
        out=np.zeros_like(te_overlaps, dtype='float')
    )
    check_density_shape(densities, transposon_data)
    return densities


# TODO remove
def rho_intra(gene_data, gene_name, transposon_data):
    """Intra density for one gene wrt transposable elements.

    The relevant ares is the gene for intra density.

    Args:
        gene_data (transponson.data.GeneData): gene container
        gene_name (hashable): name of gene to use
        transposon_data (transponson.data.TransposonData): transposon container
    """
    transposon_data.check_shape()
    gene_datum = gene_data.get_gene(gene_name)
    g_start, g_stop, g_length = gene_datum.start_stop_len
    lower = np.minimum(g_stop, transposon_data.stops)
    upper = np.maximum(g_start, transposon_data.starts)
    te_overlaps = np.maximum(0, (lower - upper + 1))

    densities = np.divide(
        te_overlaps,
        g_length,
        out=np.zeros_like(te_overlaps, dtype='float'),
        where=g_length != 0
    )

    check_density_shape(densities, transposon_data)
    return densities


# TODO remove
def rho_right_window(gene_data, gene_name, transposon_data, window):
    """Density to the right (upstream) of a gene.
    When TE is between gene and window

    Relevant things: No. TE bases / window
    GeneStop so that we may calculate the WindowStop
    WindowStop is unique to each gene. Calculated via Gstop + Window
    WindowStop is the rightmost window

    Args:
        window (int): integer value of the current window
        gene_data (transponson.data.GeneData): gene container
        gene_name (hashable): name of gene to use
        transposon_data (transponson.data.TransposonData): transposon container
    """
    transposon_data.check_shape()
    gene_datum = gene_data.get_gene(gene_name)
    g_start, g_stop, g_length = gene_datum.start_stop_len

    # Define windows
    win_length = gene_datum.win_length(window)
    win_start = gene_datum.right_win_start
    win_stop = gene_datum.right_win_stop(window)

    # Set bounds and perform density calculation
    lower_bound = np.maximum(win_start, transposon_data.starts)
    upper_bound = np.minimum(win_stop, transposon_data.stops)
    te_overlaps = np.maximum(0, (upper_bound - lower_bound + 1))
    densities = np.divide(
        te_overlaps,
        win_length,
        out=np.zeros_like(te_overlaps, dtype='float')
    )
    check_density_shape(densities, transposon_data)
    return densities




