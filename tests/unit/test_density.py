#!/usr/bin/env python3

"""
Unit test density calculations.

"""

__author__ = "Michael Teresi, Scott Teresi"

import pytest

import os
import time
import argparse
import logging
#from multiprocessing import Process
#import multiprocessing
#from threading import Thread

import numpy as np
import pandas as pd

from transposon.density import rho_intra
from transposon.density import is_inside_only
from transposon.density import in_left_window_only
from transposon.density import rho_left_window
from transposon.density import rho_right_window

from transposon.data import GeneData, TransposonData

WINDOWS = [2, 4, 8]


def mock_gene_data(start_stop=np.array([[0, 9], [10, 19], [20, 29]]) ):
    """Gene data for given the start/stop indices.

    Args:
        start_stop (np.array): N gene x (start_idx, stop_idx).
    """

    n_genes = start_stop.shape[0]
    data = []
    for gi in range(n_genes):
        g0 = start_stop[gi, 0]
        g1 = start_stop[gi, 1]
        gL = g1 - g0
        name = "gene_{}".format(gi)
        datum = [name, g0, g1, gL]
        data.append(datum)

    frame = pd.DataFrame(data, columns=['Gene_Name', 'Start', 'Stop', 'Length'])
    frame.set_index('Gene_Name', inplace=True)
    return GeneData(frame)

def mock_te_data(start_stop):
    """Transposon data for given the start/stop indices.

    Creates one

    Args:
        start_stop (np.array): N gene x (start_idx, stop_idx).
    """

    n_genes = start_stop.shape[0]
    data = []
    family = "Family_0"  # FUTURE may want to parametrize family name later
    # NB overall order is not important but the names are
    columns = ['Start', 'Stop', 'Length', 'Family', 'SubFamily']
    for gi in range(n_genes):
        g0 = start_stop[gi, 0]
        g1 = start_stop[gi, 1]
        gL = g1 - g0
        # FUTURE may want to parametrize sub family name later
        subfam_suffix = "A" if gi%2 else "B"
        subfamily = "SubFamily_{}".format(subfam_suffix)
        datum = [g0, g1, gL, family, subfamily]
        data.append(datum)

    frame = pd.DataFrame(data, columns=columns)
    return TransposonData(frame)

def test_rho_only_inside_congruent():
    """Does ONLY INSIDE work when TE is the same size as gene?"""

    genes = np.array([[0, 9], [10, 19], [20, 29]])
    transposons = genes  # the TEs start stop at the same places
    n_genes = genes.shape[0]  # MAGIC NUMBER each row is a gene

    for t_i in range(transposons.shape[0]):
        #print(transposons[t_i, :])
        hit = is_inside_only(genes, transposons[t_i, :])
        expected_hit = np.zeros((n_genes))
        # MAGIC NUMBER since TE == gene, only the curent TE should overlap
        expected_hit[t_i] = True
        assert np.all(expected_hit == hit)


def test_rho_only_inside_subset():
    """Does ONLY INSIDE work when TE is purely a subset of gene?"""
    genes = gene_array()
    transposons = np.array([(3,9), (13,18), (22, 27), (35,37)])

    for t_i in range(transposons.shape[0]):
        hit = is_inside_only(genes, transposons[t_i, :])
        assert hit[t_i]
        hit[t_i] = False
        assert np.all(np.invert(hit))

def test_rho_only_inside_one_on_gene():
    """Does ONLY INSIDE work when one end of the TE is on edge of a gene?"""
    genes = gene_array()
    transposons = genes
    transposons[:,1] = np.subtract(transposons[:,1], 1)

    for t_i in range(transposons.shape[0]):
        hit = is_inside_only(genes, transposons[t_i, :])
        assert hit[t_i]
        hit[t_i] = False
        assert np.all(np.invert(hit))

def test_rho_in_left_window_only():
    # TODO should rename test
    """Does ONLY INSIDE work when TE is on window?"""
    genes = gene_array()
    transposons = np.copy(genes)
    transposons[:,1] = np.subtract(transposons[:,1], 11) # put TEs in left
    transposons[:,0] = np.subtract(transposons[:,0], 2)
    window = 5

    for t_i in range(transposons.shape[0]):
        #print(genes)
        #print(transposons[t_i, :])
        hit = in_left_window_only(genes, transposons[t_i, :], window)
        assert hit[t_i]
        hit[t_i] = False
        assert np.all(np.invert(hit))

def test_intra_density_congruent():
    """Does the intra density return 1 when the TE is the same as gene?"""

    genes = np.array([[0, 9], [10, 19], [20, 29]])
    transposon = genes[1, :]
    g0 = genes[:,0]
    g1 = genes[:,1]
    gl = g1 - g0 # MAGIC NUMBER the data uses an inclusive stop
    t0 = transposon[0]
    t1 = transposon[1]
    rhos = rho_intra(g0, g1, gl, t0, t1)
    expected_rhos = np.array([0, 1.0, 0])
    assert np.all(rhos == expected_rhos)

def test_intra_density_partial():
    """Does the intra density return 1 when the TE is the same as gene?"""
    genes = np.array([[0, 10], [200, 300], [1000, 1500]])
    transposon = np.array([225,275])
    g_start = genes[:,0]
    g_stop = genes[:,1]
    g_length = g_stop - g_start # MAGIC NUMBER the data uses an inclusive stop
    t_start = transposon[0]
    t_stop = transposon[1]
    rhos = rho_intra(g_start, g_stop, g_length, t_start, t_stop)
    expected_rhos = np.array([0, 0.5, 0])
    assert np.all(rhos == expected_rhos)

@pytest.mark.parametrize("window", [100, 200])
@pytest.mark.parametrize("extend", [50, 150, 250])
def test_rho_left_window(window, extend):
    """
    Does the correct window density return when TE is in window?
    SCOTT: I have tested this and it works with TEs that are:
        Completely spanning the gene.
        Starting partially in the gene
        Only in the window
        It successfully captures just the part in the window
    """

    # TODO SCOTT, good job testing the different cases
    # pls encode that into pytest parameterizations to automate it (pytest.mar.parametrize)
    # TODO make sure the edge cases work when the window should be negative and
    # reset to 0 for that instance
    start_stop = np.array([[1000, 2000], [5000, 6500], [2225, 3000]])
    genes = mock_gene_data(start_stop)
    # move the TE start to the left using the 'extend' param
    transposons = mock_te_data(start_stop[0,:] - np.array([[extend, 0]]))
    gene_name = list(genes.names)[0]  # MAGIC NUMBER just use the first one
    rhos = rho_left_window(genes, gene_name, transposons, window)
    expected = min(extend / window, 1)
    expected_rhos = np.array([expected])
    print(" expected {}".format(expected_rhos))
    print(" rhos {}".format(rhos))
    assert np.all(rhos == expected_rhos)

@pytest.mark.parametrize("window", [100, 200])
@pytest.mark.parametrize("extend", [50, 150, 250])
def test_rho_right_window(window, extend):
    """
    Does the correct window density return when TE is in window?
    SCOTT: I have tested this and it works with TEs that are:
        Completely spanning the gene.
        Starting partially in the gene
        Only in the window
        It successfully captures just the part in the window
        I am pretty sure this test and function is fully working
    MIKE: Scott, you can use a pytest fixture to parametrize these multiple
        conditions rather than (presumably) changing the inputs by hand
    """

    # TODO SCOTT, good job testing the different cases
    # pls encode that into pytest parameterizations to automate it (pytest.mar.parametrize)
    start_stop = np.array([[1000, 2000], [5000, 6500], [2225, 3000]])
    genes = mock_gene_data(start_stop)
    # move the TE start to the right using the 'extend' param
    transposons = mock_te_data(start_stop[0,:] + np.array([[0, extend]]))
    gene_name = list(genes.names)[0]  # MAGIC NUMBER just use the first one
    rhos = rho_right_window(genes, gene_name, transposons, window)
    expected = min(extend / window, 1)
    expected_rhos = np.array([expected])
    assert np.all(rhos == expected_rhos)


if __name__ == "__main__":

    pytest.main(['-s', __file__])  # for convience
