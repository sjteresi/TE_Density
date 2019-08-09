#!/usr/bin/env python3

"""
Unit test density calculations.

"""

__author__ = "Michael Teresi"

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

from transposon.density import is_inside_only
from transposon.density import in_left_window_only

WINDOWS = [2, 4, 8]

def gene_array(n_genes=4, width=10, seperation=0):
    """Return an Nx2 (N x (start, stop)) dummy gene.


    """

    genes = np.zeros((n_genes, 2))
    g0 = 0
    g1 = 0
    for gi in range(n_genes):
        g0 = g1 + seperation
        g1 = g0 + width
        #print(g0, g1)
        genes[gi, :] = [g0, g1]
    return genes


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


if __name__ == "__main__":

    pytest.main(['-s', __file__])  # for convience
