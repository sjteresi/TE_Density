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
        print(g0, g1)
        genes[gi, :] = [g0, g1]
    return genes



#@pytest.mark.parametrize("window", WINDOWS)
def test_rho_only_inside_congruent():
    """Does ONLY INSIDE work when TE is the same size as gene?"""

    genes = gene_array()
    transposons = genes

    for t_i in range(transposons.shape[0]):
        hit = is_inside_only(genes, transposons[t_i, :])
        assert hit[t_i]
        hit[t_i] = False
        assert np.all(np.invert(hit))



def test_rho_only_inside_subset():
    """Does ONLY INSIDE work when TE is a subset of gene?"""
    genes = gene_array()
    transposons = np.array([(3,9), (13,18), (22, 27), (35,37)])

    for t_i in range(transposons.shape[0]):
        hit = is_inside_only(genes, transposons[t_i, :])
        assert hit[t_i]
        hit[t_i] = False
        assert np.all(np.invert(hit))


def test_rho_only_inside_on_window():
    """Does ONLY INSIDE work when TE is on window?"""
    pass



if __name__ == "__main__":

    pytest.main(['-s', __file__])  # for convience
