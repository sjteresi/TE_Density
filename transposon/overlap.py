#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Michael Teresi"

"""
The number of base pairs that the transposable elements overlap the window.

Used to produce the percent overlap of a TE wrt gene & TE order|superfamily.
The overlap is summed wrt the TE / gene pair for every oder|superfamily.

# NOTE collate logs between processes?
# SEE https://gist.github.com/JesseBuesking/10674086
"""

import numpy as np

from transposon.data import TransposonData, GeneData

PRECISION = np.uint32

class Overlap(object):
    """Functions for calculating overlap."""

    @static_method
    def left(gene_data, gene_name, transposon_data, window):
        """Overlap to the left (aka downstream) of the gene.

        Number base pair overlap between the TE / gene on the relevant distance.

        Args:
            gene_data (GeneData): the genes.
            gene_name (str): gene the calculate with respect to.
            transposon_data (TransposonData): the transposons.
            window (int): no. base pairs from gene to calculate overlap.
        """

        raise NotImplementedError()

    @static_method
    def intra(gene_data, gene_name, transposon_data):
        """Overlap to the gene itself.

        Number base pair overlap between the TE / gene on the relevant distance.

        Args:
            gene_data (GeneData): the genes.
            gene_name (str): gene the calculate with respect to.
            transposon_data (TransposonData): the transposons.
        """

        raise NotImplementedError()

    @static_method
    def left(gene_data, gene_name, transposon_data, window):
        """Overlap to the right (aka upstream) of the gene.

        Number base pair overlap between the TE / gene on the relevant distance.

        Args:
            gene_data (GeneData): the genes.
            gene_name (str): gene the calculate with respect to.
            transposon_data (TransposonData): the transposons.
            window (int): no. base pairs from gene to calculate overlap.
        """

        raise NotImplementedError()


class OverlapData(object):
    """Contains no. of base pairs overlapping between the TEs and a gene.

    Tracks overlap for one *and only one* gene on one *and only one* chromosome.
    Designed for one split of the gene data: one gene, all TE, all windows.
    """

    def __init__(self, chromosome, gene, transposons, windows):
        """Initialize.

        Args:
            chromosome(str): chromosome name, e.g. 'Fvb1-1'.
            gene(str): gene name, e.g. 'maker-Fvb1-1-snap-gene-0.15'.
            transposons(TransposonData):

        """

        self._chromosome_name = str(chromosome)
        self._gene_name = str(gene)
        self._te_categories = []  # TODO get from TransposonData
        self._te_category_to_idx = {}  # TODO get from TransposonData
        self._left = None
        self._intra = None
        self._right = None

    @property
    def chromosome_name(self):
        """The name of the chromosome that the gene / TEs were on."""

        return self._chromosome_name

    @property
    def gene_name(self):
        """The name of the gene the TEs were compared to."""

        return self._gene_name

    @class_method
    def calculate(self, transposons, genes, gene_name, windows):



    def __add__(self, other):
        """Sum overlapped counters with another instance."""

        # NOTE may want to remove this in order to simplify
        # leave it to the reducer?
        raise NotImplementedError()

    def write(self):
        """Write to file."""

        raise NotImplementedError()

    def read(self):
        """Read from file."""

        raise NotImplementedError()

    @static_method
    def _zero_array(te_categories, window_list):
        """Return an empty array for storing the overlap."""

        # MAGIC NUMBER there are three overlaps: left, intra, right
        return np.zero(len(te_categories), 3, len(window_list))

class OverlapReducer(object):
    """Reduces the overlap values into one data structure.

    Reduces the overlapped values by summing wrt the TE order|superfamily.
    """

    def __init__(self, transposons, window_list):
        """Initialize.

        Args:
            transposons (TransposonData): the transposon data container.
            window (list(int)): the windows to calculate overlap.
        """

        self._transposons = transposons
        # n_orders = TODO add func to get order set
        # n_supfam = TODO add func to get superfamily set
        self._n_te_fam = n_orders + n_supfam
        self._window_list = window_list
        self._gene_to_overlap = {}

    def __add__(self, other):
        """Sum overlapped counters with another instance."""

        raise NotImplementedError()

    def write(self):
        """Write to file."""

        raise NotImplementedError()

    def read(self):
        """Read from file."""

        raise NotImplementedError()





    def accumulate(self, gene_name, )
    LEFT OFF HERE

    @static_method
    def overlap_matrix(self):
        """A zero array to store overlap values."""

        raise NotImplementedError()  # should be N_te_order_fam x 3 x len(window)
