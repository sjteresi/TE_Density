#!/usr/bin/env python3

__author__ = "Michael Teresi"

"""
The no. base pairs that the transposable elements overlap the window for a gene.
"""

from enum import Enum, unique
import logging

import numpy as np

from transposon.data import TransposonData, GeneData


class Overlap(object):
    """Functions for calculating overlap."""

    @unique
    class Direction(Enum):
        """Where the overlap is calculated relevant to the gene."""

        LEFT = 0
        INTRA = 1
        RIGHT = 2

    @staticmethod
    def left(gene_datum, transposons, window):
        """Overlap to the left (aka downstream) of the gene.

        Number base pair overlap between the TE / gene on the relevant distance.

        Args:
            gene_datum (GeneDatam): the relevant gene.
            transposons (TransposonData): the transposable elements.
            window (int): no. base pairs from gene to calculate overlap.
        """

        w_len = gene_datum.win_length(window)
        w_0 = gene_datum.left_win_start(w_len)
        w_1 = gene_datum.left_win_stop()
        #w_l = validate_window(w0, transposons.start, w_len)
        lower_bound = np.maximum(w_0, transposons.starts)
        upper_bound = np.minimum(w_1, transposons.stops)
        te_overlaps = np.maximum(0, (upper_bound - lower_bound + 1))
        return te_overlaps

    @staticmethod
    def intra(gene_datum, transposons):
        """Overlap to the gene itself.

        Number base pair overlap between the TE / gene on the relevant distance.

        Args:
            gene_datum (GeneDatam): the relevant gene.
            transposons (TransposonData): the transposable elements.
            window (int): no. base pairs from gene to calculate overlap.
        """

        g_0 = gene_datum.start
        g_len = gene_datum.length
        g_1 = gene_datum.stop
        lower_bound = np.minimum(g_0, transposons.stops)
        upper_bound = np.maximum(g_1, transposons.starts)
        te_overlaps = np.maximum(0, (lower_bound - upper_bound + 1))
        return te_overlaps

    @staticmethod
    def right(gene_datum, transposons, window):
        """Overlap to the right (aka upstream) of the gene.

        Number base pair overlap between the TE / gene on the relevant distance.

        Args:
            gene_datum (GeneDatam): the relevant gene.
            transposons (TransposonData): the transposable elements.
            window (int): no. base pairs from gene to calculate overlap.
        """

        w_0 = gene_datum.right_win_start()
        w_len = gene_datum.win_length(window)
        w_1 = gene_datum.right_win_stop(w_len)
        lower_bound = np.maximum(w0, transposons.starts)
        upper_bound = np.minimum(w1, transposons.stops)
        te_overlaps = np.maximum(0, (upper_bound - lower_bound + 1))
        return te_overlaps


class OverlapData(object):
    """Contains the overlap between the genes and transposable elements."""

    PRECISION = np.uint32

    def __init__(self, gene_data, transposon_data, logger=None):
        """Initialize.

        Args:
            gene_data(GeneData): the genes.
            transposon_data(TransposonData): the transposons.
            windows(list(int)):
        """

        self._logger = logger or logging.getLogger(__name__)
        self._transposons = transposon_data
        self._genes = gene_data
        self._left = None
        self._intra = None
        self._right = None
        self._windows = None
        self._gene_names = None
        self._gene_name_2_idx = None
        self._window_2_idx = None

    def calculate(self, windows, gene_names, window_update=None):
        """Calculate the overlap for the genes and windows.

        N.B. the number of runs of this multi level loop is intended to be
        limited by the number of input genes / windows rather than refactoring
        into a multiprocess solution.

        Args:
            windows (list(int)): iterable of window sizes to use.
            gene_names (list(str)): iterable of the genes to use.
        """

        self._reset(windows, gene_names)
        transposons = self._transposons
        for gene_name in self._gene_names:
            gene = self._genes.get_gene(gene_name)
            g_idx = self._gene_name_2_idx[gene_name]
            self._intra[:, g_idx] = Overlap.intra(gene, transposons)
            for window in self._windows:
                left = Overlap.left(gene, transposons, window)
                right = Overlap.left(gene, transposons, window)
                w_idx = self._window_2_idx[window]
                self._left[:, g_idx, w_idx] = left
                self._right[:, g_idx, w_idx] = right
            if window_update is not None:
                window_update()

    def write(self, filename):
        """Write to disk."""
        raise NotImplementedError()

    def read(self, filename):
        """Read from disk."""
        raise NotImplementedError()

    def _filter_gene_names(self, gene_names):
        """Yield only the valid gene names."""

        for name in gene_names:
            # could use list comprehension but we need to log if it fails
            if name not in self._genes.names:
                self._logger.error(" invalid gene name: {}".format(name))
            else:
                yield name

    def _filter_windows(self, windows):
        """Yield only the valid windows."""

        windows = list(windows)
        for window in windows:
            if window < 0:
                self._logger.error(" invalid window:  {}".format(window))
            else:
                yield window

    def _reset(self, windows, gene_names):
        """Initialize overalap data; mutates self."""

        # TODO remove list of gene names, replace with generator, e.g. sum(1 for n in names)...
        self._gene_names = list(self._filter_gene_names(gene_names))
        self._windows = list(self._filter_windows(windows))
        self._gene_name_2_idx = {name: idx for idx, name in enumerate(self._gene_names)}
        self._window_2_idx = {win: idx for idx, win in enumerate(windows)}
        n_win = len(self._windows)
        n_gene = len(self._gene_names)
        n_te = self._transposons.number_elements
        self._left = np.zeros((n_te, n_gene, n_win), self.PRECISION)
        self._intra = np.zeros((n_te, n_gene), self.PRECISION)
        self._right = np.zeros((n_te, n_gene, n_win), self.PRECISION)
