#!/usr/bin/env python3

"""
Collates density results.

Each density is independent, so the processing order is not important.
Densities are summed for one gene/window pair wrt one TE type.
DensityAccumulator receives the densities and produces
"""

__author__ = "Michael Teresi, Scott Teresi"

from multiprocessing import Process, Queue

import numpy as np


class LoggerAccumulator(object):
    """Collates logs from the workers."""

    # SEE https://gist.github.com/JesseBuesking/10674086
    pass


class DensityAccumulator(object):
    """Sums density values provided results."""

    # NOTE maybe use a multiprocessing manager to do the merging?
    # https://stackoverflow.com/questions/8640367/python-manager-dict-in-multiprocessing

    def __init__(self, genes, transposons, window):
        """Initialize.

        Args:
            genes (transposon.GeneData): gene container.
            transposons (transposon.TransposonData): transposon container.
            windows (list(int)): the windows to accumulate.

        """

        self.gene_data = genes
        n_genes = sum(1 for n in genes.names)
        self.te_data = transposons
        self.window = window
        self.left_overlap = np.zeros()


    def accumulate(self, density_result):
        pass
