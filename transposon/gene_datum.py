#!/usr/bin/env python3

"""
Wrappers for input data, a single gene.

Used to provide a common interface and fast calculations with numpy.
"""

__author__ = "Michael Teresi, Scott Teresi"

import logging
import numpy as np


class GeneDatum(object):
    """Wraps data frame containing one gene.

    Provides an interface, attribute access, and to/from disk functionality.

    Note the numpy views are not necessarily a no-copy (SEE pandas.DataFrame.to_numpy).

    Expects certain column identifiers (SEE self.__init__).
    GeneData subclasses should conform to these column names or redefine the properties.
    """

    def __init__(self, gene_dataframe, gene_id, logger=None):

        self._logger = logger or logging.getLogger(__name__)
        self.name = str(gene_id)
        self.start = gene_dataframe.Start[str(gene_id)]
        self.stop = gene_dataframe.Stop[str(gene_id)]
        self.length = gene_dataframe.Length[str(gene_id)]
        self.chromosome = gene_dataframe.Chromosome[str(gene_id)]

        # For the left side, window stop is one digit before gene start
        self.left_win_stop = np.subtract(self.start, 1)
        # For the right side, window start is one digit after gene stop
        self.right_win_start = np.add(self.stop, 1)

    def write(self, filename):
        """Write to disk."""

        raise NotImplementedError()

    def read(self, filename):
        """Read from disk."""

        raise NotImplementedError()

    def win_length(self, window):
        # SCOTT why the plus 1?
        # SCOTT why is this not an off by one error?
        # If it's just semantics then it would be helpful if this func was unecessary

        # MICHAEL I am going to keep the function for now because this will be
        # necessary for the density calculations later when you need to divide
        # by the window length, but for now I have simplified the notation for
        # the window calculations so that they are clearer.
        return np.add(window, 1)

    def left_win_start(self, window):
        win_start = np.subtract(self.left_win_stop, window)
        win_start = np.clip(win_start, 0, None)
        return win_start

    def right_win_stop(self, window):
        return np.add(self.right_win_start, window)

    @property
    def start_stop_len(self):
        return (self.start, self.stop, self.length)

    def __repr__(self):
        """Printable representation."""

        info = """
               GeneDatum name: {self.name}
               GeneDatum chromosome: {self.chromosome}
               GeneDatum start_stop_len: {self.start_stop_len}
               """
        return info.format(self=self)
