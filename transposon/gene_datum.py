#!/usr/bin/env python3

"""
Wrappers for input data, a single gene.

Used to provide a common interface and fast calculations with numpy.
"""

__author__ = "Michael Teresi, Scott Teresi"

import functools
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
        """Initialize

        Args:
            gene_dataframe (DataFrame): gene data frame.
            gene_id (str): a string of the gene name.

        """

        self._logger = logger or logging.getLogger(__name__)
        self.name = str(gene_id)
        self.start = gene_dataframe.Start[str(gene_id)]
        self.stop = gene_dataframe.Stop[str(gene_id)]
        self.length = gene_dataframe.Length[str(gene_id)]
        self.chromosome = gene_dataframe.Chromosome[str(gene_id)]

        # For the left side, window stop is one digit before gene start
        self.left_win_stop = self.start - 1
        # For the right side, window start is one digit after gene stop
        self.right_win_start = self.stop + 1

    def write(self, filename):
        """Write to disk."""

        raise NotImplementedError()

    def read(self, filename):
        """Read from disk."""

        raise NotImplementedError()

    @functools.lru_cache(maxsize=None)
    def win_length(self, window):
        """Returns window length which is window val + 1

        Used in calculations of density. The value 1 is added to the window
        size , given by the user, the start stop and step of windowing
        operations, in order to correctly calculate the length because
        [0, inf).

        Args:
            window (int): integer representing the current window value with
            which we are performing the windowing procedures

        Returns:
            win_length (int): integer representing the modified window value
            with which we will calculate density.
        """
        win_length = window + 1
        return win_length  # plus 1 b/c it's [0...N]

    @functools.lru_cache(maxsize=None)
    def left_win_start(self, window):
        """Calculate left window start and clip the window to non-negative
        values

        Args:
            window (int): integer representing the original window value

        Returns:
            left_win_start (int): integer representing the first base-pair
            postion of the left-hand window for a given gene
        """
        left_win_start = np.subtract(self.left_win_stop, window)
        left_win_start = np.clip(left_win_start, 0, None)
        return left_win_start

    @functools.lru_cache(maxsize=None)
    def right_win_stop(self, window):
        """Calculate the stop base-pair of the right-window

        Args:
            window (int): integer representing the original window value

        Returns:
            right_win_stop (int): integer representing the final base-pair
            position of the right hand window for a given gene
        """
        right_win_stop = np.add(self.right_win_start, window)
        return right_win_stop

    @property
    def start_stop_len(self):
        return (self.start, self.stop, self.length)

    @functools.lru_cache(maxsize=None)
    def validate_window(self, left_win_start, window_length):
        """Make sure clipped windows produce the correct value for density
        division

        Input arg of window_length serves as  divisor if no clipping past 0
        was performed in left_win_start(), if clipping was performed then
        the divisor is calculated by recalculating "window_length" from 0 to
        left_win_stop

        Args:
            left_win_start (int): integer value for the left window start value
            window_length (int): integer value for the window length

        Returns:
            divisor (int): Returns the relevant area for density calculation
        """
        divisor = window_length
        if left_win_start < 0:
            msg = "window_start is not 0 or a positive value"
            logger.critical(msg)
            raise ValueError(msg)
        if left_win_start == 0:
            divisor = self.left_win_stop - left_win_start + 1
        return divisor

    def divisor_left(self, window):
        """Left density divisor for a gene for a window.

        Left-hand calculations: Window length most of the time, in the edge cases where
        the window extends past 0, turning negative, we need to clip the value to 0
        and make sure the relevant area is between 0 and the left window stop, one
        base-pair before gene start.

        Args:
            window (int): number of base pairs given by the windowing
            operation, used to calculate the relevant area (divisor).

        Returns:
            divisor (int): the relevant area for normalization required for
            density. This is the denominator for the density calculation

        Example:
            200 BP Gene Start with a 500 BP window moving to the left, the
            window cannot go negative past 0, so it is clipped to 0, and
            relevant area calculated to be 200 (199 <-> 0, 199 - 0 + 1)
        """
        left_win_start = self.left_win_start(window)
        win_length = self.win_length(window)
        divisor = self.validate_window(left_win_start, win_length)
        return divisor

    def divisor_intra(self, window):
        """Intronic density divisor for a gene for a window

        Intronic calculations: Gene length, here we are considering intragenic
        TEs, so the relevant area is the gene length.

        Args:
            window (None): normally the # of base pairs given by the windowing
            operation, here for intragenic calculations it is None

        Returns:
            self.length (int): the relevant area for normalization required for
            density. This is the denominator for the density calculation. Gene
            length for intragenic calculations
        """
        if window is None:
            return self.length
        else:
            msg = "window must be None for intragenic calculations"
            self._logger.critical(msg)
            raise ValueError(msg)

    def divisor_right(self, window):
        """Right density divisor for a gene for a window.

        Value will be window length for each gene because there is no clipping.
        BP positions go from [0, inf).

        Args:
            window (int): number of base pairs given by the windowing
            operation, used to calculate the relevant area (divisor).

        Returns:
            divisor (int): the relevant area for normalization required for
            density. This is the denominator for the density calculation
        """
        divisor = self.win_length(window)
        return divisor

    def __repr__(self):
        """Printable representation."""

        info = """
               GeneDatum name: {self.name}
               GeneDatum chromosome: {self.chromosome}
               GeneDatum start_stop_len: {self.start_stop_len}
               """
        return info.format(self=self)
