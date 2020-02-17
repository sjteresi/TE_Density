#!/usr/bin/env python3

"""
Wrappers for the data.

Design:
    Downstream work operates with one gene combined with many TEs, so:
        - attribute access is provided for genes wrt one gene id.
        - attribue access is provided for TEs wrt all TEs at once.
    Numpy arrays are provided for the TE series for speed calculating density.

Future:
    - could define class level keys to use with getattr, then replace the
      dot notation for accessing the columns, to make this inheritable
"""

__author__ = "Michael Teresi, Scott Teresi"

import logging

import numpy as np


class GeneData(object):
    """Wraps a gene data frame.

    Provides attributes for numpy views of the data frame entries.
    Delineates the interface for the data frame
    """

    def __init__(self, gene_dataframe):
        """Initialize.

        Args:
            gene_dataframe (DataFrame): gene data frame.
        """

        self.data_frame = gene_dataframe
        self._names = self.data_frame.index
        self.starts = self.data_frame.Start.to_numpy(copy=False)
        self.stops = self.data_frame.Stop.to_numpy(copy=False)
        self.lengths = self.data_frame.Length.to_numpy(copy=False)
        self.chromosomes = self.data_frame.Chromosome.to_numpy(copy=False)

    def get_gene(self, gene_id):
        """Return a GeneDatum for the gene identifier."""
        return GeneDatum(self.data_frame, gene_id)

    @property
    def names(self):
        """Yields the names for each gene."""

        return (name for name in self._names)


class GeneDatum(object):
    """Wraps a single gene data frame.

    Provides attribute access for a single gene.
    """

    # SCOTT implement, __init__, add 'public' attributes for the columns
    # basically, use the [] syntax from the functions in GeneData that you are deleting

    def __init__(self, gene_dataframe, gene_id):
        self.start = gene_dataframe.Start[str(gene_id)]
        self.stop = gene_dataframe.Stop[str(gene_id)]
        self.length = gene_dataframe.Length[str(gene_id)]
        self.chromosome = gene_dataframe.Chromosome[str(gene_id)]

    def win_length(self, window):
        return np.add(window, 1)

    def left_win_start(self, win_length):
        win_start = np.subtract(self.start, win_length)
        win_start = np.clip(win_start, 0, None)
        return win_start

    def left_win_stop(self):
        # For the left side, window stop is one digit before gene start
        return np.subtract(self.start, 1)

    def right_win_start(self):
        # For the right side, window start is one digit after gene stop
        return np.add(self.stop, 1)

    def right_win_stop(self, win_length):
        return np.add(self.stop, win_length)

    @property
    def start_stop_len(self):
        return (self.start, self.stop, self.length)


class TransposonData(object):
    """Wraps a transposable elements data frame.

    Provides attribute access for the transposons as a whole.
    Returns numpy arrays to the underlying data, intended to be views.
    NB: the 'copy' flag of pandas.DataFrame.to_numpy does not necessarily
        ensure that the output is a (no-copy) view.
    """

    # FUTURE filter out irrelevant TEs for a window by comparing the
    # max gene stop with the TE start? for optimization
    # FUTURE split one TE file into multiple if memory becomes and issue?

    def __init__(self, transposable_elements_dataframe, logger=None):
        """Initialize.

        Args:
            transposable_elements (DataFrame): transposable element data frame.
        """

        self._logger = logger or logging.getLogger(__name__)
        self.data_frame = transposable_elements_dataframe
        self.indices = self.data_frame.index.to_numpy(copy=False)
        self.starts = self.data_frame.Start.to_numpy(copy=False)
        self.stops = self.data_frame.Stop.to_numpy(copy=False)
        self.lengths = self.data_frame.Length.to_numpy(copy=False)
        self.orders = self.data_frame.Order.to_numpy(copy=False)
        self.superfamilies = self.data_frame.SuperFamily.to_numpy(copy=False)
        self.chromosomes = self.data_frame.Chromosome.to_numpy(copy=False)

    @property
    def number_elements(self):
        """The number of transposable elements."""

        # TODO verify
        return self.indices.size[0]  # MAGIC NUMBER it's one column

    def check_shape(self):
        """Checks to make sure the columns of the TE data are the same size.

        If the shapes don't match then there are records that are incomplete,
            as in an entry (row) does have all the expected fields (column).
        """

        start = self.starts.shape
        stop = self.stops.shape
        if start != stop:
            msg = ("Input TE missing fields: starts.shape {}  != stops.shape {}"
                   .format(start, stop))
            logger.critical(msg)
            raise ValueError(msg)

        length = transposon_data.lengths.shape
        if start != length:
            msg = ("Input TE missing fields: starts.shape {}  != lengths.shape {}"
                   .format(start, stop))
            logger.critical(msg)
            raise ValueError(msg)

    def __add__(self, other):
        """Combine transposon data."""

        # SCOTT this may be useful for testing, would you take a look?
        # for example, it's easy to make a mocked TE data for a specific family
        # so if we wanted to handle multiple sub/families we could parametrize
        # the function to produce one TE wrt sub/family and combine them after
        raise NotImplementedError()

class OverlapData(object):
    """Stores the number of base pairs that overlap the gene name / TE.

    This data is used to calculate density.
    First, the number of base pairs that overlap the gene name / TE accumulates.
    Second, the accumulated overlaps are divided by the relevant area,
        the window size for left/right and the gene length for intra.
    """

    def __init__(self, genes, transposons, window):
        """Initialize.

        Args:
            genes (GeneData): the gene data container.
            transposons (TransposonData): the transposon data container.
            window (int): the window to accumulate to.
        """
