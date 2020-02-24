#!/usr/bin/env python3

"""
Wrappers for the input data.

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
        self.unique_genes = gene_dataframe.index.unique()

    def write(self, filename):
        """Write to disk."""

        # NOTE consider for map reduce?
        raise NotImplementedError()

    def read(self, filename):
        """Read from disk."""

        # NOTE consider for map reduce?
        raise NotImplementedError()

    def get_gene(self, gene_id):
        """Return a GeneDatum for the gene identifier."""
        return GeneDatum(self.data_frame, gene_id)

    @property
    def names(self):
        """Yields the names for each gene."""

        return (name for name in self._names)

    def __repr__(self):
        info = """
               Wrapped Gene DataFrame: {self.data_frame}
               """
        return info.format(self=self)


class GeneDatum(object):
    """Wraps a single gene data frame.

    Provides attribute access for a single gene.
    """

    def __init__(self, gene_dataframe, gene_id):
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

    def write(self, filename):
        """Write to disk."""

        raise NotImplementedError()

    def read(self, filename):
        """Read from disk."""

        raise NotImplementedError()

    @property
    def number_elements(self):
        """The number of transposable elements."""

        # TODO verify
        return self.indices.shape[0]  # MAGIC NUMBER it's one column

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
            self._logger.critical(msg)
            raise ValueError(msg)

        # TODO verify that this isn't weird
        length = self.lengths.shape
        if start != length:
            msg = ("Input TE missing fields: starts.shape {}  != lengths.shape {}"
                   .format(start, stop))
            self._logger.critical(msg)
            raise ValueError(msg)

    def __repr__(self):
        """Printable representation."""

        info = "Wrapped TE DataFrame: {self.data_frame}"
        return info.format(self=self)

    def __add__(self, other):
        """Combine transposon data."""

        # SCOTT this may be useful for testing, would you take a look?
        # for example, it's easy to make a mocked TE data for a specific family
        # so if we wanted to handle multiple sub/families we could parametrize
        # the function to produce one TE wrt sub/family and combine them after
        raise NotImplementedError()
