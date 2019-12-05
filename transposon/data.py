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

__author__ = "Michael Teresi"


class GeneData(object):
    """Wraps a gene data frame.

    Provides attribute access for a single gene.
    """

    def __init__(self, gene_dataframe):
        """Initialize.

        Args:
            gene_dataframe (DataFrame): gene data frame.
        """

        self.data_frame = gene_dataframe
        self._names = self.data_frame.index

    def start(self, gene_id):
        """The start index for the given gene."""

        return self.data_frame.Start[str(gene_id)]

    def stop(self, gene_id):
        """The stop index for the given gene."""

        return self.data_frame.Stop[str(gene_id)]

    def length(self, gene_id):
        """The number of base pairs for the given gene."""

        return self.data_frame.Length[str(gene_id)]

    def start_stop_len(self, gene_id):
        """Helper to return start, stop, and length for the given gene."""

        return (self.start(gene_id), self.stop(gene_id), self.length(gene_id))

    def chromosome(self, gene_id):
        """The chromosome identifier for the given gene."""

        return self.data_frame.Chromosome[str(gene_id)]

    @property
    def names(self):
        """Yields the names for each gene."""

        return (name for name in self._names)


class TransposonData(object):
    """Wraps a transposable elements data frame.

    Provides attribute access for the transposons as a whole.
    Returns numpy arrays to the underlying data, intended to be read-only views.
    NB: the 'copy' flag of pandas.DataFrame.to_numpy does not necessarily
        ensure that the output is a (no-copy) view.
    """

    def __init__(self, transposable_elements_dataframe):
        """Initialize.

        Args:
            gene_dataframe (DataFrame): transposable element data frame.
        """
        self._data_frame = transposable_elements_dataframe
        self.indices = self._data_frame.index.to_numpy(copy=False)
        self.starts = self._data_frame.Start.to_numpy(copy=False)
        self.stops = self._data_frame.Stop.to_numpy(copy=False)
        self.lengths = self._data_frame.Length.to_numpy(copy=False)
        self.orders = self._data_frame.Family.to_numpy(copy=False)
        self.superfamilies = self._data_frame.SubFamily.to_numpy(copy=False)
        # TODO: change input data frame names, Family-->Order, SubFamily-->SuperFamily

        # TODO run the validate_chromosome
        self.chromosome_unique_id = ""  # TODO get the 'Chromosome' entry
        self.order_superfamily_set = set()  # TODO get the set

    def validate_chromosome(self):
        """Checks if the given data is specific to one chromosome."""

        raise NotImplementedError()

    def scrape_order_superfamily_set(self):
        """Set of order / superfamilies for the transposable elements."""

        raise NotImplementedError()
