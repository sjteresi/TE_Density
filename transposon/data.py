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
    Returns numpy arrays to the underlying data, intended to be views.
    NB: the 'copy' flag of pandas.DataFrame.to_numpy does not necessarily
        ensure that the output is a (no-copy) view.
    """

    def __init__(self, transposable_elements_dataframe):
        """Initialize.

        Args:
            gene_dataframe (DataFrame): transposable element data frame.
        """
        self.data_frame = transposable_elements_dataframe
        self.indices = self.data_frame.index.to_numpy(copy=False)
        self.starts = self.data_frame.Start.to_numpy(copy=False)
        self.stops = self.data_frame.Stop.to_numpy(copy=False)
        self.lengths = self.data_frame.Length.to_numpy(copy=False)
        self.families = self.data_frame.Family.to_numpy(copy=False)
        self.sub_families = self.data_frame.SubFamily.to_numpy(copy=False)

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
