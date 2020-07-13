#!/usr/bin/env python3

"""
Wrappers for input data, the transposable elements.

Used to provide a common interface and fast calculations with numpy.
"""

__author__ = "Michael Teresi, Scott Teresi"

import logging
import numpy as np
import pandas as pd


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
            transposable_elements (pandas.DataFrame): transposable element data frame.
            genome_id (str): a string of the genome name, provided via arg
            parser in density.py
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
        self.genome_id = None

    @classmethod
    def mock(cls,
             start_stop=np.array([[0, 9], [10, 19], [20, 29]]), chromosome='Chr_ID',
             genome_id="fake_genome_id"):
        """Mocked data for testing.

        Args:
            start_stop (numpy.array): N gene x (start_idx, stop_idx)
        """

        n_genes = start_stop.shape[0]
        data = []
        family = "Family_0"  # FUTURE may want to parametrize family name later
        # NB overall order is not important but the names are
        columns = ['Start', 'Stop', 'Length', 'Order', 'SuperFamily', 'Chromosome']
        for gi in range(n_genes):
            g0 = start_stop[gi, 0]
            g1 = start_stop[gi, 1]
            gL = g1 - g0 + 1
            # FUTURE may want to parametrize sub
            # family name later
            subfam_suffix = "A" if gi % 2 else "B"
            subfamily = "SubFamily_{}".format(subfam_suffix)
            datum = [g0, g1, gL, family, subfamily, chromosome]
            data.append(datum)
        frame = pd.DataFrame(data, columns=columns)
        return TransposonData(frame, genome_id)


    def write(self, filename, key='default'):
        """Write a Pandaframe to disk.

        Args:
            filename (str): a string of the filename to write.
            key (str): identifier for the group (dataset) in the hdf5 obj.
        """
        self.data_frame.to_hdf(filename, key=key, mode='w')

    @classmethod
    def read(cls, filename, key='default'):
        """Read from disk. Returns a wrapped Pandaframe from an hdf5 file

        Args:
            filename (str): a string of the filename to write.
            genome_id (str): identifier for the genome name.
            key (str): identifier for the group (dataset) in the hdf5 obj.
        """
        new_instance = cls(pd.read_hdf(filename, key=key))
        new_instance.add_genome_id(genome_id)
        return new_instance

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
        """Printable representation for developers."""

        info = "TransposonData({self.data_frame})"
        return info.format(self=self)

    def add_genome_id(self, genome_id):
        """Add the genome_id as an extra column to the gene_dataframe"""
        self.data_frame.loc[:,'Genome_ID'] = genome_id
        self.genome_id = genome_id


    @property
    def genome_unique_id(self):
        """Return the unique genome identifier for all the TEs.

        Returns:
            str: the unique identifier.
        Raises:
            RuntimeError: if multiple TEs have different genome_id labels
            (not unique), I do not know how this could occur but want to be
            safe.
        """

        genome_id_list = self.data_frame.Genome_ID.unique().tolist()
        if not genome_id_list:
            raise RuntimeError("column 'Geneome_ID' is empty")
        elif len(genome_id_list) > 1:
            raise RuntimeError("Genome IDs are are not unique: %s" %
                               genome_id_list)
        else:
            return genome_id_list[0]  # MAGIC NUMBER list to string

    @property
    def chromosome_unique_id(self):
        """Unique chromosome identifier for all the genes available.

        This will raise f the genes are not from the same chromosome,
        for example you you didn't split the dataset wrt this data.

        Returns:
            str: the unique identifier.
        Raises:
            RuntimeError: if multiple chromosomes are in the data frame (i.e. no unique).
        """

        chromosome_list = self.data_frame.Chromosome.unique().tolist()
        if not chromosome_list:
            raise RuntimeError("column 'Chromosome' is empty")
        elif len(chromosome_list) > 1:
            raise RuntimeError("chromosomes are not unique: %s" % chromosome_list)
        else:
            return chromosome_list[0]  # MAGIC NUMBER list to string

    @property
    def superfamily_name_set(self):
        """The set of superfamily keys.

        Returns:
            set(str): superfamily grouping.
        """

        superfamilies = pd.unique(self.superfamilies)
        return set(superfamilies)

    @property
    def order_name_set(self):
        """The set of order keys.

        Returns:
            set(str): order grouping.
        """

        orders = pd.unique(self.orders)
        return set(orders)
