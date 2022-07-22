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

    def __init__(self, transposable_elements_dataframe, genome_id, logger=None):
        """Initialize.

        Args:
            transposable_elements (pandas.DataFrame): transposable element data frame.
            genome_id (str): a string of the genome name, provided via arg
            parser in density.py
        """

        self._logger = logger or logging.getLogger(__name__)
        self.data_frame = transposable_elements_dataframe.copy()
        self.indices = self.data_frame.index.to_numpy(copy=False)
        self.starts = self.data_frame.Start.to_numpy(copy=False)
        self.stops = self.data_frame.Stop.to_numpy(copy=False)
        self.lengths = self.data_frame.Length.to_numpy(copy=False)
        self.orders = self.data_frame.Order.to_numpy(copy=False)
        self.superfamilies = self.data_frame.SuperFamily.to_numpy(copy=False)
        self.chromosomes = self.data_frame.Chromosome.to_numpy(copy=False)
        self.genome_id = genome_id
        self.add_genome_id()

    @classmethod
    def mock(
        cls,
        start_stop=np.array([[0, 9], [10, 19], [20, 29]]),
        chromosome="Chr_ID",
        genome_id="fake_genome_id",
    ):
        """Mocked data for testing.

        Args:
            start_stop (numpy.array): N gene x (start_idx, stop_idx)
        """

        n_genes = start_stop.shape[0]
        data = []
        family = "Family_0"  # FUTURE may want to parametrize family name later
        # NB overall order is not important but the names are
        columns = ["Start", "Stop", "Length", "Order", "SuperFamily", "Chromosome"]
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

    @classmethod
    def mock_v2(
        cls,
        start_stop=np.array([[0, 9], [10, 19], [20, 29]]),
        order_ids=["LTR", "TIR", "LINE"],
        superfamily_ids=["Gypsy", "Mutator", "L1"],
        chromosome_ids=["Chrom_1", "Chrom_2", "Chrom_3"],
        genome_id="fake_genome_id",
    ):
        """Mocked data for testing. TODO refactor with main mock.

        Args:
            start_stop (numpy.array): N gene x (start_idx, stop_idx)
        """
        # TODO check to make sure all the arrays are same length

        n_tes = start_stop.shape[0]
        data = []
        columns = ["Start", "Stop", "Length", "Order", "SuperFamily", "Chromosome"]
        for ti in range(n_tes):
            g0 = start_stop[ti, 0]
            g1 = start_stop[ti, 1]
            gL = g1 - g0 + 1
            chromosome_id = chromosome_ids[ti]
            order_id = order_ids[ti]
            superfamily_id = superfamily_ids[ti]
            datum = [g0, g1, gL, order_id, superfamily_id, chromosome_id]
            data.append(datum)
        frame = pd.DataFrame(data, columns=columns)
        return TransposonData(frame, genome_id)

    def write(self, filename, key="default"):
        """Write a Pandaframe to disk.

        Args:
            filename (str): a string of the filename to write.
            key (str): identifier for the group (dataset) in the hdf5 obj.
        """
        self.data_frame.to_hdf(filename, key=key, mode="w")

    @classmethod
    def read(cls, filename, key="default"):
        """Read from disk. Returns a wrapped Pandaframe from an hdf5 file

        Args:
            filename (str): a string of the filename to write.
            key (str): identifier for the group (dataset) in the hdf5 obj.
        """
        data_frame = pd.read_hdf(filename, key=key)
        genome_id_list = data_frame["Genome_ID"].unique().tolist()
        if not genome_id_list:
            raise RuntimeError("column 'Genome_ID' is empty")
        elif len(genome_id_list) > 1:
            raise RuntimeError("Genome IDs are are not unique: %s" % genome_id_list)
        else:
            genome_id = genome_id_list[0]  # MAGIC NUMBER list to string
        new_instance = cls(data_frame, genome_id)
        return new_instance

    def write_human_readable(self, filename):
        """
        Write the Pandaframe to disk as a tsv.

        Args:
            filename (str): a string of the filename to write.
        """

        self.data_frame.to_csv(filename, header=True, index=False, sep="\t")

    def write_gff(self, filename):
        """
        Write the Pandaframe to disk as a GFF.

        Args:
            filename (str): a string of the filename to write.
        """
        # TODO Come back to this post-v1.0 build and flesh out.
        self.gff_data_frame = self.data_frame.copy(deep=True)
        self.gff_data_frame.rename(columns={"Stop": "End"}, inplace=True)
        self.gff_data_frame["Source"] = "Test"
        self.gff_data_frame["Score"] = "."
        self.gff_data_frame["Phase"] = "1"
        self.gff_data_frame["Attributes"] = "Test"
        self.gff_data_frame["Feature"] = (
            self.gff_data_frame["Order"] + "/" + self.gff_data_frame["SuperFamily"]
        )

        self.gff_data_frame = self.gff_data_frame.astype(
            {"Start": int, "End": int, "Phase": int}
        )

        self.gff_data_frame = self.gff_data_frame[
            [
                "Chromosome",
                "Source",
                "Feature",
                "Start",
                "End",
                "Score",
                "Strand",
                "Phase",
                "Attributes",
            ]
        ]
        self.gff_data_frame = self.gff_data_frame.sort_values(by=["Start"])
        self.gff_data_frame.to_csv(filename, header=False, index=False, sep="\t")

    @property
    def number_elements(self):
        """The number of transposable elements."""
        return self.indices.shape[0]  # MAGIC NUMBER it's one column

    def subset_by_superfam(self):
        """
        Return a list of dataframes containing only one type of superfamily
        at a time
        """
        pass
        # for superfamily in self.superfamily_name_set:

    def check_shape(self):
        """Checks to make sure the columns of the TE data are the same size.

        If the shapes don't match then there are records that are incomplete,
            as in an entry (row) does have all the expected fields (column).
        """

        start = self.starts.shape
        stop = self.stops.shape
        if start != stop:
            msg = "Input TE missing fields: starts.shape {}  != stops.shape {}".format(
                start, stop
            )
            self._logger.critical(msg)
            raise ValueError(msg)

        # TODO verify that this isn't weird
        length = self.lengths.shape
        if start != length:
            msg = (
                "Input TE missing fields: starts.shape {}  != lengths.shape {}".format(
                    start, stop
                )
            )
            self._logger.critical(msg)
            raise ValueError(msg)

    def __repr__(self):
        """Printable representation for developers."""

        info = "TransposonData({self.data_frame})"
        return info.format(self=self)

    def add_genome_id(self):
        """Add the genome_id as an extra column to the gene_dataframe"""
        self.data_frame["Genome_ID"] = self.genome_id

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
