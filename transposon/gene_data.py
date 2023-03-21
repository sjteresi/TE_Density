#!/usr/bin/env python3

"""
Wrappers for input data, multiple genes.

Used to provide a common interface and fast calculations with numpy.
"""

__author__ = "Michael Teresi, Scott Teresi"

import logging
import numpy as np
import pandas as pd

from transposon.gene_datum import GeneDatum


class GeneData(object):
    """Wraps a data frame containing many genes.
    Provides an interface, attribute access, and to/from disk functionality.

    Note the numpy views are not necessarily a no-copy (SEE pandas.DataFrame.to_numpy).

    Expects certain column identifiers (SEE self.__init__).
    GeneData subclasses should conform to these column names or redefine the properties.
    """

    def __init__(self, gene_dataframe, genome_id, logger=None):
        """Initialize.

        Args:
            gene_dataframe (DataFrame): gene data frame.
            genome_id (str): a string of the genome name.
        """

        self._logger = logger or logging.getLogger(__name__)
        self.data_frame = gene_dataframe.copy(deep=True)
        self._names = self.data_frame.index  # names of genes list(str)
        self.starts = self.data_frame.Start.to_numpy(copy=False)
        self.stops = self.data_frame.Stop.to_numpy(copy=False)
        self.lengths = self.data_frame.Length.to_numpy(copy=False)
        self.chromosomes = self.data_frame.Chromosome.to_numpy(copy=False)
        self.genome_id = genome_id
        self.add_genome_id()

    @classmethod
    def mock(
        cls,
        start_stop=np.array([[0, 9], [10, 19], [20, 29]]),
        genome_id="fake_genome_id",
    ):
        """Mocked data for testing.

        Args:
            start_stop (numpy.array): N gene x (start_idx, stop_idx)
        """

        n_genes = start_stop.shape[0]
        data = []
        for gi in range(n_genes):
            g0 = start_stop[gi, 0]
            g1 = start_stop[gi, 1]
            gL = g1 - g0 + 1
            name = "gene_{}".format(gi)
            chromosome = genome_id
            datum = [name, g0, g1, gL, chromosome]
            data.append(datum)

        col_names = ["Gene_Name", "Start", "Stop", "Length", "Chromosome"]
        frame = pd.DataFrame(data, columns=col_names)
        frame.set_index("Gene_Name", inplace=True)
        return GeneData(frame, genome_id)

    @classmethod
    def mock_v2(
        cls,
        start_stop=np.array([[500, 1000], [800, 1500], [1600, 2000]]),
        chromosome_ids=["Chrom_1", "Chrom_2", "Chrom_3"],
        genome_id="fake_genome_id",
    ):
        """Mocked data for testing. TODO refactor with main mock.

        Args:
            start_stop (numpy.array): N gene x (start_idx, stop_idx)
            chromosome_id (list): List of string, len(N gene)
        """

        n_genes = start_stop.shape[0]
        data = []
        for gi in range(n_genes):
            g0 = start_stop[gi, 0]
            g1 = start_stop[gi, 1]
            gL = g1 - g0 + 1
            name = "gene_{}".format(gi)
            chromosome_id = chromosome_ids[gi]
            datum = [name, g0, g1, gL, chromosome_id]
            data.append(datum)

        col_names = ["Gene_Name", "Start", "Stop", "Length", "Chromosome"]
        frame = pd.DataFrame(data, columns=col_names)
        frame.set_index("Gene_Name", inplace=True)
        return GeneData(frame, genome_id)

    def write(self, filename):
        """Write a Pandaframe to disk.

        Args:
            filename (str): a string of the filename to write.

        """
        # Begin refactor
        self.data_frame.to_csv(filename, sep="\t", header=True, index=True)

    @classmethod
    def read(cls, filename):
        """Read from disk. Returns a wrapped Pandaframe from an hdf5 file

        Args:
            filename (str): a string of the filename to write.
        """
        # NOTE this is a little duplicate with import_filtered_genes
        # NOTE I don't have logger obj here, so talk to Mike if we just want to have
        # duplicate code.
        data_frame = pd.read_csv(
            filename,
            header="infer",
            sep="\t",
            index_col="Gene_Name",
            dtype={
                "Start": "float64",
                "Stop": "float64",
                "Length": "float64",
                "Chromosome": str,
                "Strand": str,
                "Feature": str,
                "Genome_ID": str,
            },
        )
        data_frame.sort_values(by=["Chromosome", "Start"], inplace=True)
        genome_id_list = data_frame["Genome_ID"].unique().tolist()
        if not genome_id_list:
            raise RuntimeError("column 'Genome_ID' is empty")
        elif len(genome_id_list) > 1:
            raise RuntimeError("Genome IDs are are not unique: %s" % genome_id_list)
        else:
            genome_id = genome_id_list[0]  # MAGIC NUMBER list to string

        data_frame.drop(columns=["Genome_ID"], inplace=True)  # NOTE, have to do
        # this to avoid it getting the column twice, perhaps refactor?

        new_instance = cls(data_frame, genome_id)
        return new_instance

    def get_gene(self, gene_id):
        """Return a GeneDatum for the gene identifier."""

        return GeneDatum(self.data_frame, gene_id)

    def add_genome_id(self):
        """Add the genome_id as an extra column to the gene_dataframe"""
        self.data_frame["Genome_ID"] = self.genome_id

    @property
    def names(self):
        """Yields the names for each gene."""

        return (name for name in self._names)

    @property
    def chromosome_unique_id(self):
        """Unique chromosome identifier for all the genes available.

        This will raise if the genes are not from the same chromosome,
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

    def __repr__(self):
        """String representation for developer."""

        return "GeneData{}".format(self.data_frame)
