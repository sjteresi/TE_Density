#!/usr/bin/env python

"""
Contains transposon density data.
    Read and write to files, access the values.

DESIGN:
This is a refactor of MergeData.
The design goals are:
    - [x] update string storage for latest h5py
    - [ ] decouple density calculation from storage
    - [ ] support for 'reducing' multiple sources of density calculations
        - [x] support modification of existing files;
              support read / write, check if set exists and create it or,
              check if it's the wrong file
        - [ ] ingest and insert density results
        - [ ] add a 'bitmap' field designating which entries are valid
        - [ ] support start / restart (iterate over missing entries)

Other goals
    - [ ] add 'strand' column to gene_data and/or entry gene_datum
    - [ ] add 'strand' column to new density data (boolean)
    - [ ] add  'flipped' field for when antisense genes (-) have / have-not been flipped
        (all antisense genes must be swapped (left/right values) after calculation)
"""


__author__ = "Michael Teresi"


import logging
from dataclasses import dataclass, field
from functools import cache, partial

import numpy as np
import h5py


_STR_DTYPE = h5py.string_dtype("utf-8")  # MAGIC h5 convention for strings


@dataclass
class _DensitySubsetConfig:
    """Contains configuration parameters for reading / writing density data.

    Args:
        compression: HDF5 compression type: gzip | lzf | szip,
                     lzf recommended, see HDF5 'filter pipeline'
        windows: window sizes in number of base pairs
        te_names: the TE categories, names to either the superfamily or order
        gene_names: the names for each gene
    """

    compression: str = "gzip"
    windows: list[int] = field(default_factory=list)
    te_names: list[str] = field(default_factory=list)
    gene_names: list[str] = field(default_factory=list)


class _DensitySubset:
    """Contains the arrays for either a superfamily XOR order.

    This is an abstraction b/c both the superfamily/order have the same format.
    As such this is intended to be reused by the caller to handle the superfamily, order,
        or other category associated with Transposable Element density.

    Stores data in an HDF5 file, using an HDF5 'group'.
    Data is written to the group given the user ID provided for the 'prefix'.
    Data is initialized to the various values provided in the configuration.
    If the file does not contain the data, the fields will be initialized.
    If the file already contains the data,
        the fields will be read and compared to the configuration,
        and an exception will be thrown if the file data mismatches the configuration.

    NB this structure is *not* intended to manage the HDF5 file, rather it depends
        on an already open file and left to the caller to manage the resource on disk.
    """

    _LEFT = "_LEFT"
    _INTRA = "_INTRA"
    _RIGHT = "_RIGHT"

    _BITMAP = "_BITMAP"

    _GENES = "GENE_NAMES"
    _TRANSPOSONS = "TE_NAMES"  # NOTE this is the superfam | order, *needs confirmation*
    _WINDOWS = "WINDOWS"

    def __init__(self, file, prefix, config, logger=None):
        """Initialize.

        Args:
            file(h5py.File): an open HDF5 file
            prefix(str): user identifier, e.g. one of: {superfamily, order}
            config(_DensitySubsetConfig): options for dataset, e.g. compression type
            logger(logging.Logger): logger instance
        """

        self._logger = logger or logging.getLogger(self.__class__.__name__)
        self.file = file
        self._prefix = prefix
        self.cfg = config
        self._group = file.require_group(prefix)  # our convention, e.g. /superfamily

        self._init_gene_names()
        self._init_te_names()
        self._init_windows()

        self._init_densities()
        self.left = self._group[self._LEFT]
        self.intra = self._group[self._INTRA]
        self.right = self._group[self._RIGHT]
        self.bitmap = self._group[self._BITMAP]

    @property
    def filename(self):
        """Filename of the contained H5 data.

        Returns:
            str: path of the file
        """

        return self.file.filename

    @property
    def transposon_names(self):
        """The name of each transposon.

        Returns:
            generator(str): string format of each transposable element
        """

        # NB by h5py convention, str are stored as utf-8 (variable length)
        # NB using a generator b/c we expect a large number of TEs
        names_encoded = self._group[self._TRANSPOSONS][:]
        return (name.decode("utf-8") for name in names_encoded)

    # NB there should be < 10 instances of self per chromosome
    @property
    @cache  # pylint: disable=W1518
    def n_transposons(self):
        """Total number of transposable elements, e.g. superfamily."""

        return sum(1 for _ in self.transposon_names)

    @property
    def windows(self):
        """The window values.

        NB this could be an instance variable but is a property b/c
            these should not be changed once created.

        Returns:
            list(int): window sizes in number of base pairs
        """

        return self._group[self._WINDOWS][:]

    # NB there should be < 10 instances of self per chromosome
    @property
    @cache  # pylint: disable=W1518
    def n_windows(self):
        """Total number of windows."""

        return len(self.windows)

    @property
    def gene_names(self):  # NOTE should be gene_names?
        """The name of each gene.

        Returns:
            generator(str): string format of each gene
        """

        # NB by h5py convention, str are stored as utf-8 (variable length)
        # NB using a generator b/c we expect a large number of genes
        names_encoded = self._group[self._GENES][:]
        return (name.decode("utf-8") for name in names_encoded)

    # NB there should be < 10 instances of self per chromosome
    @property
    @cache  # pylint: disable=W1518
    def n_genes(self):
        """Total number of genes."""

        return sum(1 for _ in self.gene_names)

    def reduce(self, density_job):
        """Assign the left, intra, right densities and update the write bitmap."""

        raise NotImplementedError()

    def _read_dataset(self, key, shape, datatype, data=None, exact=False):
        """Read or reserve the dataset within the HDF5.

        If the data does NOT exist, create an empty dataset.
        Else, check the type and shape of the existing dataset.
        If the type and shape is compatible, load the data.

        Args:
            key(str): name of the dataset
            shape(numpy.shape): number or shape of elements
            datatype: numpy data type
        Raises:
            TypeError: dataset size requested to does not match existing dataset in file
        """

        comp = self.cfg.compression
        try:
            # NB this will read the data OR initialize if it DNE
            # SEE h5py.Group.create_dataset for other kwargs
            self._group.require_dataset(
                key, shape, datatype, data=data, exact=exact, compression=comp
            )
        except TypeError as err:
            msg = "input dataset for %s shape %s but the file is %s"
            self._logger.error(msg, key, shape, self._group[key].shape)
            raise err

    def _init_strings(self, key, data_from_file, data_from_config):
        """Initialize a string dataset in the file or validate the values if it exists.

        Args:
            key(str): name of the dataset
            data_from_file(list(string)): data read from the file
            data_from_config(list(string)): data from the user
        Raises:
            ValueError: data mismatch, cannot use existing file
        """

        # NB by h5py convention, the elements will be "" if group didn't exist for string
        data_from_file = list(data_from_file)
        data_from_config = list(data_from_config)
        is_empty = any(i for i in map(lambda i: i == "", data_from_file))
        n_elements = len(data_from_config)
        if is_empty:
            msg = "initialize %i elements for %s in file %s"
            self._logger.debug(msg, n_elements, key, self.filename)
            self._group[key][:] = data_from_config
            return

        msg = "read %i elements for %s in file %s"
        self._logger.debug(msg, n_elements, key, self.filename)
        if data_from_file != data_from_config:
            msg = (
                "input data for key '{}' does not match file: input v.s. file:\n{}\n{}"
            )
            self._logger.error(msg.format(key, data_from_config, data_from_file))
            raise ValueError(msg.format(key, data_from_config, data_from_file))

    def _init_array(self, key, data_in, dtype):
        """Initialize an array of numbers, reading if exists.

        Args:
            key(str):
            data_in(numpy.array): use the data to initalize if DNE
        Raises:
            ValueError: input data does not match the data in the file if it existed
        """

        self._read_dataset(key, len(data_in), dtype, exact=True, data=data_in)
        data_out = self._group[key][:]
        if data_in is not None and (data_in != data_out).any():
            msg = "input {} does not match file, \n{}\n{}"
            self._logger.error(msg.format(key, data_in, data_out))
            raise ValueError(msg.format(key, data_in, data_out))

    def _init_gene_names(self):
        """Write or read the gene names.

        NB this could be an instance variable but is a property b/c
            these should not be changed once created.

        Returns:
            list(int): window sizes in number of base pairs
        """

        key = self._GENES
        genes_in = self.cfg.gene_names
        self._read_dataset(key, len(genes_in), _STR_DTYPE)

        genes_out = self.gene_names
        self._init_strings(key, genes_out, genes_in)

    def _init_te_names(self):
        """Write or read the transposable element names, e.g. superfamily."""

        key = self._TRANSPOSONS
        tes_in = self.cfg.te_names
        self._read_dataset(key, len(tes_in), _STR_DTYPE)

        tes_out = self.transposon_names
        self._init_strings(key, tes_out, tes_in)

    def _init_windows(self):
        """Write or read the window values."""

        key = self._WINDOWS
        data_in = self.cfg.windows
        # MAGIC arbitrary integer, window matrix is small
        self._init_array(key, data_in, np.uint)

    def _init_densities(self):
        """Initialize the density containers and it's bitmap."""

        n_tes = self.n_transposons
        n_win = self.n_windows
        n_gen = self.n_genes
        # MAGIC intra does not use a window parameter so n_win is 1
        comp = self.cfg.compression

        # MAGIC experimental, so far needed more than f16 for density
        require = partial(self._group.require_dataset, compression=comp)
        require(self._LEFT, (n_tes, n_win, n_gen), np.double, exact=True)
        require(self._INTRA, (n_tes, 1, n_gen), np.double, exact=True)
        require(self._RIGHT, (n_tes, n_win, n_gen), np.double, exact=True)
        require(self._BITMAP, (n_tes, n_win, n_gen), bool, exact=True)
