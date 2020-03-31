#!/usr/bin/env python

"""
Combines overlap data with respect to the genes / transposons.
"""

__author__ = "Michael Teresi"

from collections import namedtuple
import logging
import os
from functools import partial

import h5py
import numpy as np

import transposon


_MergeConfigSink = namedtuple('_MergeConfigSink',
                             ['transposons', 'gene_names', 'windows', 'filepath', 'ram_bytes'])
_MergeConfigSource = namedtuple('_MergeConfigSource',
                               ['filepath'])
_Density = namedtuple('_Density', ['left', 'intra', 'right'])

""" Scott
please implement the items marked TODO SCOTT in this file
please add tests for said items
please implement the items marked TODO SCOTT in transposon_data.py & gene_data.py
"""

class _SumOverlaps():
    """Adds the overlaps with respect to each order || superfamily entry."""

    def __init__(self, overlap, transposon_group, transposon_subset, density, logger):
        """

        Args:
            overlap(OverlapData): active overlap file
            transposon_group(numpy.array): superfamily || order array of TransposonData
            transposon_subset(list(str)): superfamily || order set of identifiers
            density(_Density): container to output sums
        """

        self._overlap = overlap
        self._te_group = transposon_group
        self._te_subset = transposon_subset
        self.density = density
        self._logger = logger

    def _zip_input_outputs(self, windows):

        output_arrays = [
            self.density.left,
            self.density.intra,
            self.density.right
        ]
        input_arrays = [
            self._overlap.left,
            self._overlap.intra,
            self._overlap.right
        ]
        window_list = [
            windows,
            [1],  # TODO need gene len? or something else?
            windows
        ]
        where = [  # compare the ID to the array of names
            *(partial(numpy.equal, self._te_group),)*3,
        ]

    def sum(self, windows, genes, progress_bar):

        pass
        # for input, output, win, where in self._zip_input_outputs(windows):
        #   for te_group_id in self._te_subset
        #       for w_ in win  # track progress here? max of n_te_subset * windows?
        #           for gene in genes
        #               numpy.sum using where=where(te_group_id) out=output


class MergeData():
    """Contains the sum and density of the overlap wrt superfamilies and orders.

    Merges from one or more OverlapData into two sets: superfamilies, orders.
    The result has the following dimensions (not necessarily in this order):
        superfamilies x windows x genes x 3
        orders        x windows x genes x 3
    Where the three is for the left, intra, and right densities.

    Usage:
        Use the factory methods to create an instance (`from_param` | `from_file`).
        Use the context manager to activate the buffer.
        Use the slice methods to provide slices for the array public instance variables.
        Use the sum
    """

    DTYPE = np.float32  # MAGIC NUMBER experimental, depends on your data
    _GENE_NAMES = 'GENE_NAMES'
    _WINDOWS = 'WINDOWS'
    _CHROME_ID = 'CHROMOSOME_ID'
    _SUPERFAMILY_NAMES = 'SUPERFAMILY_NAMES'
    _ORDER_NAMES = 'ORDER_NAMES'
    _RHO_SUPERFAMILY_DATASET_KEY = 'RHO_SUPERFAMILIES'
    _S_LEFT = _RHO_SUPERFAMILY_DATASET_KEY + '_LEFT'
    _S_RIGHT = _RHO_SUPERFAMILY_DATASET_KEY + '_RIGHT'
    _S_INTRA = _RHO_SUPERFAMILY_DATASET_KEY + '_INTRA'
    _RHO_ORDER_DATASET_KEY = 'RHO_ORDERS'
    _O_LEFT = _RHO_ORDER_DATASET_KEY + '_LEFT'
    _O_RIGHT= _RHO_ORDER_DATASET_KEY + '_RIGHT'
    _O_INTRA = _RHO_ORDER_DATASET_KEY + '_INTRA'

    def __init__(self, configuration, compression='lzf', logger=None):

        self._logger = logger or logging.getLogger(__name__)
        self._config = configuration
        self._compression = self._check_compression(compression, self._logger)
        self._h5_file = None

        self.chromosome_id = None
        self.windows = None
        self.gene_names = None
        self.superfamily_names = None
        self.order_names = None
        self._window_2_idx = None
        self._gene_2_idx = None
        self._superfam_2_idx = None
        self._order_2_idx = None

        self.superfamily = None
        self.order = None

    @staticmethod
    def _check_compression(filter, logger):
        """Raise if compression type is invalid."""

        LOSSLESS_FILTERS = ['gzip', 'lzf', 'szip']  # SEE h5py datasets
        if filter not in LOSSLESS_FILTERS:
            msg = "'%s' compression invalid, options are: %s" % (filter, LOSSLESS_FILTERS)
            logger.critical(msg)
            raise ValueError(msg)
        else:
            return filter

    def __add__(self, other):
        """Combine density files."""
        raise NotImplementedError()

    @property
    def filepath(self):
        """Return the active filepath, None if no open file."""

        return self._h5_file.filename if self._h5_file is not None else None

    @classmethod
    def from_param(cls, transposon_data, gene_names, windows, output_dir, ram=2, logger=None):
        """Writable sink for a new file.

        Args:

        """

        logger = logger or logging.getLogger(__name__)
        ram_bytes= int(ram * 1024. ** 3)  # MAGIC NUMBER bytes to gigabytes
        transposon.check_ram(ram_bytes, logger)
        chrome = str(transposon_data.chromosome_unique_id)
        genome = str(transposon_data.genome_id)
        filename = genome + '_' + chrome + '.h5'
        filepath = os.path.join(output_dir, filename)
        config = _MergeConfigSink(
            transposons = transposon_data,
            gene_names = gene_names,
            windows = windows,
            filepath = filepath,
            ram_bytes = ram_bytes
        )
        return cls(config, logger=logger)

    @classmethod
    def from_file(cls, filepath, logger=None):
        """Read only source for an existing file."""

        raise NotImplementedError()

    def __enter__(self):
        """Context manager begin."""

        self._open_dispatcher()
        return self

    def __exit__(self, exc_type, exc_val, exc_traceback):
        """Context manager end."""

        self._h5_file.flush()
        self._h5_file.close()
        self._h5_file = None

        self.chromosome_id = None
        self.windows = None
        self.gene_names = None
        self.superfamily_names = None
        self.order_names = None
        self._window_2_idx = None
        self._gene_2_idx = None
        self._superfam_2_idx = None
        self._order_2_idx = None

        self.superfamily = None
        self.order = None


    def _open_dispatcher(self):
        """Open the file.

        NOTE when upgrading to python 3.8 refactor to use functools.singledispatchmethod.
        """

        if isinstance(self._config, _MergeConfigSink):
            self._open_new_file(self._config)
        elif isinstance(self._config, _MergeConfigSource):
            self._open_existing_file(self._config)
        else:
            raise TypeError("expecting {} or {} but got {}".format(
                type(_MergeConfigSink), type(_MergeConfigSource), type(self._config)))

    def _open_new_file(self, cfg):
        """Opens a new H5 file for writing and initializes the data sets.

        Args:
            cfg(_MergeConfigSink): configuration parameters
        """

        self.windows = list(cfg.windows)
        self.gene_names = list(cfg.gene_names)
        self.chromsome_id = cfg.transposons.chromosome_unique_id
        self.superfamily_names = cfg.transposons.superfamily_set
        self.order_names = cfg.transposons.order_set
        self._h5_file = h5py.File(cfg.filepath, 'w', rdcc_nbytes=cfg.ram_bytes)
        self._create_sets(self._h5_file, cfg)
        transposon.write_vlen_str_h5py(self._h5_file, self.windows, self._WINDOWS)
        transposon.write_vlen_str_h5py(self._h5_file, self.gene_names, self._GENE_NAMES)
        transposon.write_vlen_str_h5py(self._h5_file, self.chromsome_id, self._CHROME_ID)
        self._window_2_idx = {w: i for i, w in enumerate(self.windows)}
        self._gene_2_idx = {g: i for i, g in enumerate(self.gene_names)}


    def _open_existing_file(self, cfg):

        self._h5_file = h5py.File(cfg.filepath, 'r')
        transposon.read_vlen_str_h5py(self._h5_file, )
        # read windows
        # read chromosome id
        # read gene names
        # read superfams
        # read orders

    def _create_sets(self, h5_file, cfg):
        """Add data sets to the file.

        Args:
            h5_file(h5py.File): open h5 file
            cfg(_MergeConfigSink): tuple of configuration parameters
        """

        # FUTURE resizable / unlimited datasets are possible, consider for streaming
        # OverlapData into the merge? (since it can be split across many OverlapData)
        # then one does not need to know at this point how many genes / windows there are
        # that could be more complicated but more robust...
        create_set = partial(h5_file.create_dataset,
                             dtype=self.DTYPE,
                             compression=self._compression)

        n_win = len(self.windows)
        n_genes = sum(1 for n in self.gene_names)
        lr_shape = (n_win, n_genes)

        superfamilies = cfg.transposons.superfamily_set
        n_sfams = len(superfamilies)
        s_left = create_set(self._S_LEFT, (n_sfams, *lr_shape))
        s_intra = create_set(self._S_INTRA, (n_sfams, *lr_shape))
        s_right = create_set(self._S_RIGHT, (n_sfams, *lr_shape))
        self.superfamily = _Density(left=s_left, intra=s_intra, right=s_right)

        orders = cfg.transposons.order_set
        n_orders = len(orders)
        o_left = create_set(self._O_LEFT, (n_orders, *lr_shape))
        o_intra = create_set(self._O_INTRA, (n_orders, *lr_shape))
        o_right = create_set(self._O_RIGHT, (n_orders, *lr_shape))
        self.order = _Density(left=o_left, intra=o_intra, right=o_right)

    def _sum_superfam(self, overlap):
        raise NotImplementedError()

    def _sum_order(self, overlap):
        raise NotImplementedError()

    def sum(self, overlap, progress_bar):
        """Sum across the superfamily / order dimension."""

        self._validate_chromosome(overlap)
        self._validate_windows(overlap)
        self._validate_gene_names(overlap)

        o_windows = overlap.windows
        o_genes = overlap.gene_names

        sfam_set = self._config.transposons.superfamily_set
        n_sfam = len(sfam)


        order_set = self._config.transposons.order_set

        # for in, out, windows, genes
        #

    def _zip_density_inputs(subset_names, transposons, te_subset, density):
        """Zip the names for calculating density.

        Args:
            subset_names(list(str)): list of superfamily or order identifiers
            transposons(TransposonData): the transposon container
            te_subset(numpy.array): superfamily or order column of transposon container
            density(_Density): density output container
        Returns:

        """

        output_arrays = [
            density.left,
            density.intra,
            density.right
        ]
        input_arrays = [
            overlap.left,
            overlap.intra,
            overlap.right
        ]
        where = [
            *(partial(numpy.equal, te_subset),)*3,  #
        ]
        return zip(subset_names, output_arrays, input_arrays, where)



    def _validate_chromosome(self, overlap):
        """ValueError if the chromosome ID of the overlap data does not match this.

        The OverlapData cannot be marged into a container created for another chromosome.
        This is because the genes may be different, and more importantly the density is
            not defined between chromosomes as they are different locations.
        """

        pass  # TODO SCOTT

    def _validate_windows(self, overlap):
        """ValueError if the overlap windows are not in the destination."""

        pass  # TODO SCOTT

    def _validate_gene_names(self, overlap):
        """ValueError if the overlap genes are not in the destination."""

        pass  # TODO SCOTT







class MergeWorker():
    """Orchestrates merging multiple OverlapData."""
    pass


    # scrape the windows from the overlap files, pass to MergeData
    # scrape genes from the overlap files, pass to MergeData
