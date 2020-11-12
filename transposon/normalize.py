#!/usr/bin/env python3

__author__ = "Scott Teresi"

"""
Output the relevant area (divisor) for density calculation.
"""

from collections import namedtuple
import logging
from enum import Enum, unique
import os
from functools import partial

import numpy as np
import pandas as pd

from transposon.density import validate_window

from transposon.merge_data import MergeData
from transposon import MAX_SYSTEM_RAM_GB, check_ram

# NB SCOTT
# we need:
# NormalizedData: class for storing the normalized results (e.g. same as OverlapData / MergeData)
#                 stores / provides access to the normalized data
# NormalizeWorker: class to do the actual calculation (or you can add this to Normalize) (similar to OverlapWorker)
#                  performs calculations to create normalized datum (which is then stored in above)
#                  Use the slicing notation that Michael has to interact with the specific window you need.
#                  Normalization matrix is in shape (windows, genes), Michael's slice matrix will have to be maniuplated to interact with this
#                  Emulate Michael's class factories (named tuples) at the top of Overlap.py
# NormalizeManager: class to do the multiprocess stuff, the same as the OverlapManager
#                   multiprocessing to execute the jobs


# 11/10/2020:
# Go over the data structures, make sure NormalizeData makes sense
# Discuss with Michael best way to store the data in H5 format. Maybe get a
# small implementation of actual slicing and reformat the normalize data to do
# math
# Tie everything together.
# Normalize has the appropriate functions, would it be easier to refactor those
# functions so that it takes in a GeneDatum (single gene) rather than the
# entire annotation? That was it could work more like the overlap calculations.


# 2020/11/10
# feedback 
# move NormalizeData to normalize_data.py
# move NormalizeWorker to normalize_worker.py

_NormalizeConfigSink = namedtuple(
    "_NormalizeConfigIn", ["genes", "te_data", "windows", "filepath", "ram_bytes"]
)
_NormalizeConfigSource = namedtuple("_NormalizeConfigSource", ["filepath"])

NormalizeResult = namedtuple("NormalizeResult", ["genes_processed", "filepath", "exception", "gene_file"])

NormalizeResult.__new__.__defaults__ = (
    0,
    None,
    None,
    "",
)  # REFACTOR to dataclass in 3.7+




class NormalizedData():
    """Contains normalized density values of genes wrt window, direction,
    superfamilies and orders.


    INCOMPLETE:


    Contains a file resource which must be acquired before writing | reading.
    Use the context manager or start/stop methods to manage the file resource.


    SCOTT:
        Needs to be able to open and interact with h5 files and write them.

    Usage:
        Use the factory methods to create an instance (`from_param` | `from_file`).
        Use the context manager to activate the buffer.
        Use the slice methods to provide slices for the array public instance variables.
    """

    # NOTES:
        # NormalizeData should be similar in structure to OverlapData, just
        # with different values. The main things that need to change between
        # the two classes are the set creations methods 

    DTYPE = np.float32  # MAGIC NUMBER experimental, depends on your data
    _LEFT = Normalize.Direction.LEFT.name
    _RIGHT = Normalize.Direction.RIGHT.name
    _INTRA = Normalize.Direction.INTRA.name
    _GENE_NAMES = "GENE_NAMES"
    _WINDOWS = "WINDOWS"
    _CHROME_ID = "CHROMOSOME_ID"
    _GENOME_ID = "GENOME_ID"

    def __init(self, configuration, compression="lzf", logger=None):
        """Initializer.

        Args:
            configuration (tuple): `_NormalizeConfigSink` || `_NormalizeConfigSource`.
            compression (str): hdf5 compression type
            logger (logging.Logger): logger instance, creates one if None.
        """
        self._logger = logger or logging.getLogger(__name__)
        self._config = configuration
        self._h5_file = None
        #self._compression = self._check_compression(compression, self._logger)

        self.chromosome_id = None


        # TODO should I also have something to accomodate the TE data?
        # TODO how do I accomodate the transposon data in the config? I guess
        # that is already done due to being in a named tuple and we set the
        # configuration above?

        self.left = None
        self.intra = None
        self.right = None
        self.gene_names = None
        self.chromosome_id = None
        self.genome_id = None
        self.windows = None

    @classmethod
    def from_file(cls, filepath, logger=None):
        """Read only source for an existing file."""

        file_abs = os.path.abspath(filepath)
        if not os.path.isfile(file_abs):
            raise ValueError("input filepath not a file: %s" % file_abs)
        config = _NormalizeConfigSource(filepath=file_abs)
        normalize = cls(config, logger)
        return normalize

    @classmethod
    def from_param(
        cls, genes, te_data, windows, output_dir, ram=1.2, logger=None
    ):
        """Writable sink for a new file.

        Args:
            genes (GeneData): genes container
            te_data (TransposonData): transposon container
            windows (list(int)): window sizes
            output_dir (str): directory for output files
            ram (int): upper limit for caching in gigabytes
        """
        # NOTE probably do not need n_transposons, just the types of
        # transposons, so something like a list of superfamilies and orders, or
        # if you give it the TransposonData for that specific chromosome we can
        # get it that way.

        logger = logger or logging.getLogger(__name__)
        ram_bytes = int(ram * 1024.0 ** 3)  # MAGIC NUMBER bytes to gigabytes
        check_ram(ram_bytes, logger)

        # NOTE, question, is there a way to collate all of these together?
        # Maybe that can be done after the fact but having an h5 file for each
        # chromosome isn't that bad, it just may be better to have all of the
        # chromosomes together in one massive h5 file.

        chromosome = genes.chromosome_unique_id
        filename = chromosome + "_density_" + ".h5"  # name files by chromosome
        filepath = os.path.join(output_dir, filename)

        config = _NormalizeConfigSink(
            genes=genes,
            te_data=te_data,
            windows=windows,
            filepath=filepath,
            ram_bytes=ram_bytes,
        )
        return cls(config, logger)


    @staticmethod
    def _check_compression(filter, logger):
        """Raise if compression type is invalid."""
        # NOTE Seems complete, brought over from MergeData, may not be needed

        LOSSLESS_FILTERS = ["gzip", "lzf", "szip"]  # SEE h5py datasets
        if filter not in LOSSLESS_FILTERS:
            msg = "'%s' compression invalid, options are: %s" % (
                filter,
                LOSSLESS_FILTERS,
            )
            logger.critical(msg)
            raise ValueError(msg)
        else:
            return filter

    @property
    def filepath(self):
        """Return the active filepath, None if no open file."""
        # NOTE Seems complete

        return self._h5_file.filename if self._h5_file is not None else None

    def start(self):
        """Obtain the resource by opening the file."""

        self._open_dispatcher()

    def stop(self):
        """Release resources by closing the file."""

        self._h5_file.flush()  # method inherent to h5 to write
        self._h5_file.close() # method inherent to h5 to close
        self._h5_file = None
        self._config = None

        self.left = None
        self.right = None
        self.intra = None
        self.gene_names = None
        self.chromosome_id = None
        self.windows = None

    def __enter__(self):
        """Context manager start."""

        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_traceback):
        """Context manager stop."""

        self.stop()

    def _create_sets(self, h5_file, cfg):
        """Initialize instance variables from input configuration, mutates self.

        Args:
            h5_file (hdf5.File): the open file.
            cfg (_NormalizeConfigSink): the configuration for a new file.
        """
        self.chromosome_id = cfg.genes.chromosome_unique_id
        self.genome_id = cfg.genes.genome_id
        create_set = partial(
            h5_file.create_dataset, dtype=self.DTYPE, compression=self.COMPRESSION
        )
        self.gene_names = list(cfg.genes.names)
        n_genes = len(self.gene_names)
        self.windows = list(cfg.windows)
        n_win = len(self.windows)
        n_tes = cfg.n_transposons  # TODO modify this usage

        # TODO check to see if the following is right

        # N.B. numpy uses row major by default, iterate over data accordingly
        # this is coupled with OverlapWorker.calculate which iterates
        # this is coupled with the slicing methods of self
        left_right_shape = (
            n_genes,
            n_win,
            n_tes,# TODO modify this usage
        )
        # TODO validate chunk dimensions, 4 * nwin * ntes may be too big pending inputs
        # (and probably already is, but it worked ok...)
        gene_chunks = min(32, n_genes)
        chunks = tuple((gene_chunks, n_win, n_tes))  # MAGIC NUMBER experimental
        # TODO modify this usage
        self.left = create_set(self._LEFT, left_right_shape, chunks=chunks)
        self.right = create_set(self._RIGHT, left_right_shape, chunks=chunks)
        intra_shape = (n_genes, n_tes)
        self.intra = create_set(self._INTRA, intra_shape)

    def _open_dispatcher(self):
        """Open the file.

        NOTE when upgrading to python 3.8 refactor to use functools.singledispatchmethod.
        """

        if isinstance(self._config, _NormalizeConfigSink):
            self._open_new_file(self._config)
        elif isinstance(self._config, _NormalizeConfigSource):
            self._open_existing_file(self._config)
        else:
            raise TypeError(
                "expecting {} or {} but got {}".format(
                    _NormalizeConfigSink, _NormalizeConfigSource, self._config
                )
            )

    def _open_existing_file(self, cfg):
        """Open the file, mutates self."""

        self._h5_file = h5py.File(cfg.filepath, "r")
        self.gene_names = self._read_gene_names()
        self.windows = self._read_windows()
        self.chromosome_id = self._read_chromosome_id()
        self.genome_id = self._read_genome_id()
        self.left = self._h5_file[self._LEFT]
        self.intra = self._h5_file[self._INTRA]
        self.right = self._h5_file[self._RIGHT]

    def _open_new_file(self, cfg):
        """Initialize a new file for writing."""

        self._h5_file = h5py.File(cfg.filepath, "w", rdcc_nbytes=cfg.ram_bytes)
        self._create_sets(self._h5_file, cfg)
        self._write_gene_names()
        self._write_windows()
        self._write_chromosome_id()
        self._write_genome_id()

    def _write_gene_names(self):
        """Assign list of gene names to the file."""
        raise NotImplementedError

    def _read_gene_names(self):
        """Return list of gene names."""

        return self._h5_file[self._GENE_NAMES][:].tolist()

    def _write_windows(self):
        """Assign list of windows to the file."""

        self._h5_file.create_dataset(
            self._WINDOWS, data=np.array(self.windows), dtype=np.int64
        )

    def _read_windows(self):
        """Return list of windows."""

        return self._h5_file[self._WINDOWS][:].tolist()

    def _write_chromosome_id(self):
        """Assign chromosome identifier to the file."""

        # FUTURE consider ASCII (fixed len) for storage if non-python clients exist
        vlen = h5py.special_dtype(vlen=str)
        # MAGIC NUMBER there can only be one unique ID
        dset = self._h5_file.create_dataset(self._CHROME_ID, (1,), dtype=vlen)
        dset[:] = self.chromosome_id

    def _read_chromosome_id(self):
        """Return the unique chromosome identifier."""

        return self._h5_file[self._CHROME_ID][:].tolist()[0]  # MAGIC NUMBER only one ID

    def _write_genome_id(self):
        """Assign genome identifier to file."""

        vlen = h5py.special_dtype(vlen=str)
        # MAGIC NUMBER there can only be one unique ID
        dset = self._h5_file.create_dataset(self._GENOME_ID, (1,), dtype=vlen)
        dset[:] = self.genome_id

    def _read_genome_id(self):
        """Read genome identifier."""

        return self._h5_file[self._GENOME_ID][:].tolist()[0]  # MAGIC NUMBER only one ID



class NormalizeWorker():
    """Calculates the density data.

    Divides the corresponding matrixes of MergeData and the normalized values
    """


    def __init__(self, normalize_dir, logger=None):
        """Initialize.

        Args:
            normalize_dir (str): A string path of the directory to output te
            normalized files. This comes from the ArgParser obj in density.py and
            defaults to tmp.
        """
        self.root_dir = overlap_dir
        self._logger = logger or logging.getLogger(__name__)

        self._data = None
        self._windows = None
        self._gene_names = None

        self._gene_name_2_idx = None
        self._window_2_idx = None

    def calculate(
        self, genes, windows, gene_names, progress=None
    ):
        """Calculate the normalization matrix for the genes and windows.

        N.B. the number of runs of this multi level loop is intended to be
        limited by the number of input genes / windows rather than refactoring
        into a multiprocess solution.

        Args:
            genes (GeneData): the genes; the chromosome with gene annotations.
            windows (list(int)): iterable of window sizes to use.
            gene_names (list(str)): iterable of the genes to use.
            progress (Callable): callback for when a gene name is processed.
        """
        # TODO the gene_names variable may not be needed
        # TODO flesh out args

        self._reset()  # TODO

        summed_overlaps = MergeData.from_file(filepath)  # mergedata for a
        # given chromosome, but is this where we should call from_file?
        # Probably not.

        # TODO Reshape the MergeData summed matrix into a shape amenable for
        # matrix math with shape (n windows X n genes) which is the divisors

        # For each window we only need ONE row of the divisor matrix. So if I
        # can just index for all groups, for all genes, at a specific window
        # index on the MergeData, I only need to then index the divisor frame
        # on the appropriate window index and division can happen because they
        # will share a same shape (genes).

        # This is where I do the invoke the slicing?








        for gene_name in genes.names:
            # TODO check stop event


    def divisors():
        raise NotImplementedError

        # NOTE, Michael do you think these ought to be attached to an instance
        # of the class? They only need to be calculated once per chromosome
        # i.e (group of unique genes). 

        # shape (n windows X n genes), requires GeneData and iterable of windows
        left_divisors = NormalizationDivisors.divisor_left(genes, windows)
        intra_divisors = NormalizationDivisors.divisor_intra(genes, windows)
        right_divisors = NormalizationDivisors.divisor_intra(genes, windows)

    def rho():
        raise NotImplementedError


    def _reset(self, transposons, genes, windows, gene_names):
        """Initialize overlap data; mutates self."""
        raise NotImplementedError
        # Something like?
        gene_names_filtered = self._filter_gene_names(gene_names, genes.names)
        self._gene_names = list(gene_names_filtered)
        self._gene_name_2_idx = self._map_gene_names_2_indices(self._gene_names)
        self._windows = list(self._filter_windows(windows))
        self._window_2_idx = self._map_windows_2_indices(self._windows)

        # TODO fill out
        self._data = NormalizedData.from_param()


    def _filter_gene_names(self, gene_names_requested, gene_names_available):
        """Yield only the valid gene names."""

        for name in gene_names_requested:
            # could use list comprehension but we need to log if it fails
            if name not in gene_names_available:
                self._logger.error(" invalid gene name: %s", name)
            else:
                yield name


class NormalizeManager():
    """Orchestrate multiple NormalizeWorkers to calculate normalized values
    (density)


    """
    def __init__(self):
        raise NotImplementedError

