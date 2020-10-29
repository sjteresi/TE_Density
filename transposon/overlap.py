#!/usr/bin/env python3

__author__ = "Michael Teresi"

"""
The no. base pairs that the transposable elements overlap the window for a gene.
"""

from collections import namedtuple
import copy
from enum import Enum, unique
import logging
import os
from functools import partial
import tempfile

import numpy as np
import h5py

from transposon import MAX_SYSTEM_RAM_GB, check_ram


_OverlapConfigSink = namedtuple(
    "_OverlapConfigIn", ["genes", "n_transposons", "windows", "filepath", "ram_bytes"]
)
_OverlapConfigSource = namedtuple("_OverlapConfigSource", ["filepath"])
OverlapResult = namedtuple(
    "OverlapResult", ["genes_processed", "filepath", "exception", "gene_file"]
)
OverlapResult.__new__.__defaults__ = (
    0,
    None,
    None,
    "",
)  # REFACTOR to dataclass in 3.7+


class Overlap:
    """Functions for calculating overlap."""

    @unique
    class Direction(Enum):
        """Where the overlap is calculated relevant to the gene."""

        LEFT = 0
        INTRA = 1
        RIGHT = 2

    @staticmethod
    def left(gene_datum, transposons, window):
        """Overlap to the left of the gene.

        Number base pair overlap between the TE / gene on the relevant distance.

        Args:
            gene_datum (GeneDatam): the relevant gene.
            transposons (TransposonData): the transposable elements.
            window (int): no. base pairs from gene to calculate overlap.
        """

        w_start = gene_datum.left_win_start(window)
        w_stop = gene_datum.left_win_stop
        # TODO implement validate window to fix the window length for the
        # density caluclations, because w_len is not needed for the overlaps.
        # w_len = gene_datum.win_length(window)
        # w_len = validate_window(w_start, gene_datum.start, w_len)
        lower_bound = np.maximum(w_start, transposons.starts)
        upper_bound = np.minimum(w_stop, transposons.stops)
        te_overlaps = np.maximum(0, (upper_bound - lower_bound + 1))
        return te_overlaps

    @staticmethod
    def intra(gene_datum, transposons):
        """Overlap to the gene itself.

        Number base pair overlap between the TE / gene on the relevant distance.

        Args:
            gene_datum (GeneDatam): the relevant gene.
            transposons (TransposonData): the transposable elements.
            window (int): no. base pairs from gene to calculate overlap.
        """

        g_start = gene_datum.start
        g_stop = gene_datum.stop
        lower_bound = np.minimum(g_stop, transposons.stops)
        upper_bound = np.maximum(g_start, transposons.starts)
        te_overlaps = np.maximum(0, (lower_bound - upper_bound + 1))
        return te_overlaps

    @staticmethod
    def right(gene_datum, transposons, window):
        """Overlap to the right of the gene.

        Number base pair overlap between the TE / gene on the relevant distance.

        Args:
            gene_datum (GeneDatam): the relevant gene.
            transposons (TransposonData): the transposable elements.
            window (int): no. base pairs from gene to calculate overlap.
        """

        w_start = gene_datum.right_win_start
        w_stop = gene_datum.right_win_stop(window)
        lower_bound = np.maximum(w_start, transposons.starts)
        upper_bound = np.minimum(w_stop, transposons.stops)
        te_overlaps = np.maximum(0, (upper_bound - lower_bound + 1))
        return te_overlaps


class OverlapData:
    """Contains overlap values buffered to disk.

    Contains a file resource which must be acquired before writing | reading.
    Use the context manager or start/stop methods to manage the file resource.

    Usage:
        Use the factory methods to create an instance (`from_param` | `from_file`).
        Use the context manager to activate the buffer.
        Use the slice methods to provide slices for the array public instance variables.
    """

    DTYPE = np.float32  # MAGIC NUMBER experimental, depends on your data
    COMPRESSION = "lzf"  # MAGIC NUMBER experimental, fast w/ decent compression ratio
    _LEFT = Overlap.Direction.LEFT.name
    _RIGHT = Overlap.Direction.RIGHT.name
    _INTRA = Overlap.Direction.INTRA.name
    _GENE_NAMES = "GENE_NAMES"
    _WINDOWS = "WINDOWS"
    _CHROME_ID = "CHROMOSOME_ID"
    _GENOME_ID = "GENOME_ID"

    def __init__(self, configuration, logger=None):
        """Initializer.

        Args:
            configuration (tuple): `_OverlapConfigSink` || `_OverlapConfigSource`.
            logger (logging.Logger): logger instance, creates one if None.
        """

        self._logger = logger or logging.getLogger(__name__)
        self._config = configuration
        self._h5_file = None

        self.left = None
        self.intra = None
        self.right = None
        self.gene_names = None
        self.chromosome_id = None
        self.genome_id = None
        self.windows = None

        # TODO add public proteced variable for gene name to index
        # specify that the indices in the overlap matrices are the indices
        # in the name list

    def __str__(self):
        """User representation."""

        return "{}".format(self._config)

    @property
    def filepath(self):
        """Return the active filepath, None if no open file."""

        return self._h5_file.filename if self._h5_file is not None else None

    @classmethod
    def from_param(
        cls, genes, n_transposons, windows, output_dir, ram=1.2, logger=None
    ):
        """Writable sink for a new file.

        Args:
            genes (GeneData): genes container
            n_transposons (int): transposon count
            windows (list(int)): window sizes
            output_dir (str): directory for output files
            ram (int): upper limit for caching in gigabytes
        """

        logger = logger or logging.getLogger(__name__)
        ram_bytes = int(ram * 1024.0 ** 3)  # MAGIC NUMBER bytes to gigabytes
        check_ram(ram_bytes, logger)

        chromosome = genes.chromosome_unique_id
        filename = chromosome + "_overlap_" + ".h5"  # name files by chromosome
        # NOTE old method below
        # filename = next(tempfile._get_candidate_names()) + ".h5"

        filepath = os.path.join(output_dir, filename)
        config = _OverlapConfigSink(
            genes=genes,
            n_transposons=n_transposons,
            windows=windows,
            filepath=filepath,
            ram_bytes=ram_bytes,
        )
        return cls(config, logger)

    @classmethod
    def from_file(cls, filepath, logger=None):
        """Read only source for an existing file."""

        file_abs = os.path.abspath(filepath)
        if not os.path.isfile(file_abs):
            raise ValueError("input filepath not a file: %s" % file_abs)
        config = _OverlapConfigSource(filepath=file_abs)
        overlap = cls(config, logger)
        return overlap

    @staticmethod
    def left_right_slice(window_idx, gene_idx):
        """Slice for left or right overlap for one window / gene."""

        return (gene_idx, window_idx, slice(None))

    @staticmethod
    def intra_slice(gene_idx):
        """Slice for intra overlap for one window / gene.

        Args:
            gene_idx(int): zero based index for the gene
        """

        # TODO really should just have one window so the interfaces are consistent...
        return (gene_idx, slice(None))

    def start(self):
        """Obtain the resource by opening the file."""

        self._open_dispatcher()

    def stop(self):
        """Release resources by closing the file."""

        self._h5_file.flush()
        self._h5_file.close()
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
            cfg (_OverlapConfigSink): the configuration for a new file.
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
        n_tes = cfg.n_transposons
        # N.B. numpy uses row major by default, iterate over data accordingly
        # this is coupled with OverlapWorker.calculate which iterates
        # this is coupled with the slicing methods of self
        left_right_shape = (
            n_genes,
            n_win,
            n_tes,
        )
        # TODO validate chunk dimensions, 4 * nwin * ntes may be too big pending inputs
        # (and probably already is, but it worked ok...)
        gene_chunks = min(32, n_genes)
        chunks = tuple((gene_chunks, n_win, n_tes))  # MAGIC NUMBER experimental
        self.left = create_set(self._LEFT, left_right_shape, chunks=chunks)
        self.right = create_set(self._RIGHT, left_right_shape, chunks=chunks)
        intra_shape = (n_genes, n_tes)
        self.intra = create_set(self._INTRA, intra_shape)

    def _open_dispatcher(self):
        """Open the file.

        NOTE when upgrading to python 3.8 refactor to use functools.singledispatchmethod.
        """

        if isinstance(self._config, _OverlapConfigSink):
            self._open_new_file(self._config)
        elif isinstance(self._config, _OverlapConfigSource):
            self._open_existing_file(self._config)
        else:
            raise TypeError(
                "expecting {} or {} but got {}".format(
                    _OverlapConfigSink, _OverlapConfigSource, self._config
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

        # FUTURE consider ASCII (fixed len) for storage if non-python clients exist
        # TODO refactor into generic write / read variable len text?
        vlen = h5py.special_dtype(vlen=str)
        dset = self._h5_file.create_dataset(
            self._GENE_NAMES, (len(self.gene_names),), dtype=vlen
        )
        dset[:] = self.gene_names

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


class OverlapWorker:
    """Calculates the overlap values."""

    PROGRESS_CHUNKS = 32  # MAGIC report progress on N genes processed

    # TODO simplify initializer / calculate methods (move args)
    # TODO the args here don't match the docstring
    def __init__(self, overlap_dir, logger=None):
        """Initialize.

        Args:
            overlap_dir (str): A string path of the directory to output te
            overlap files. This comes from the ArgParser obj in density.py and
            defaults to tmp.
            gene_data(GeneData): the genes.
            transposon_data(TransposonData): the transposons.
            windows(list(int)):
        """
        self.root_dir = overlap_dir
        self._logger = logger or logging.getLogger(__name__)

        self._data = None
        self._windows = None
        self._gene_names = None

        self._gene_name_2_idx = None
        self._window_2_idx = None

    def calculate(
        self, genes, transposons, windows, gene_names, stop=None, progress=None
    ):
        """Calculate the overlap for the genes and windows.

        N.B. the number of runs of this multi level loop is intended to be
        limited by the number of input genes / windows rather than refactoring
        into a multiprocess solution.

        Args:
            genes (GeneData): the genes; the chromosome with gene annotations.
            transponsons (TransposonData): the transposable elements.
            windows (list(int)): iterable of window sizes to use.
            gene_names (list(str)): iterable of the genes to use.
            progress (Callable): callback for when a gene name is processed.
        """

        # FUTURE just take in an OverlapJob as input? (plus event, progress bar)

        self._reset(transposons, genes, windows, gene_names)

        # N.B. the order of the loop nesting effects calculation speed
        # it is coupled with the order that the data is stored in OverlapData
        # NOTE consider decoupling?
        # OverlapData is decently sized already but it wouldn't be that much more...
        # OverlapWorker worker would then be empty but needs additions for multiproc
        path = None
        with self._data as sink:
            path = copy.copy(self._data.filepath)
            for gene_name in self._gene_names:
                # TODO check stop event
                gene_datum = genes.get_gene(gene_name)
                g_idx = self._gene_name_2_idx[gene_name]
                out_slice = sink.intra_slice(g_idx)
                sink.intra[out_slice] = Overlap.intra(gene_datum, transposons)
                # NOTE consider iterating on one array at a time? (L / I / R)
                # for dest in l,r: for w in window...
                for window in self._windows:
                    w_idx = self._window_2_idx[window]
                    out_slice = sink.left_right_slice(w_idx, g_idx)
                    sink.left[out_slice] = Overlap.left(gene_datum, transposons, window)
                    sink.right[out_slice] = Overlap.right(
                        gene_datum, transposons, window
                    )
                if progress:
                    progress()

        return path

    def _filter_gene_names(self, gene_names_requested, gene_names_available):
        """Yield only the valid gene names."""

        for name in gene_names_requested:
            # could use list comprehension but we need to log if it fails
            if name not in gene_names_available:
                self._logger.error(" invalid gene name: %s", name)
            else:
                yield name

    def _filter_windows(self, windows):
        """Yield only the valid windows."""

        for window in windows:
            # could use list comprehension but we need to log if it fails
            if window < 0:
                self._logger.error(" invalid window:  %i", window)
            else:
                yield window

    @staticmethod
    def _map_gene_names_2_indices(gene_names):
        """Return the gene names mapped to an index for the array."""

        return {name: idx for idx, name in enumerate(gene_names)}

    @staticmethod
    def _map_windows_2_indices(windows):
        """Return the window values mapped to an index for the array."""

        return {win: idx for idx, win in enumerate(windows)}

    def _reset(self, transposons, genes, windows, gene_names):
        """Initialize overlap data; mutates self."""

        gene_names_filtered = self._filter_gene_names(gene_names, genes.names)
        self._gene_names = list(gene_names_filtered)
        self._gene_name_2_idx = self._map_gene_names_2_indices(self._gene_names)
        self._windows = list(self._filter_windows(windows))
        self._window_2_idx = self._map_windows_2_indices(self._windows)

        n_te = transposons.number_elements
        self._data = OverlapData.from_param(genes, n_te, self._windows, self.root_dir)
