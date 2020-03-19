#!/usr/bin/env python3

__author__ = "Michael Teresi"

"""
The no. base pairs that the transposable elements overlap the window for a gene.
"""

from enum import Enum, unique
import logging
import os
from functools import partial
import tempfile

import numpy as np
import h5py

class Overlap():
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
        lower_bound = np.minimum(g_start, transposons.stops)
        upper_bound = np.maximum(g_stop, transposons.starts)
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

class _OverlapDataSink():
    """Destination for overlap calculations.

    Contains a numpy array for each `Overlap.Direction` buffered to a file.
    The file is h5py and contains data sets named by the keys in `Overlap.Direction`.
    """

    # FUTURE make implementation that just keeps it in ram?

    DTYPE = np.float32
    COMPRESSION = 'lzf'  # MAGIC NUMBER experimental
    LEFT = Overlap.Direction.LEFT.name
    RIGHT = Overlap.Direction.RIGHT.name
    INTRA = Overlap.Direction.INTRA.name
    GENE_NAMES = 'GENE_NAMES'

    def __init__(self, gene_names, n_transposons, n_win, out_dir, logger=None):

        self._logger = logger or logging.getLogger(__name__)
        filename = next(tempfile._get_candidate_names()) + '.h5'
        self.filepath = os.path.join(out_dir, filename)
        self.h5_file = None
        self.left = None
        self.intra = None
        self.right = None
        self.gene_names = list(gene_names)  # FUTURE require generator func instead
        self._n_genes = len(gene_names)
        self._n_tes = int(n_transposons)
        self._n_win = int(n_win)

    @classmethod
    def from_param(cls, gene_names, n_transposons, n_win, ram=2, logger=None):
        """Empty container when writing data to a new file."""
        raise NotImplementedError()

    @classmethod
    def from_file(cls, filepath, logger=None):
        """Target an existing data file."""
        raise NotImplementedError()

    @staticmethod
    def left_right_slice(window, gene):
        """Slice for left or right overlap for one window / gene."""

        return (gene, window, slice(None))

    @staticmethod
    def intra_slice(gene):
        """Slice for intra overlap for one window / gene."""

        return (gene, slice(None))

    def _write_gene_names(self, h5_file, gene_names):
        """Write list of gene names the file."""

        # FUTURE consider ASCII (fixed len) for storage if non-python clients exist
        vlen = h5py.special_dtype(vlen=str)
        dset = h5_file.create_dataset(self.GENE_NAMES, (self._n_genes,), dtype=vlen)
        dset[:] = gene_names

    @classmethod
    def _read_gene_names(cls, h5_file):
        """Numpy ndarray of S type to generator of string."""

        gene_array = h5_file[cls.GENE_NAMES][:]
        return (gene_array[i] for i in range(gene_array.size))

    def _write(self):
        raise NotImplementedError()

    def _read(self):
        raise NotImplementedError()

    def __enter__(self):
        # FUTURE do either reading or writing depending on file mode

        # TODO parametrize, allow one to specify pending their ram availability
        # MAGIC NUMBER about 2GB for ram cache, effects overall speed
        self.h5_file = h5py.File(self.filepath, 'w', rdcc_nbytes=2*1024*1024**2)
        # MAGIC NUMBER lz4 compression looks good so far...
        create_set = partial(self.h5_file.create_dataset,
                             dtype=self.DTYPE,
                             compression=self.COMPRESSION)
        # N.B. numpy uses row major by default
        lr_shape = (self._n_genes, self._n_win, self._n_tes)
        # TODO validate chunk dimensions, 4 * nwin * ntes may be too big pending inputs
        # (and probably already is, but it worked ok...)
        chunks = tuple((32, self._n_win, self._n_tes))  # MAGIC NUMBER experimental
        self.left = create_set(self.LEFT, lr_shape, chunks=chunks)
        self.right = create_set(self.RIGHT, lr_shape, chunks=chunks)
        i_shape = (self._n_genes, self._n_tes)
        self.intra = create_set(self.INTRA, i_shape)

        self._write_gene_names(self.h5_file, self.gene_names)

        # TODO output window map
        # TODO output chromosome ID

        return self

    def __exit__(self, exc_type, exc_val, exc_traceback):
        self.h5_file.flush()
        self.h5_file.close()
        self.left = None
        self.right = None
        self.intra = None


# REFACTOR rename to OverlapDataCalculator
class OverlapData():
    """Contains the overlap between the genes and transposable elements."""

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

        self._gene_name_2_idx = None  # NOTE Michael can you explain sometime?
        self._window_2_idx = None

    def calculate(self, genes, transposons, windows, gene_names, progress=None):
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

        self._reset(transposons, genes, windows, gene_names)

        with self._data as sink:
            for gene_name in self._gene_names:
                gene_datum = genes.get_gene(gene_name)
                g_idx = self._gene_name_2_idx[gene_name]
                destination = sink.intra_slice(g_idx)
                sink.intra[destination] = Overlap.intra(gene_datum, transposons)
                for window in self._windows:
                    left = Overlap.left(gene_datum, transposons, window)
                    right = Overlap.right(gene_datum, transposons, window)
                    w_idx = self._window_2_idx[window]
                    destination = sink.left_right_slice(w_idx, g_idx)
                    sink.left[destination] = left
                    sink.right[destination] = right
                if progress is not None:
                    progress()

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
        """Initialize overalap data; mutates self."""

        gene_names_filtered = self._filter_gene_names(gene_names, genes.names)
        self._gene_names = list(gene_names_filtered)
        self._gene_name_2_idx = self._map_gene_names_2_indices(self._gene_names)
        self._windows = list(self._filter_windows(windows))
        self._window_2_idx = self._map_windows_2_indices(self._windows)

        n_win = len(self._windows)
        n_gene = len(self._gene_names)
        n_te = transposons.number_elements
        self._data = _OverlapDataSink(self._gene_names, n_te, n_win, self.root_dir)
