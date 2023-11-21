

from dataclasses import dataclass
from contextlib import ExitStack
from enum import Enum
from tempfile import NamedTemporaryFile
import logging

import h5py

from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.density2 import DensitySubset, _DensitySubsetConfig



class DensityData(ExitStack):
    """Read and write the transposable element densities

    Wraps an HDF5 file.
    Is a context manager.
        If no filepath is provided, a temporary file will be used.
        NB the temporary file will be deleted on exit, so flush and copy as necessary.
    Provides methods to read/update densities.

    Intialize your instance, enter, read/write, exit.
    """

    # TODO take in OverlapData and NOT windows, b/c that has the windows in it
    def __init__(self, gene_data, te_data, windows, filepath=None, logger=None):
        """Initialize.

        NB if reading an existing file, the GeneData, TransposonData, and windows,
            _must_ match the data in the existing file or the program will abort,
        Args:
            gene_data(): the genes container
            te_data(): the transposable elements container
            windows(): the window counts
            filepath(str): path to a file to read/write, creates a TemporaryFile if None
            logger(Logger): logging instance, creates one under class name if None
        """

        super().__init__()
        self._logger: logging.Logger = logger or logging.getLogger(__name__)
        self._path: str = filepath
        self._hdf5: h5py.File = None

        self.gene_data: GeneData = gene_data
        self.te_data: TransposonData = te_data
        self.windows: list(str) = windows

        self.rho_superfamily: DensitySubset = None
        self.rho_order: DensitySubset = None

    @property
    def chromosome_id(self):
        """Unique chromosome identifier for all the genes available."""

        return str(self.te_data.chromosome_unique_id)

    @property
    def path(self):
        """The current filepath."""

        return self._path

    def flush(self):
        """Write buffers to file."""

        self._hdf5.flush()

    def generate_jobs(self):
        """Yield DensityJob instances, for each missing entry."""

        pass  # TODO

    def generate_workers(self, overlap_data):

        chromosome_uid = self.chromosome_id

        transposons = self.te_data.superfamilies
        rho_sub = self.rho_superfamily  # contains the windows/te_names/gene_names
        overlap = overlap_data

        transposons = self.te_data.oders
        rho_sub = self.rho_order
        overlap = overlap_data

    def reduce_result(self):
        """Save the DensityResult in the matrix."""

        pass  # TODO

    def __enter__(self):
        """Acquire resources.

        Open the file, initialize the superfamily and order data.
        """

        self._open()

        self.rho_superfamily = self._init_subset(
            "superfamily",
            self.te_data.superfamily_name_set,
        )
        self.rho_order = self._init_subset(
            "order",
            self.te_data.order_name_set,
        )

        return self

    def __exit__(self, *args):
        """Release resources."""

        super().__exit__(*args)
        self._hdf5 = None
        self.superfamily = None
        self.order = None

    def _init_subset(self, prefix, te_names):
        """Initialize a subset of TE Densities, e.g. order or superfamily.

        Args:
            prefix(str):
            te_name(list(str)):
        """

        gene_names = self.gene_data.names
        cfg = _DensitySubsetConfig(
            windows=list(self.windows),
            te_names=te_names,
            gene_names=list(n for n in gene_names)
        )
        return DensitySubset(
            self._hdf5,
            str(prefix),
            cfg
        )

    def _open(self):
        """Open the H5File.

        If path is None, creates a temporary file.
        Mutates the hdf5 instance.
        """

        if self._path is None:
            msg = "path is None, using a temporary file, data will be lost on exit!"
            self._logger.warning(msg)
            # NOTE it's more secure to first create a temporary directory
            # then create the file
            handle = NamedTemporaryFile()
            self.enter_context(handle)
            self._path = handle.name
        # MAGIC a is: Read/write if exists, create otherwise
        # MAGIC latest: use most performant version over backwards compatibility
        hfile = h5py.File(self._path, 'a', libver='latest')
        self._logger.debug("opened %s", self._path)
        # NB close the h5py.File _before_ the file obj, so enter it second
        self.enter_context(hfile)
        self._hdf5 = hfile
