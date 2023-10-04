

from dataclasses import dataclass
from contextlib import ExitStack
from enum import Enum
from tempfile import NamedTemporaryFile
import logging

import h5py

from transposon.density2 import _DensitySubset, _DensitySubsetConfig


@dataclass
class DensityJob:

    chrome_id: str
    group_type: str  # superfam or order, may not need this one?
    windows: list[int]
    direction: Enum  # left, intra, right
    gene_name: str  # we can get the gene idx from this
                    # also needed to get the divisor value


class DensityResult:

    job: DensityJob


class DensityWorker:


    def __init__(self):

        # static
        self.chrome_id = None  # pseudo-molecule unique id
        self.group_type = None  # superfamily or order or etc

        self.windows = list()
        self.gene_names = list()  # str
        self.group_names = list() # str, either superfamily or order

        #self.gene_datum
        #self.overlap_data


class DensityData(ExitStack):
    """Read and write the transposable element densities

    Wraps an HDF5 file.
    Is a context manager.
        If no filepath is provided, a temporary file will be used.
    Provides methods to read/update densities.

    Intialize your instance, enter, read/write, exit.
    """

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
        self._logger = logger or logging.getLogger(__name__)

        self._path = filepath
        self._hdf5 = None

        self.gene_data = gene_data
        self.te_data = te_data
        self.windows = windows

        self.superfamily = None
        self.order = None

    @property
    def path(self):
        """The current filepath."""

        return self._path

    def flush(self):
        """Write buffers to file."""

        self._hdf5.flush()

    def write(self, path):
        """Write the data to a new filepath."""

        raise NotImplementedError()

    def __enter__(self):
        """Acquire resources.

        Open the file, initialize the superfamily and order data.
        """

        self._open()

        self.superfamily = self._init_subset(
            "superfamily",
            self.te_data.superfamily_name_set,
        )
        self.order = self._init_subset(
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
        """Initialize the type of TE density data"""

        gene_names = self.gene_data.names
        cfg = _DensitySubsetConfig(
            windows=list(self.windows),
            te_names=te_names,
            gene_names=[n for n in gene_names]
        )
        return _DensitySubset(
            self._hdf5,
            str(prefix),
            cfg
        )

    def _open(self):
        """Open the H5File.

        If path is None, creates a temporary file without a filepath.
        """

        if self._path is None:
            # NOTE it's more secure to first create a temporary directory
            # then create the file
            handle = NamedTemporaryFile()
            self.enter_context(handle)
            self._path = handle.name
        # MAGIC a is: Read/write if exists, create otherwise
        # MAGIC latest: use most performant version over backwards compatibility
        hfile = h5py.File(self._path, 'a', libver='latest')
        self.enter_context(hfile)
        self._hdf5 = hfile

    def generate_jobs(self):
        """Yield DensityJob instances, for each missing entry."""

        pass  # TODO

    def reduce_result(self):
        """Save the DensityResult in the matrix."""

        pass  # TODO
