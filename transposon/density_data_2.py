

from dataclasses import dataclass
from contextlib import ExitStack
from enum import Enum
from tempfile import NamedTemporaryFile

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
    """Contains the calculated transposable element density.

    Is a context manager; wraps an HDF5 file.

    Contains the densities for the superfamilies and orders.
    """

    def __init__(self, gene_data, te_data, windows, filepath=None, logger=None):
        """Initialize.

        Args:
            gene_data():
            te_data():
            windows():
            filepath():
            logger():
        """

        super().__init__()

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

    def __exit__(self, exc_type, exc, exc_tb):
        pass
        # TODO close file

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
            # NB SpooledTemporaryFile doesn't appear to work w/ h5py by default
            tf = NamedTemporaryFile()
            self.enter_context(tf)
            self._path = tf.name
        # MAGIC a is: Read/write if exists, create otherwise
        hfile = h5py.File(self._path, 'a')
        self.enter_context(hfile)
        self._hdf5 = hfile

    def generate_jobs(self):
        """Yield DensityJob instances, for each missing entry."""
        pass  # TODO

    def reduce_result(self):
        """Save the DensityResult in the matrix."""
        pass  # TODO


