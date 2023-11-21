#!/usr/bin/env python3

"""
Unit test refactor of DensityData.
"""

import logging
import coloredlogs
from copy import deepcopy
from contextlib import contextmanager
import shutil
import tempfile

import numpy as np
import pytest
import h5py

from transposon.density2 import DensitySubset, _DensitySubsetConfig


WINDOWS = [10, 100, 1000]
CHROME_ID = "fake_chromosome_id"
TE_NAMES = ["fake_te_0", "fake_te_1", "fake_te_2"]
GENE_NAMES = ["fake_gene_0", "fake_gene_1"]
LOGGER = logging.getLogger(__name__)
coloredlogs.install(level=logging.DEBUG)


@pytest.fixture
def temp_dir():
    """Yields a temporary directory.

    NB this is used to run once per file exec, i.e fixture level 'module',
    so that the temporary files are removed automatically after all tests complete,
    and so that we can re-open the files that got written.
    """

    with tempfile.TemporaryDirectory() as dir:
        yield dir


@contextmanager
def temp_h5_context(temp_dir):
    """Yields an open HDF5 file, writes on close."""

    with tempfile.NamedTemporaryFile(dir=temp_dir, delete=False) as temp:
        # MAGIC h5py convention, 'a' is append
        with h5py.File(temp.name, 'a') as file:
            yield file


@pytest.fixture
def default_config():
    """A configuration for testing."""

    conf = _DensitySubsetConfig(
        windows=WINDOWS,
        te_names=TE_NAMES,
        gene_names=GENE_NAMES
    )
    return conf


@pytest.fixture
def default_density(temp_dir, default_config):

    with temp_h5_context(temp_dir) as file:
        yield DensitySubset(file, CHROME_ID, default_config)


@contextmanager
def read_density_from_other(density):
    """Create a new DensitySubset from another.

    Intended to simplify writing / reading data from an H5.
    """

    density.file.flush()
    # NB copy the file so that changes to one file aren't applied to the other
    with tempfile.TemporaryDirectory() as dir:
        with tempfile.NamedTemporaryFile(dir=dir, delete=True) as temp:
            shutil.copyfile(density.filename, temp.name)
            # MAGIC h5py convention, 'a' is append
            with h5py.File(temp, 'a') as new_h5:
                new_density = DensitySubset(new_h5, CHROME_ID, density.cfg)
                yield new_density


def test_init(default_density):
    """Can a new file be created?"""

    pass  # pass b/c we just want the fixture to run


def test_write_read_shared(default_density):
    """Can we change one file without affecting the other?

    NB we want to be able to write and read data from files reliably.
    We do this by writing / modifying / reading files in varous orders.
    If we read from an existing H5 file however, the data is shared.
    This precludes us from writing, reading, and then modifying *one* of the files,
        b/c the data is shared between the first and second.

    This test confirms that our helper function is valid,
        b/c we copy the file and we can reliably change the data.
    """

    with read_density_from_other(default_density) as other:
        # NB we should be able to change the first w/o it affecting the other
        default_density._group[default_density._WINDOWS][0] += 42
        assert (default_density.windows != other.windows).any()


def test_write_read_genes(default_density):
    """Can we read gene names back?"""

    with read_density_from_other(default_density) as other:
        default_density.gene_names == other.gene_names


def test_read_check_mismatch_genes(default_density):
    """Does it raise if the genes we expect aren't there?"""

    other_cfg = deepcopy(default_density.cfg)
    other_cfg.gene_names.append("unexpected_gene")
    with pytest.raises(TypeError):
        other_density = DensitySubset(default_density.file, CHROME_ID, other_cfg)


def test_gene_value_mismatch(default_density):
    """Does it raise if the genes have the right shape but wrong values?"""

    other_cfg = deepcopy(default_density.cfg)
    other_cfg.gene_names[0] = "unexpected_gene"
    with pytest.raises(ValueError):
        other_density = DensitySubset(default_density.file, CHROME_ID, other_cfg)


def test_write_read_transposons(default_density):
    """Can we read transposon names back?"""

    with read_density_from_other(default_density) as other:
        default_density.transposon_names == other.transposon_names

def test_read_check_mismatch_transposons(default_density):
    """Does it raise if the transposons we expect aren't there?"""

    other_cfg = deepcopy(default_density.cfg)
    other_cfg.te_names.append("unexpected_transposon")
    with pytest.raises(TypeError):
        other_density = DensitySubset(default_density.file, CHROME_ID, other_cfg)


def test_te_value_mismatch(default_density):
    """Does it raise if the tes have the right shape but wrong values?"""

    other_cfg = deepcopy(default_density.cfg)
    other_cfg.te_names[0] = "unexpected_te"
    with pytest.raises(ValueError):
        other_density = DensitySubset(default_density.file, CHROME_ID, other_cfg)


def test_write_read_windows(default_density):
    """Can we read gene names back?"""

    with read_density_from_other(default_density) as other:
        assert (default_density.windows == other.windows).all()

def test_read_check_mismatch_windows(default_density):
    """Does it raise if the genes we expect aren't there?"""

    other_cfg = deepcopy(default_density.cfg)
    other_cfg.windows.append(9191)
    with pytest.raises(TypeError):
        other_density = DensitySubset(default_density.file, CHROME_ID, other_cfg)


def test_window_value_mismatch(default_density):
    """Does it raise if the genes have the right shape but wrong values?"""

    other_cfg = deepcopy(default_density.cfg)
    other_cfg.windows[0] = 9999
    with pytest.raises(ValueError):
        DensitySubset(default_density.file, CHROME_ID, other_cfg)


def test_write_read_left(default_density):
    """Can we write / read the left densities?"""

    with read_density_from_other(default_density) as other:
        assert (default_density.left[:] == other.left[:]).all()


def test_write_read_intra(default_density):
    """Can we write / read the intra densities?"""

    with read_density_from_other(default_density) as other:
        assert (default_density.intra[:] == other.intra[:]).all()


def test_write_read_right(default_density):
    """Can we write / read the right densities?"""

    with read_density_from_other(default_density) as other:
        assert (default_density.right[:] == other.right[:]).all()


def test_write_read_bitmap(default_density):
    """Can we write / read the bitmap?"""

    with read_density_from_other(default_density) as other:
        assert (default_density._bitmap[:] == other._bitmap[:]).all()


def test_write_read_bitmap_noequals(default_density):
    """If we flip one bitmap, is the other unchanged?"""

    with read_density_from_other(default_density) as other:
        other._bitmap[:] = np.invert(default_density._bitmap[:])
        assert (default_density._bitmap[:] != other._bitmap[:]).all()
