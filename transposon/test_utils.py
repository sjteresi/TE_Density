#!/usr/bin/env python3

"""
Sundry functions intended for testing and development.
"""

import h5py

import pytest
import tempfile


@pytest.fixture
def temp_dir():
    """Yields a temporary directory."""

    # NOTE using a scope of 'module' doesn't appear to work
    # if you import this, so default scope is used
    with tempfile.TemporaryDirectory() as dir:
        yield dir


@pytest.fixture
def temp_h5_file(temp_dir):  # FUTURE could be a contextmanager rather than fixture
    """Yields a temporary HDF5 file."""

    file = tempfile.NamedTemporaryFile(dir=temp_dir, suffix=".h5")
    with file as temp:
        yield temp.name


@pytest.fixture
def temp_h5_context():  # FUTURE could be a contextmanager rather than fixture
    """Yields an open HDF5 file, writes on close."""

    # MAGIC 1KB as a reasonable default
    with tempfile.SpooledTemporaryFile(max_size=1024*1024) as temp:
        # MAGIC h5py convention, 'a' is append
        with h5py.File(temp, 'a') as file:
            yield file
            file.flush()  # NOTE likely don't need this
            file.close()  # NOTE likely don't need this
