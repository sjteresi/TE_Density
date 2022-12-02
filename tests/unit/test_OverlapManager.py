#!/usr/bin/env python3

"""
Unit test OverlapManager.
"""

__author__ = "Michael Teresi"

import pytest
import os
import tempfile

from transposon.overlap_manager import OverlapManager
from transposon.transposon_data import TransposonData
from transposon.gene_data import GeneData

N_SUBGENES = 4  # MAGIC arbitrary, enough for testing
GENOME_ID = "FAKE_GENOME_ID"


@pytest.fixture(scope="session")
def temp_dir():
    """Temporary directory."""

    with tempfile.TemporaryDirectory() as dir:
        yield dir


@pytest.fixture(scope="session")
def temp_filenames(temp_dir):
    names = [next(tempfile._get_candidate_names()) for _ in range(N_SUBGENES)]
    paths = [os.path.join(temp_dir, n) for n in names]
    return paths


@pytest.fixture(scope="session")
def sub_genes(temp_filenames):

    genes = [GeneData.mock(genome_id=GENOME_ID) for _ in range(len(temp_filenames))]
    return [gene.write(path) for gene, path in zip(genes, temp_filenames)]


@pytest.fixture(scope="session")
def sub_transposons(temp_filenames):

    genes = [
        TransposonData.mock(genome_id=GENOME_ID) for _ in range(len(temp_filenames))
    ]
    return [gene.write(path) for gene, path in zip(genes, temp_filenames)]


# def test_init_nothrow(sub_genes, sub_transposons):
# pass


if __name__ == "__main__":
    pytest.main(["-s", __file__])  # for convenience
