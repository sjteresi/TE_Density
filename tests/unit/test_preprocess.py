#!/usr/bin/env python3

"""
Test preprocess
"""

__author__ = "Scott Teresi"

import pytest
import tempfile
import pandas as pd

from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.preprocess import PreProcessor


@pytest.fixture
def gene_file(request):
    """path to simple gene annotation file for testing."""
    return request.param


@pytest.fixture
def transposon_file(request):
    """
    path to simple TE annotation file for testing. This file is special in
    that it does not have corresponding chromosome IDs to the gene file.
    """
    return request.param


@pytest.fixture()
def temp_dir():
    """Temporary directory."""

    with tempfile.TemporaryDirectory() as dir:
        yield dir


@pytest.mark.parametrize(
    "gene_file",
    [
        "tests/input_data/Test_Preprocess_Cleaned_Genes.tsv",
    ],
)
@pytest.mark.parametrize(
    "transposon_file",
    [
        "tests/input_data/Test_Preprocess_Cleaned_TEs_Unequal_Chrom_IDs.tsv",
        "tests/input_data/Test_Preprocess_Cleaned_TEs_Unequal_Number_Chrom.tsv",
    ],
)
def test_validate_split(gene_file, transposon_file, temp_dir):
    """
    Does the program fail elegantly when you give it BAD input annotations that
    do not have the right corresponding chromosomes?

    Args:
        gene_file (str): string to filepath of a "cleaned" gene annotation
        transposon_file (str): string to filepath of a "cleaned" TE annotation
        temp_dir (str): string to filepath for temporary output dir
    """
    preprocessor_obj = PreProcessor(
        gene_file,
        transposon_file,
        temp_dir,
        False,  # reset h5 arg
        "fake_genome_id",  # MAGIC genome ID arg
        True,  # revise TE arg
    )
    gene_frame = preprocessor_obj._filter_genes()
    transposon_frame = preprocessor_obj._filter_transposons()
    transposon_frame = preprocessor_obj._revise_transposons(transposon_frame)
    with pytest.raises(ValueError) as exc:
        preprocessor_obj._split_wrt_chromosome(gene_frame, transposon_frame)


if __name__ == "__main__":
    pytest.main(["-s", __file__])
