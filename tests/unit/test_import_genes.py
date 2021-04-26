#!/usr/bin/env python3

"""
Unit test GeneDatum
"""

__author__ = "Scott Teresi"

import os
import pytest
import numpy as np
import pandas as pd
from io import StringIO

from import_scripts.import_strawberry_gene_anno import import_genes


col_names = [
    "Chromosome",
    "Software",
    "Feature",
    "Start",
    "Stop",
    "Score",
    "Strand",
    "Frame",
    "FullName",
]

col_to_use = [
    "Chromosome",
    "Software",
    "Feature",
    "Start",
    "Stop",
    "Strand",
    "FullName",
]

gene_anno_path = "tests/input_data/Test_Gene_Anno_Float_Conversion.tsv"


@pytest.fixture
def import_as_float32():
    gene_anno = pd.read_csv(
        gene_anno_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        dtype={"Start": "float32", "Stop": "float32"},
        comment="#",
    )
    gene_anno = gene_anno[gene_anno.Feature == "gene"]  # drop non-gene rows
    # clean the names and set as the index (get row wrt name c.f. idx)
    gene_anno["Gene_Name"] = gene_anno["FullName"].str.extract(r"ID=(.*?);")
    gene_anno.set_index("Gene_Name", inplace=True)
    gene_anno = gene_anno.drop(columns=["FullName", "Software"])
    gene_anno.Strand = gene_anno.Strand.astype(str)
    gene_anno["Length"] = gene_anno.Stop - gene_anno.Start + 1
    return gene_anno


@pytest.fixture
def import_as_float64():
    gene_anno = pd.read_csv(
        gene_anno_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        dtype={"Start": "float64", "Stop": "float64"},
        comment="#",
    )
    gene_anno = gene_anno[gene_anno.Feature == "gene"]  # drop non-gene rows
    # clean the names and set as the index (get row wrt name c.f. idx)
    gene_anno["Gene_Name"] = gene_anno["FullName"].str.extract(r"ID=(.*?);")
    gene_anno.set_index("Gene_Name", inplace=True)
    gene_anno = gene_anno.drop(columns=["FullName", "Software"])
    gene_anno.Strand = gene_anno.Strand.astype(str)
    gene_anno["Length"] = gene_anno.Stop - gene_anno.Start + 1
    return gene_anno


true_start_list = [
    41.0,
    5556.0,
    8487.0,
    9361.0,
    11127.0,
    84598.0,
    117287120.0,
    118974314397.0,
    22456307315831.0,
    88877765432319026.0,
]
true_stop_list = 2


def test_small_numbers_to_float32(import_as_float32):
    """
    Using float64 the import script produce the right start and stop values when
    importing a gene annotation with small initial values?
    """
    assert import_as_float32.Start.to_list()[0:5] == true_start_list[0:5]


def test_large_numbers_to_float32(import_as_float32):
    """
    Using float64 the import script produce the right start and stop values when
    importing a gene annotation with large initial values?
    """
    with pytest.raises(AssertionError):
        assert import_as_float32.Start.to_list()[4:] == true_start_list[4:]


def test_large_numbers_to_float64(import_as_float64):
    """
    Using float64 the import script produce the right start and stop values when
    importing a gene annotation with large initial values?
    """
    assert import_as_float64.Start.to_list()[4:] == true_start_list[4:]


if __name__ == "__main__":
    pytest.main(["-svv", __file__])  # for convenience
