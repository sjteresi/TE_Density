#!/usr/bin/env python

"""
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import re
import pandas as pd


# import matplotlib.pyplot as plt
# import seaborn as sns

# from examples.Arabidopsis.src.import_Arabidopsis_gene_anno import import_genes
from transposon.density_data import DensityData
from transposon.gene_data import GeneData
from transposon.import_filtered_genes import import_filtered_genes

def extract_centromeric_pericentromeric_genes():
    # NOTE
    # Ultimately, I need to compare the density values of certain lists(?) of
    # genes; there will be two categories of genes. Genes that are
    # centromeric/pericentromeric, and "regular" genes.
    # I have the MAGIC numbers that signify the base-pair boundaries of
    # centromeric genes, so every gene NOT in that boundary will be a "regular"
    # gene.
    # The issue will be extracting genes from the HDF5.
    # Steps:
        # Extract gene names (in list form?) from GeneData that have the
        # correct position value.

        # Extract the indices of the genes from the HDF5
        # This may require a new function in DensityData because I currently
        # only have a function that does one gene at a time...

        # Unsure:
            # Make a bar plot similar to the bar plots in the Blueberry example?
            # Get a gene expression dataset and verify with Pat?

if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.join(dir_main, "..", "results/graphs")
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "gene_input_file", type=str, help="Parent path of cleaned gene annotation file"
    )

    parser.add_argument(
        "density_data_folder",
        type=str,
        help="Parent path of folders containing TE Density results",
    )

    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=output_default,
        help="parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.gene_input_file = os.path.abspath(args.gene_input_file)
    args.density_data_folder = os.path.abspath(args.density_data_folder)
    args.output_dir = os.path.abspath(args.output_dir)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -----------------------------------

    cleaned_genes = import_filtered_genes(args.gene_input_file, logger)
    gene_dataframe_list = [
        dataframe for k, dataframe in cleaned_genes.groupby("Chromosome")
    ]
    gene_data_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in gene_dataframe_list
    ]

    processed_dd_data = DensityData.from_list_gene_data_and_hdf5_dir(
        gene_data_list, args.density_data_folder, "Arabidopsis_(.*?).h5", logger
    )
