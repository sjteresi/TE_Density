#!/usr/bin/env/python

"""
Barebones initialization of DensityData class intended for demonstration
purposes.
"""

__author__ = "Scott Teresi"

import argparse
import pandas as pd
import os
import logging
import coloredlogs

from transposon.gene_data import GeneData
from transposon.density_data import DensityData
from transposon.import_filtered_genes import import_filtered_genes
from transposon.density_utils import (
    add_hdf5_indices_to_gene_data_from_list_hdf5,
    add_te_vals_to_gene_info_pandas_from_list_hdf5,
    info_of_gene,
)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="generate graphs")

    parser.add_argument(
        "cleaned_gene_annotation",
        type=str,
        help="""parent path to the cleaned gene data file that was derived from
        preprocessing""",
    )

    parser.add_argument(
        "density_data_folder",
        type=str,
        help="Parent path of folder containing ONLY the TE Density results",
    )

    parser.add_argument(
        "hdf5_string",
        type=str,
        help="""Regex string rule for extracting the genome name and chromosome
        ID from the h5 file. Please refer to the
        from_list_gene_data_and_hdf5_dir method of DensityData for more
        documentation""",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.cleaned_gene_annotation = os.path.abspath(args.cleaned_gene_annotation)
    args.density_data_folder = os.path.abspath(args.density_data_folder)

    # NB just for logging arguments to import_filtered command and DensityData
    # initialization
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    cleaned_genes = import_filtered_genes(args.cleaned_gene_annotation, logger)

    gene_dataframe_list = [
        dataframe for k, dataframe in cleaned_genes.groupby("Chromosome")
    ]

    # NB this is how I get a list of GeneData objects so that I may
    # initialize DensityData further on
    # NB MAGIC int to get chromosome ID
    gene_data_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in gene_dataframe_list
    ]

    # NOTE this object is a list of DensityData instances
    # Each instance corresponds to one pseudmolecule of TE Density data.
    # If you are having trouble, please examine the documentation of the
    # classmethod 'from_list_gene_data_and_hdf5_dir'.
    processed_dd_data = DensityData.from_list_gene_data_and_hdf5_dir(
        gene_data_list, args.density_data_folder, args.hdf5_string, logger
    )

    # NOTE, this adds the indices of the genes in the HDF5 datasets to a pandas
    # dataframe, this is used later on to access the density data for each gene
    gene_frame_with_hdf5_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_genes, processed_dd_data
    )

    # NOTE, now we can add columns to the pandas dataframe which are the
    # density values for a specific TE type and window combo. This right here
    # adds a column called "LTR_400_Upstream" to the pandas dataframe, this
    # column represents the TE density values for each gene, for the 400 bp
    # upstream window. 400 bp window is used here because I use a smaller
    # testing set of windows for development.
    gene_frame_with_te_vals_of_interest = (
        add_te_vals_to_gene_info_pandas_from_list_hdf5(
            gene_frame_with_hdf5_indices,
            processed_dd_data,
            "Order",
            "LTR",
            "Upstream",
            500,
        )
    )

    print(gene_frame_with_te_vals_of_interest)
    print(gene_frame_with_te_vals_of_interest)
    print(info_of_gene(processed_dd_data[0], "AT1G01ad10", 0))

    # ------------------------------------------------------------
    # NOTE now the user may begin analyzing their data. For more examples
    # please examine the code used for the Arabidopsis, Blueberry, Human, and
    # Rice datasets. But here are some select ones:
