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
    parser = argparse.ArgumentParser(description="start to analyze TE Density data")

    parser.add_argument(
        "gene_data_cache_folder",
        type=str,
        help="""parent path to the folder containing the gene data cache files""",
    )

    parser.add_argument(
        "density_data_folder",
        type=str,
        help="Parent path of folder containing ONLY the TE Density results",
    )

    parser.add_argument(
        "cleaned_gene_annotation",
        type=str,
        help="path to your cleaned gene annotation (.tsv)",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.gene_data_cache_folder = os.path.abspath(args.gene_data_cache_folder)
    args.density_data_folder = os.path.abspath(args.density_data_folder)
    args.cleaned_gene_annotation = os.path.abspath(args.cleaned_gene_annotation)

    # NB just for logging arguments to import_filtered command and DensityData
    # initialization
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # --------------------------------------------------
    # NOTE this object is a list of DensityData instances
    processed_dd_data = DensityData.from_list_genedata_dir_and_hdf5_dir(
        args.gene_data_cache_folder, args.density_data_folder, logger
    )
    # Note this is a pandas object of all of your genes
    cleaned_genes = import_filtered_genes(args.cleaned_gene_annotation, logger)

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

    # TODO give more explicit examples so the users don't have to go digging as
    # much through the example code.

    # Example 1: Get the quick values for a gene
    # Here we will use the first DensityData object, which is probably
    # chromosome 1 because it was sorted alphabetically
    window_idx = 0  # with default windows, index 1 this should be 1KB
    your_gene_name = "AT1G01010"
    print(processed_dd_data[0])
    info_for_your_favorite_gene = processed_dd_data[0].info_of_gene(
        your_gene_name, window_idx
    )
    print(info_for_your_favorite_gene)

    # Example 2: Add the HDF5 indices and TE vals to the gene pandaframe
    # TODO
