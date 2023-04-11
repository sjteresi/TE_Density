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
