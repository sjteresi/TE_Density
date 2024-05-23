#!/usr/bin/env/python

"""
Barebones initialization of DensityData class intended for demonstration
purposes.
"""

__author__ = "Scott Teresi"

import pandas as pd
import numpy as np
import os
import argparse
import logging
import coloredlogs

from transposon.import_filtered_genes import import_filtered_genes
from transposon.gene_data import GeneData
from transposon.density_data import DensityData

from transposon.density_utils import (
    add_hdf5_indices_to_gene_data_from_list_hdf5,
    add_te_vals_to_gene_info_pandas_from_list_hdf5,
    add_te_vals_to_gene_info_pandas,
    get_specific_slice,
    add_hdf5_indices_to_gene_data,
    info_of_gene,
)


def get_gene_data_as_list(cleaned_genes):
    """
    Take a cleaned genes annotation file from TE Density
    (import_filtered_genes) and break it into a list of GeneData objects by
    chromosome ID. This is used to initialize all of the DensityData objects in
    a list.

    Args:
        cleaned_genes (pandas.DataFrame)
            Index:
                Name: Gene_Name, strings of gene names
            Columns:
                Name: Chromosome, object
                Name: Feature, object
                Name: Start, float64
                Name: Stop, float64
                Name: Strand, object
                Name: Length, float64

    Returns:
        genedata_list (list of GeneData)
    """
    # MAGIC group by column Chromosome
    gene_dataframe_list = [
        dataframe for k, dataframe in cleaned_genes.groupby("Chromosome")
    ]

    # MAGIC initialize GeneData iteratively using the magic unique chromosome
    genedata_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in gene_dataframe_list
    ]
    return genedata_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="start to analyze TE Density data")

    parser.add_argument(
        "cleaned_gene_annotation",
        type=str,
        help="path to your cleaned gene annotation (.tsv)",
    )

    parser.add_argument(
        "density_data_dir",
        type=str,
        help="Parent path of folder containing ONLY the TE Density results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.cleaned_gene_annotation = os.path.abspath(args.cleaned_gene_annotation)
    args.density_data_dir = os.path.abspath(args.density_data_dir)

    # NB just for logging arguments to import_filtered command and DensityData
    # initialization
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # --------------------------------------------------
    # Read cleaned genes for the given genome as pandas
    cleaned_genes = import_filtered_genes(args.cleaned_gene_annotation, logger)

    # NOTE users will need to edit this as needed!
    # This is the regular expression that is used to extract the chromosome IDs
    # and initialize the DensityData objects. This is specific to the naming of
    # your chromosomes in the HDF5 files.
    # Please consult the docstring of from_list_gene_data_and_hdf5_dir in
    # density_data.py for more information, typically you will only need to
    # edit the part of the string before the underscore. Here this was specific
    # to a naming convention of "DN_(chromosome_id).h5"
    chromosome_string = "DN_(.*?).h5"

    # Get list of GeneData for each genome to enable initialization of
    # DensityData
    genedata_list = get_gene_data_as_list(cleaned_genes)

    # Initialize DensityData for each genome
    # NOTE this object is a list of DensityData instances
    processed_dd_data = DensityData.from_list_gene_data_and_hdf5_dir(
        genedata_list, args.density_data_dir, chromosome_string, logger
    )

    gene_frame_with_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_genes, processed_dd_data
    )

    # NOTE, this adds the indices of the genes in the HDF5 datasets to a pandas
    # dataframe, this is used later on to access the density data for each gene
    gene_frame_with_hdf5_indices = add_hdf5_indices_to_gene_data_from_list_hdf5(
        cleaned_genes, processed_dd_data
    )

    # NOTE, now we can add columns to the pandas dataframe which are the
    # density values for a specific TE type and window combo. This right here
    # adds a column called "LTR_500_Upstream" to the pandas dataframe, this
    # column represents the TE density values for each gene, for the 500 bp
    # upstream window. 500 bp window is used here because I use a smaller
    # testing set of windows for development.

    # NOTE users can loop this over a list of TE types and window sizes if they
    # want
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
