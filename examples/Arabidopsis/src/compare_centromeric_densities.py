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
from collections import namedtuple


import matplotlib.pyplot as plt
import seaborn as sns

from transposon.density_data import DensityData
from transposon.gene_data import GeneData
from transposon.import_filtered_genes import import_filtered_genes

Boundary_Tuple = namedtuple("Boundary_Tuple", ["start", "stop", "chrom"])


def define_boundaries(cleaned_gene):
    """
    Args:
        cleaned_genes (pandas.Data.Frame):

    Returns:

    """
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

    # MAGIC centromeric/pericentromeric boundaries derived from publication
    arabidopsis_chr_1_boundary = Boundary_Tuple(10.5, 17.8, "Chr1")
    arabidopsis_chr_2_boundary = Boundary_Tuple(1.0, 6.0, "Chr2")
    arabidopsis_chr_3_boundary = Boundary_Tuple(10.2, 17.1, "Chr3")
    arabidopsis_chr_4_boundary = Boundary_Tuple(
        1.0, 2.0, "Chr4"
    )  # TODO verify with Pat,
    # temporary place holder
    arabidopsis_chr_5_boundary = Boundary_Tuple(9.0, 16.0, "Chr5")
    chromosome_boundary_obj_list = [
        arabidopsis_chr_1_boundary,
        arabidopsis_chr_2_boundary,
        arabidopsis_chr_3_boundary,
        arabidopsis_chr_4_boundary,
        arabidopsis_chr_5_boundary,
    ]

    # NB convert to BP
    chromosome_boundary_obj_list = [
        Boundary_Tuple(start * 1000000, stop * 1000000, chrom)
        for start, stop, chrom in chromosome_boundary_obj_list
    ]
    return chromosome_boundary_obj_list


def get_genes_within_boundary(chromosome_boundary_obj, cleaned_genes):
    """

    Returns:
        (pandaframe)

    """
    cleaned_genes_subsetted = cleaned_genes.loc[
        cleaned_genes["Chromosome"] == chromosome_boundary_obj.chrom
    ]
    genes_within = cleaned_genes_subsetted.loc[
        (cleaned_genes_subsetted["Start"] >= chromosome_boundary_obj.start)
        & (cleaned_genes_subsetted["Stop"] <= chromosome_boundary_obj.stop)
    ].copy(deep=True)
    return genes_within


def get_genes_outside_boundary(chromosome_boundary_obj, cleaned_genes):
    """

    Returns:
        (pandaframe)

    """
    cleaned_genes_subsetted = cleaned_genes.loc[
        cleaned_genes["Chromosome"] == chromosome_boundary_obj.chrom
    ]
    genes_outside = cleaned_genes_subsetted.loc[
        (cleaned_genes_subsetted["Stop"] <= chromosome_boundary_obj.start)
        | (cleaned_genes_subsetted["Start"] >= chromosome_boundary_obj.stop)
    ].copy(deep=True)
    return genes_outside


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

    # NB constructed in ascending chromosome order
    chromosome_boundary_obj_list = define_boundaries(cleaned_genes)
    gene_frames_within_boundary = [
        get_genes_within_boundary(chromosome_boundary_obj, cleaned_genes)
        for chromosome_boundary_obj in chromosome_boundary_obj_list
    ]
    gene_frames_outside_boundary = [
        get_genes_outside_boundary(chromosome_boundary_obj, cleaned_genes)
        for chromosome_boundary_obj in chromosome_boundary_obj_list
    ]

    gene_names_within = [
        gene_frame.index.to_list() for gene_frame in gene_frames_within_boundary
    ]
    gene_names_outside = [
        gene_frame.index.to_list() for gene_frame in gene_frames_outside_boundary
    ]

    # Get their indices
    # TODO consider refactoring to function because the gene list is the only
    # thing changing
    gene_indices_within = [
        [processed_dd_datum._index_of_gene(gene_name) for gene_name in gene_list]
        for processed_dd_datum, gene_list in zip(processed_dd_data, gene_names_within)
    ]
    gene_indices_outside = [
        [processed_dd_datum._index_of_gene(gene_name) for gene_name in gene_list]
        for processed_dd_datum, gene_list in zip(processed_dd_data, gene_names_outside)
    ]
    # gene_indices_within.sort()

    # Extract their TE values and then compare the distros
    # x = processed_dd_data[0].left_supers[:, :, [1, 2, 3]]
    # y = processed_dd_data[0].left_supers[:, :, slice(None)]
    # print(x)
    # print(y)
    # print(np.all((x == y)))
    x = (
        processed_dd_data[0]
        .get_specific_slice("Order", "LTR", 1000, "Upstream", gene_indices_within[0])
        .slice
    )
    y = (
        processed_dd_data[0]
        .get_specific_slice("Order", "LTR", 1000, "Upstream", gene_indices_outside[0])
        .slice
    )
    # print(x)
    print(len(x))
    print(x.shape)

    # print()
    # print(len(gene_names_within[0]))
    # print(len(x))

    sns.displot(x, kind="hist")
    # TODO bin the stuff before we do the plot
    plt.yscale("log")
    plt.title("Inside")
    plt.savefig("Test_Inside.png")
    plt.gcf()
    sns.displot(y, kind="hist")
    # TODO bin the stuff before we do the plot
    # TODO create a subplot so that y-axes can be shared.
    plt.yscale("log")
    plt.title("Outside")
    plt.savefig("Test_Outside.png")

    # plt.show()
