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
from examples.Blueberry_Expression.src.compare_expression import gen_bins

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
    arabidopsis_chr_1_boundary = Boundary_Tuple(11.5, 17.8, "Chr1")
    arabidopsis_chr_2_boundary = Boundary_Tuple(1.0, 7.1, "Chr2")
    arabidopsis_chr_3_boundary = Boundary_Tuple(10.2, 17.1, "Chr3")
    arabidopsis_chr_4_boundary = Boundary_Tuple(2.8, 6.2, "Chr4")

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


def return_density_array_from_indices(
    dd_data, te_category, te_name, window_val, direction, indices
):
    density_array = dd_data.get_specific_slice(
        te_category, te_name, window_val, direction, indices
    ).slice
    return density_array


def call_new(
    gene_frames_within_boundary, gene_frames_outside_boundary, processed_dd_data
):
    gene_pandas_inside = pd.concat(gene_frames_within_boundary)
    gene_pandas_outside = pd.concat(gene_frames_outside_boundary)
    gene_names_within = gene_pandas_inside.index.to_list()
    gene_names_outside = gene_pandas_outside.index.to_list()

    # gene_pandas_inside.reset_index(inplace=True)

    #############################
    for chrom, dataframe in gene_pandas_inside.groupby(["Chromosome"]):
        for processed_dd_datum in processed_dd_data:
            if chrom == processed_dd_datum.unique_chromosome_id:
                dataframe["Index_Val"] = dataframe.apply(
                    lambda x: processed_dd_datum._index_of_gene(x["Gene_Name"]), axis=1
                )

                # dataframe["Index_Val"] = dataframe.apply(
                # lambda x: processed_dd_datum._index_of_gene(x["Gene_Name"]), axis=1
                # )
                raise ValueError

    #############################
    raise ValueError
    gene_indices_within = [
        [processed_dd_datum._index_of_gene(gene_name) for gene_name in gene_list]
        for processed_dd_datum, gene_list in zip(processed_dd_data, gene_names_within)
    ]
    gene_indices_outside = [
        [processed_dd_datum._index_of_gene(gene_name) for gene_name in gene_list]
        for processed_dd_datum, gene_list in zip(processed_dd_data, gene_names_outside)
    ]


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
    #########################################
    big_new = DensityData.add_hdf5_indices_to_gene_data(
        processed_dd_data, GeneData(cleaned_genes, "Test")
    )
    print(big_new)
    # big_new.to_csv("test.tsv", sep="\t", index=False, header=True)

    # print(
    # processed_dd_data[0].add_hdf5_indices_to_gene_data(
    # GeneData(cleaned_genes, "Test")
    # )
    # )
    # TODO proceed from here because I now have a pandaframe with all the
    # indices of the genes. Perhaps refactor all of the examples at a later
    # date.
    raise ValueError
    #########################################
    call_new(
        gene_frames_within_boundary, gene_frames_outside_boundary, processed_dd_data
    )
    raise ValueError

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

    dd_and_gene_indices_within = {
        dd_data: gene_indices_within
        for dd_data, gene_indices_within in zip(processed_dd_data, gene_indices_within)
    }
    dd_and_gene_indices_outside = {
        dd_data: gene_indices_outside
        for dd_data, gene_indices_outside in zip(
            processed_dd_data, gene_indices_outside
        )
    }

    cut_bins, cut_labels = gen_bins()

    ##################
    # NOTE I think this should be refactored...
    inside_values = np.concatenate(
        [
            return_density_array_from_indices(
                dd_data, "Order", "LTR", 1000, "Upstream", gene_indices_within
            )
            for dd_data, gene_indices_within in dd_and_gene_indices_within.items()
        ]
    )
    outside_values = np.concatenate(
        [
            return_density_array_from_indices(
                dd_data, "Order", "LTR", 1000, "Upstream", gene_indices_outside
            )
            for dd_data, gene_indices_outside in dd_and_gene_indices_outside.items()
        ]
    )

    #############
    fig, axs = plt.subplots(1, 2, figsize=(12, 8), sharey=True)
    fig.subplots_adjust(hspace=0.5)

    inside_val_length = len(inside_values)
    outside_val_length = len(outside_values)

    inside_values = pd.DataFrame.from_dict({"LTR_1000_Upstream": inside_values})
    outside_values = pd.DataFrame.from_dict({"LTR_1000_Upstream": outside_values})

    inside_values["Density_Bins"] = pd.cut(
        inside_values["LTR_1000_Upstream"],
        bins=cut_bins,
        labels=cut_labels,
        include_lowest=True,  # this makes the leftmost bin inclusive
    )
    outside_values["Density_Bins"] = pd.cut(
        outside_values["LTR_1000_Upstream"],
        bins=cut_bins,
        labels=cut_labels,
        include_lowest=True,  # this makes the leftmost bin inclusive
    )

    bin_count_series_inside = inside_values["Density_Bins"].value_counts()
    bin_count_frame_inside = (
        bin_count_series_inside.to_frame().reset_index()
    )  # Adds the Bin IDs as column, away
    # from index
    bin_count_frame_inside.rename(
        columns={"index": "Bin_IDs", "Density_Bins": "Bin_Counts"},
        inplace=True,
    )
    bin_count_frame_inside.sort_values(
        by=["Bin_IDs"], inplace=True
    )  # sorts the Bin IDs from least
    # to greatest

    my_x_tick_labels = bin_count_frame_inside["Bin_IDs"].to_list()
    #############
    # DO Y

    bin_count_series_outside = outside_values["Density_Bins"].value_counts()
    bin_count_frame_outside = (
        bin_count_series_outside.to_frame().reset_index()
    )  # Adds the Bin IDs as column, away
    # from index
    bin_count_frame_outside.rename(
        columns={"index": "Bin_IDs", "Density_Bins": "Bin_Counts"},
        inplace=True,
    )
    bin_count_frame_outside.sort_values(
        by=["Bin_IDs"], inplace=True
    )  # sorts the Bin IDs from least
    # to greatest

    sns.barplot(
        ax=axs[1],
        x="Bin_IDs",
        y="Bin_Counts",
        color="tab:blue",
        data=bin_count_frame_inside,
    )
    axs[1].set_xticklabels(my_x_tick_labels, rotation=40, ha="right")
    axs[1].set_title("Inside")
    axs[1].set(xlabel=None, ylabel=None)
    axs[1].legend(title=("Total Genes: " + str(inside_val_length)), loc="upper right")
    sns.barplot(
        ax=axs[0],
        x="Bin_IDs",
        y="Bin_Counts",
        color="tab:blue",
        data=bin_count_frame_outside,
    )
    axs[0].set_xticklabels(my_x_tick_labels, rotation=40, ha="right")
    axs[0].set_title("Outside")
    axs[0].set(xlabel=None, ylabel=None)
    axs[0].legend(title=("Total Genes: " + str(outside_val_length)), loc="upper right")
    plt.yscale("log")
    fig.suptitle("Centromeric/Pericentromeric vs 'Regular' Genes")
    fig.supylabel("Log10(No. Genes)")
    fig.supxlabel("Density Bins")
    plt.show()
    plt.savefig("Test_All.png")
