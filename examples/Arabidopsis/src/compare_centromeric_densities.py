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

    Returns: chromosome_boundary_obj_list (list of Boundary_Tuple)

    """
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


def assign_boundary_ID_to_genes(chromosome_boundary_obj, cleaned_genes):
    """
    NOTE if anything overlaps with the edges I put it in the
    non-centromeric/pericentromeric category

    Returns:
        (pandaframe)

    """
    cleaned_genes["Within_Boundary"] = cleaned_genes.apply(
        lambda x: "Y"
        if (x["Start"] >= chromosome_boundary_obj.start)
        & (x["Stop"] <= chromosome_boundary_obj.stop)
        else "N",
        axis=1,
    )
    return cleaned_genes


def return_1D_array_of_vals_for_indices(
    gene_frame,
    list_of_density_data,
    te_category_str="Order",
    te_name_str="LTR",
    window_val=1000,
    direction_str="Upstream",
):
    to_concat = []
    for chrom, dataframe in gene_frame.groupby(["Chromosome"]):
        for dd_obj in processed_dd_data:
            if chrom == dd_obj.unique_chromosome_id:
                indices = dataframe["Index_Val"].to_list()
                te_values = dd_obj.get_specific_slice(
                    te_category_str,
                    te_name_str,
                    window_val,
                    direction_str,
                    dataframe["Index_Val"].to_list(),
                ).slice
                to_concat.append(te_values)
    return np.concatenate(to_concat)


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
    #####################
    # NOTE begin analysis

    cleaned_genes = import_filtered_genes(args.gene_input_file, logger)
    gene_dataframe_list = [
        dataframe for k, dataframe in cleaned_genes.groupby("Chromosome")
    ]

    # NB MAGIC get unique chromosome ID
    gene_data_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in gene_dataframe_list
    ]

    processed_dd_data = DensityData.from_list_gene_data_and_hdf5_dir(
        gene_data_list, args.density_data_folder, "Arabidopsis_(.*?).h5", logger
    )

    # NB add the HDF5 indices to the pandaframe
    # MAGIC dummy genome name
    gene_frame_with_indices = DensityData.add_hdf5_indices_to_gene_data(
        processed_dd_data, GeneData(cleaned_genes, "dummy")
    )

    # NB constructed in ascending chromosome order
    chromosome_boundary_obj_list = define_boundaries(cleaned_genes)

    to_concat = []
    for chromosome_boundary_obj in chromosome_boundary_obj_list:
        to_assign = gene_frame_with_indices.loc[
            gene_frame_with_indices["Chromosome"] == chromosome_boundary_obj.chrom
        ].copy(deep=True)
        to_concat.append(
            assign_boundary_ID_to_genes(chromosome_boundary_obj, to_assign)
        )
    genes_w_ind_and_bndry = pd.concat(to_concat)

    gene_frame_within = genes_w_ind_and_bndry.loc[
        genes_w_ind_and_bndry["Within_Boundary"] == "Y"
    ].copy(deep=True)
    gene_frame_outside = genes_w_ind_and_bndry.loc[
        genes_w_ind_and_bndry["Within_Boundary"] == "N"
    ].copy(deep=True)

    inside_values = return_1D_array_of_vals_for_indices(
        gene_frame_within, processed_dd_data
    )
    outside_values = return_1D_array_of_vals_for_indices(
        gene_frame_outside, processed_dd_data
    )

    cut_bins, cut_labels = gen_bins()

    #############
    fig, axs = plt.subplots(1, 2, figsize=(12, 8), sharey=True)
    fig.subplots_adjust(hspace=0.5)

    inside_val_length = len(inside_values)
    outside_val_length = len(outside_values)

    percent_inside_genes_over_50 = ((inside_values >= 0.5).sum()) / len(inside_values)
    percent_outside_genes_over_50 = ((outside_values >= 0.5).sum()) / len(
        outside_values
    )
    # NOTE
    print(percent_inside_genes_over_50)
    print(percent_outside_genes_over_50)

    # Reformat as pandas for graphing purposes
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
    fig.suptitle(
        """Counts of Genes Inside and Outside the Centromere/Pericentromere by Density Bin"""
    )
    fig.supylabel("Log10(No. Genes)")
    fig.supxlabel("Density Bins")
    # plt.show()
    plt.savefig("Gene_Counts_by_Bin_Centromeric.png")
