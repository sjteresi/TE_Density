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


def read_arab_expression(arabidopsis_gene_expression):
    """ """
    gene_expression_set = pd.read_csv(
        arabidopsis_gene_expression, sep=",", header="infer"
    )
    gene_expression_set.rename(columns={"Gene": "Gene_Name"}, inplace=True)
    gene_expression_set.set_index("Gene_Name", inplace=True)
    gene_expression_set["Total_Expression_Mean"] = gene_expression_set.mean(
        axis=1
    )  # take mean across
    # row for each gene, i.e average all expression libraries for each gene

    # Drop all columns except the new 'Expression' column
    gene_expression_set.drop(
        gene_expression_set.columns.difference(["Total_Expression_Mean"]),
        1,
        inplace=True,
    )
    gene_expression_set.reset_index(inplace=True)
    return gene_expression_set


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


def gen_violin_plot(gene_frame_within, gene_frame_outside, te_string, output_dir):
    cut_bins, cut_labels = gen_bins()
    fig, axs = plt.subplots(1, 2, figsize=(16, 9), sharey=True)
    fig.subplots_adjust(wspace=0.08)

    # Convert to log 10(x+1) for visualization purposes
    # MAGIC column name
    gene_frame_within["Total_Expression_Mean"] = np.log10(
        gene_frame_within["Total_Expression_Mean"] + 1
    )
    gene_frame_outside["Total_Expression_Mean"] = np.log10(
        gene_frame_outside["Total_Expression_Mean"] + 1
    )

    inside_val_length = len(gene_frame_within)
    outside_val_length = len(gene_frame_outside)

    gene_frame_within["Density_Bins"] = pd.cut(
        gene_frame_within[te_string],
        bins=cut_bins,
        labels=cut_labels,
        include_lowest=True,  # this makes the leftmost bin inclusive
    )
    gene_frame_outside["Density_Bins"] = pd.cut(
        gene_frame_outside[te_string],
        bins=cut_bins,
        labels=cut_labels,
        include_lowest=True,  # this makes the leftmost bin inclusive
    )

    bin_count_frame_inside, bin_count_frame_outside = get_bin_labels(
        gene_frame_within, gene_frame_outside
    )
    # NOTE the labels will be same for both due to the pd.cut syntax being the
    # same
    my_x_tick_labels = bin_count_frame_inside["Bin_IDs"].to_list()

    sns.set_theme()
    sns.violinplot(
        x="Density_Bins",
        y="Total_Expression_Mean",
        data=gene_frame_within,
        ax=axs[1],
        inner="box",
        color="skyblue",
    )
    sns.violinplot(
        x="Density_Bins",
        y="Total_Expression_Mean",
        data=gene_frame_outside,
        ax=axs[0],
        inner="box",
        color="skyblue",
    )

    axs[1].set_xticklabels(my_x_tick_labels, rotation=40, ha="right")
    axs[0].set_xticklabels(my_x_tick_labels, rotation=40, ha="right")
    axs[1].set_title("Inside")
    axs[0].set_title("Outside")
    axs[1].set(xlabel=None, ylabel=None)
    axs[0].set(xlabel=None, ylabel=None)
    axs[1].legend(title=("Total Genes: " + str(inside_val_length)), loc="upper right")
    axs[0].legend(title=("Total Genes: " + str(outside_val_length)), loc="upper right")
    fig.suptitle(
        """Expression Profiles of Genes Inside and Outside the Centromere/Pericentromere by Density Bin"""
    )
    fig.supylabel("Log(Avg+1) Gene Expression in Arabidopsis")
    fig.supxlabel("Density Bins")
    plt.ylim(-1.0, 5.5)
    # plt.show()
    plt.savefig(
        os.path.join(output_dir, "Expression_Profile_Gene_Counts_Centromere.png"),
        bbox_inches="tight",
    )
    plt.close()
    # are using a specific 'library' from the RNA-seq data that represents
    # the stem


def get_bin_labels(inside_values, outside_values):
    finished_bin_count_frames = []
    for pandaframe in [inside_values, outside_values]:
        bin_count_series = pandaframe["Density_Bins"].value_counts()
        bin_count_frame = (
            bin_count_series.to_frame().reset_index()
        )  # Adds the Bin IDs as column, away
        # from index
        bin_count_frame.rename(
            columns={"index": "Bin_IDs", "Density_Bins": "Bin_Counts"},
            inplace=True,
        )
        bin_count_frame.sort_values(
            by=["Bin_IDs"], inplace=True
        )  # sorts the Bin IDs from least
        # to greatest
        finished_bin_count_frames.append(bin_count_frame)

    return finished_bin_count_frames


def gen_barplot_counts(gene_frame_within, gene_frame_outside, te_string, output_dir):
    cut_bins, cut_labels = gen_bins()
    fig, axs = plt.subplots(1, 2, figsize=(12, 8), sharey=True)
    fig.subplots_adjust(wspace=0.08)

    # Get TE values
    inside_values = gene_frame_within[te_string].to_list()
    outside_values = gene_frame_outside[te_string].to_list()
    inside_val_length = len(inside_values)
    outside_val_length = len(outside_values)

    percent_inside_genes_over_50 = ((np.array(inside_values) >= 0.5).sum()) / len(
        inside_values
    )
    percent_outside_genes_over_50 = ((np.array(outside_values) >= 0.5).sum()) / len(
        outside_values
    )
    # NOTE
    print(percent_inside_genes_over_50)
    print(percent_outside_genes_over_50)

    # Reformat as pandas for graphing purposes
    inside_values = pd.DataFrame.from_dict({te_string: inside_values})
    outside_values = pd.DataFrame.from_dict({te_string: outside_values})

    inside_values["Density_Bins"] = pd.cut(
        inside_values[te_string],
        bins=cut_bins,
        labels=cut_labels,
        include_lowest=True,  # this makes the leftmost bin inclusive
    )
    outside_values["Density_Bins"] = pd.cut(
        outside_values[te_string],
        bins=cut_bins,
        labels=cut_labels,
        include_lowest=True,  # this makes the leftmost bin inclusive
    )

    bin_count_frame_inside, bin_count_frame_outside = get_bin_labels(
        inside_values, outside_values
    )

    # NOTE the labels will be same for both due to the pd.cut syntax being the
    # same
    my_x_tick_labels = bin_count_frame_inside["Bin_IDs"].to_list()
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
    plt.savefig(
        os.path.join(output_dir, "Gene_Counts_by_Bin_Centromeric.png"),
        bbox_inches="tight",
    )
    plt.close()


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
        "arabidopsis_gene_expression",
        type=str,
        help="Parent path of arabidopsis gene expression file",
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
    args.arabidopsis_gene_expression = os.path.abspath(args.arabidopsis_gene_expression)
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

    expression_dataset = read_arab_expression(args.arabidopsis_gene_expression)

    # NB MAGIC get unique chromosome ID
    gene_data_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in gene_dataframe_list
    ]

    processed_dd_data = DensityData.from_list_gene_data_and_hdf5_dir(
        gene_data_list, args.density_data_folder, "Arabidopsis_(.*?).h5", logger
    )

    # NOTE
    cleaned_genes.reset_index(inplace=True)  # necessary
    to_concat = []
    for chrom, dataframe in cleaned_genes.groupby(["Chromosome"]):
        for processed_dd_datum in processed_dd_data:
            if processed_dd_datum.unique_chromosome_id == chrom:
                x = processed_dd_datum.add_hdf5_indices_to_gene_data(dataframe)
                to_concat.append(x)
    gene_frame_with_indices = pd.concat(to_concat)

    # NOTE
    to_concat = []
    for chrom, dataframe in gene_frame_with_indices.groupby(["Chromosome"]):
        for processed_dd_datum in processed_dd_data:
            if processed_dd_datum.unique_chromosome_id == chrom:
                x = processed_dd_datum.add_te_vals_to_gene_info_pandas(
                    dataframe, "Order", "LTR", "Upstream", 1000
                )
                to_concat.append(x)
    gene_frame_w_ind_te_vals = pd.concat(to_concat)

    # NOTE
    chromosome_boundary_obj_list = define_boundaries(cleaned_genes)
    to_concat = []
    for chromosome_boundary_obj in chromosome_boundary_obj_list:
        to_assign = gene_frame_w_ind_te_vals.loc[
            gene_frame_w_ind_te_vals["Chromosome"] == chromosome_boundary_obj.chrom
        ].copy(deep=True)
        to_concat.append(
            assign_boundary_ID_to_genes(chromosome_boundary_obj, to_assign)
        )
    genes_w_ind_bndry_te_val = pd.concat(to_concat)

    genes_w_ind_bndry_te_val = pd.merge(
        expression_dataset, genes_w_ind_bndry_te_val, how="inner", on="Gene_Name"
    )

    # NOTE
    gene_frame_within = genes_w_ind_bndry_te_val.loc[
        genes_w_ind_bndry_te_val["Within_Boundary"] == "Y"
    ].copy(deep=True)
    gene_frame_outside = genes_w_ind_bndry_te_val.loc[
        genes_w_ind_bndry_te_val["Within_Boundary"] == "N"
    ].copy(deep=True)

    # NOTE
    te_string = "LTR_1000_Upstream"  # NOTE MAGIC because I hard-coded in the
    # TE info
    gen_barplot_counts(
        gene_frame_within, gene_frame_outside, te_string, args.output_dir
    )
    gen_violin_plot(gene_frame_within, gene_frame_outside, te_string, args.output_dir)
