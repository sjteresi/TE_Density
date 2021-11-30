#!/usr/bin/env python

"""
Compare TE density values with gene expression values for the blueberry system
as an usage example for the tool
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels

from examples.Blueberry_Expression.src.import_blueberry_gene_anno import import_genes
from transposon.density_data import DensityData
from transposon.gene_data import GeneData


def read_TPM_matrix(tpm_matrix_file):
    """
    Read a tpm matrix into a pandas object

    Args:
        tpm_matrix_file (str): path to tpm matrix file, shape (genes x
        libraries)

    Returns:
        matrix (pandas.DataFrame, shape(genes x libraries)): genes are
        technically the first column
    """
    matrix = pd.read_csv(
        tpm_matrix_file, header="infer", sep="\t", index_col="Gene_Name"
    )
    matrix["Total_Expression_Mean"] = matrix.mean(axis=1)
    matrix.drop(
        matrix.columns.difference(["Total_Expression_Mean"]), axis=1, inplace=True
    )
    matrix.reset_index(inplace=True)
    return matrix


def verify_gene_order(tpm_matrix, density_data):
    """
    Reorder the TPM matrix's genes so that they correspond 1:1 to the output of
    the TE Density tool.

    Args:
        tpm_matrix (pandas.DataFrame): A pandas dataframe of the gene
        expression values in TPM format. Shape (genes X rows).

        density_data (DensityData): DensityData object containing the TE
        Density values for ONE chromosome.

    Returns:
        tpm_matrix (pandas.DataFrame, shape(genes x libraries)): genes are now
        index values, and are re-ordered to match the order prvided in the
        density data so that we may later extract values easily for comparison.
    """
    tpm_matrix = tpm_matrix.reindex(density_data.gene_list)
    if tpm_matrix.index.to_list() != density_data.gene_list.tolist():
        raise ValueError(
            """Was unable to reorder the TPM matrix's genes in the
                         same order as the density data. Will produce erroneous
                         results when graphing because each gene's information
                         will not be in the appropriate order for
                         comparison."""
        )
    return tpm_matrix


def add_TE_density_to_expression_matrix(
    tpm_matrix, processed_dd_data, output_dir, logger, show=False
):

    one_chrom_storage = []
    for dd_chromosome in processed_dd_data:
        # NOTE the regular tpm matrix represents all genes, so we must break it
        # apart into chunks to get the correct density data before we merge it
        # back together
        subsetted_tpm_matrix_by_genes = tpm_matrix[
            tpm_matrix["Gene_Name"].isin(dd_chromosome.gene_list)
        ]

        if len(subsetted_tpm_matrix_by_genes) != len(dd_chromosome.gene_list):
            raise ValueError(
                """Unable to subset the TPM matrix correctly so that
                             its genes directly correspond to what is in the
                             DensityData object."""
            )
        subsetted_tpm_matrix_by_genes.set_index(
            "Gene_Name", inplace=True
        )  # set index as gene names
        subsetted_tpm_matrix_by_genes = subsetted_tpm_matrix_by_genes.reindex(
            dd_chromosome.gene_list
        )

        te_columns_to_graph = []  # NB container to store string names, doesn't
        # matter that it is overwritten between loops, the names will be the
        # same for the last loop

        for density_slice in dd_chromosome.yield_all_slices():
            # Here we are still on a single chromosome, beginning to loop over
            # all of the windows, TEs, and directions

            # NOTE yield all slices is a function that skips the Revision
            # datasets. TODO this will likely be accomplished through other
            # means in future releases.
            window_direction_type = str(
                density_slice.te_type
                + "_"
                + str(density_slice.window_val)
                + "_"
                + density_slice.direction
            )

            subsetted_tpm_matrix_by_genes[window_direction_type] = density_slice.slice

            # Store the TE column names to graph later
            te_columns_to_graph.append(window_direction_type)

        # Now I have a matrix that contains all of the TE Data for one
        # chromosome, so we will save it for concatenating
        one_chrom_storage.append(subsetted_tpm_matrix_by_genes)

    complete_df = pd.concat(one_chrom_storage)
    return (complete_df, te_columns_to_graph)


def gen_bins():
    # TODO
    """
    Return
    """
    step = 0.1
    cut_bins = np.arange(0, 1.1, step)  # MAGIC get bins for plotting
    cut_labels = [
        "(" + str(round(x, 2)) + " " + str(round(x + step, 2)) + "]" for x in cut_bins
    ]  # this is ordered from least to greatest
    cut_labels.pop()
    cut_labels[0] = cut_labels[0].replace("(", "[")  # MAGIC to get the
    # correct label for the first bin
    return (cut_bins, cut_labels)


def plot_expression_v_density_violin(
    tpm_matrix,
    processed_dd_data,
    output_dir,
    logger,
    blueberry_ex_column,
    show=False,
):
    """
    Create a violin plot of density values for each gene vs the expression
    values for each gene

    Args:
        tpm_matrix (pandas.DataFrame): A pandas dataframe of the gene
        expression values in TPM format. Shape (genes X 1 column of
        expression values

        processed_dd_data (list of DensityData): DensityData object containing the TE
        Density values for ONE chromosome. TODO edit

        output_dir (str): path of folder to output results

        logger (logging.Logger): logger object

        blueberry_ex_column (str):


    Returns:
        None, creates graph, saves it to disk, and shows it if the keyword
        argument show is provided as True.
    """
    complete_df, te_columns_to_graph = add_TE_density_to_expression_matrix(
        tpm_matrix,
        processed_dd_data,
        output_dir,
        logger,
        blueberry_ex_column,
        show=False,
    )

    ########################
    # TODO could probably refactor the above and below sections into their own
    # places so that it is more coherent. The above portion has to do with
    # getting a fully complete dataframe and the below portion has to do with
    # actually plotting it all.
    ########################
    cut_bins, cut_labels = gen_bins

    complete_df = complete_df[complete_df[blueberry_ex_column] >= 0.1]
    # MAGIC we only want to plot genes with expression values over 0.1
    # TPM
    total_genes = len(complete_df)

    # Convert to log 10 for visualization purposes
    # MAGIC
    complete_df[blueberry_ex_column] = np.log10(complete_df[blueberry_ex_column])

    for data_column in te_columns_to_graph:
        complete_df["Density_Bins"] = pd.cut(
            complete_df[data_column],
            bins=cut_bins,
            labels=cut_labels,
            include_lowest=True,  # this makes the leftmost bin inclusive
        )

        sns.set_theme()
        sns.violinplot(
            x="Density_Bins",
            y=blueberry_ex_column,
            data=complete_df,
            inner="box",
            color="skyblue",
        )

        plt.legend(title=("Total Genes: " + str(total_genes)), loc="upper right")
        plt.ylabel("Gene Expression in Blueberry Stem log(TPM)")  # MAGIC we
        # are using a specific 'library' from the RNA-seq data that represents
        # the stem
        bin_count_series = complete_df["Density_Bins"].value_counts()
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
        bin_ids = bin_count_frame["Bin_IDs"].tolist()
        if bin_ids != cut_labels:
            raise ValueError(
                "Bin IDs not equal to cut_labels. %s != %s" % (bin_ids, cut_labels)
            )
        bin_counts = bin_count_frame["Bin_Counts"].tolist()

        new_labels = []
        for individual_bin, individual_bin_count in zip(bin_ids, bin_counts):
            new_labels.append(
                individual_bin + "\n" + " N= " + str(individual_bin_count)
            )
        normal_locs, _ = plt.xticks()
        plt.xticks(normal_locs, new_labels)
        plt.tight_layout()

        # MAGIC
        title_info = data_column.split("_")
        direction = data_column.split("_")[-1]
        bp_val = data_column.split("_")[-2]
        te_type = data_column.split("_")[0:-2]
        te_type = " ".join(te_type)

        # NB type, window, direction
        title_string = str(te_type + " " + bp_val + " BP " + direction)
        plt.title(title_string)

        figure = plt.gcf()
        figure.set_size_inches(10.5, 8)  # MAGIC get figure size to fit
        # everything nicely
        plt.xlabel("Density Bins")
        plt.ylim(-2.0, 5.5)
        logger.info("Violin plot created for %s" % data_column)
        plt.savefig(
            os.path.join(output_dir, str(data_column + "_ViolinPlot.png")),
            bbox_inches="tight",
        )
        if show:
            plt.show()

        plt.clf()

        # Call the code to plot the line plot of gene counts
        # NOTE
        # This is kinda crappy because I am calling another function to graph
        # another graph here. Ideally it would be refactored in a later
        # release.
        # plot_expression_v_density_gene_counts(
        # bin_count_frame, title_string, data_column, output_dir, logger
        # )


def plot_new_hist(
    tpm_matrix,
    processed_dd_data,
    output_dir,
    logger,
    blueberry_ex_column,
    show=False,
):
    """
    Generate a multi-panel plot of histograms.
    Histograms will be the number of genes in a given density bin.
    Subplot one will be the histogram of all genes.
    Subplot two will be the histogram of low/non-expressed genes.
    Subplot three will be the histogram of expressed genes.


    """
    blueberry_ex_column = "Total_Expression_Mean"
    complete_df, te_columns_to_graph = add_TE_density_to_expression_matrix(
        tpm_matrix,
        processed_dd_data,
        output_dir,
        logger,
        show=False,
    )
    cut_bins, cut_labels = gen_bins()

    expressed_genes = complete_df.loc[complete_df[blueberry_ex_column] >= 0.1].copy(
        deep=True
    )
    non_expressed_genes = complete_df.loc[complete_df[blueberry_ex_column] <= 0.1].copy(
        deep=True
    )
    all_genes = complete_df.copy(deep=True)

    # Make list object to iterate over to work on
    panda_dataframes = {
        "All Genes": all_genes,
        "Low/Non-Expressed Genes": non_expressed_genes,
        "Expressed Genes": expressed_genes,
    }

    fig, axs = plt.subplots(1, 3, figsize=(12, 8), sharey=True)
    fig.subplots_adjust(hspace=0.5)
    # axs = axs.ravel()  # NOTE what does this do?
    i = 0
    for name, pandaframe in panda_dataframes.items():
        # Convert to expression values to log 10 for visualization purposes
        # Perform  np.log10(x+1) on the expression column, and 1 so we can do
        # log 10 on zero values.
        pandaframe[blueberry_ex_column] = pandaframe[blueberry_ex_column].map(
            lambda x: np.log10(x + 1)
        )

        # Determine the bins and number of genes in each
        pandaframe["Density_Bins"] = pd.cut(
            pandaframe["LTR_1000_Upstream"],
            bins=cut_bins,
            labels=cut_labels,
            include_lowest=True,  # this makes the leftmost bin inclusive
        )

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

        # Do I need to log on this?
        bin_count_frame["Bin_Counts"] = bin_count_frame["Bin_Counts"].map(
            lambda x: np.log10(x)
        )
        my_x_tick_labels = bin_count_frame["Bin_IDs"].to_list()
        sns.barplot(
            ax=axs[i],
            x="Bin_IDs",
            y="Bin_Counts",
            color="tab:blue",
            data=bin_count_frame,
        )

        # Failed regression plot, didn't look very good
        # sns.regplot(
        #    ax=axs[i],
        #    x=np.arange(0, len(my_x_tick_labels)),
        #    y="Bin_Counts",
        #    color="tab:red",
        # data=bin_count_frame,
        # )

        # TODO get A B C subplot labels?
        axs[i].set_xticklabels(my_x_tick_labels, rotation=40, ha="right")
        axs[i].set_title(name)
        axs[i].set(xlabel=None, ylabel=None)
        axs[i].legend(title=("Total Genes: " + str(len(pandaframe))), loc="upper right")

        i += 1

    fig.supylabel("Log10(No. Genes)")
    fig.supxlabel("Density Bins")
    fig.suptitle(
        "No. Genes Binned by LTR TE 1KB Upstream Density According to Expression Profile"
    )
    # TODO parametrize the file name or put up a MAGIC notifier
    plt.savefig(os.path.join(output_dir, "Test.png"))
    plt.show()


def plot_expression_v_density_gene_counts(
    bin_count_frame, title_string, data_column, output_dir, logger
):
    """
    Create a histogram plot with the density bins on the x-axis and the number of
    genes in that given bin on the y-axis.

    Args:
        bin_count_frame ():
        title_string ():
        data_column ():
        output_dir ():
        logger ():

    Returns:
        None, generates a graph and saves it to disk
    """
    # MAGIC code below for graphing and getting correct columns
    # NOTE the bin counts here still represent genes that were pre-filtered
    # (i.e having expression greater than "lowly expressed")
    bin_count_frame["Bin_Counts"] = np.log10(bin_count_frame["Bin_Counts"] + 1)

    sns.barplot(x="Bin_IDs", y="Bin_Counts", color="tab:blue", data=bin_count_frame)

    # sns.lineplot(x="Bin_IDs", y="Bin_Counts", marker="o", data=bin_count_frame)
    plt.title(title_string)
    plt.ylabel("Log(x+1) Number of Genes in Density Bin")
    plt.xlabel("Density Bins")
    plt.yticks(np.arange(0.0, 4.75, 0.25))  # MAGIC set of values for y axis
    plt.ylim(0, 4.5)  # MAGIC set ylim
    plt.tight_layout()
    logger.info("Bar plot created for %s" % data_column)
    plt.savefig(
        os.path.join(output_dir, str(data_column + "_GeneCounts.png")),
        bbox_inches="tight",
    )
    plt.clf()


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.join(dir_main, "..", "results/graphs")
    parser = argparse.ArgumentParser(
        description="compare density values with gene expression values"
    )
    parser.add_argument(
        "TPM_matrix", type=str, help="parent path of TPM matrix for all genes"
    )

    parser.add_argument(
        "gene_input_file", type=str, help="Parent path of gene annotation file"
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
    args.TPM_matrix = os.path.abspath(args.TPM_matrix)
    args.gene_input_file = os.path.abspath(args.gene_input_file)
    args.output_dir = os.path.abspath(args.output_dir)
    args.density_data_folder = os.path.abspath(args.density_data_folder)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    tpm_matrix = read_TPM_matrix(args.TPM_matrix)
    cleaned_genes = import_genes(args.gene_input_file, logger)
    gene_dataframe_list = [
        dataframe for k, dataframe in cleaned_genes.groupby("Chromosome")
    ]
    gene_data_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in gene_dataframe_list
    ]

    processed_dd_data = DensityData.from_list_gene_data_and_hdf5_dir(
        gene_data_list, args.density_data_folder, "Vacc_Cory_(.*?).h5", logger
    )

    # NOTE at this point I have a list of initialized DensityData (chromosome
    # level) and one large expression matrix (all chromosomes). But the
    # expression matrix doesn't actually have the chromosome information so I
    # can't easily split by that method. So I need to use the list of genes in
    # the density data to split the 'chromosomes' of the expression matrix.

    blueberry_ex_column = "Dra_4dpo_c1_S46_L005"  # MAGIC column for a single

    # plot_expression_v_density_violin(
    # tpm_matrix,
    # processed_dd_data,
    # args.output_dir,
    # logger,
    # blueberry_ex_column,
    # show=False,
    # )

    plot_new_hist(
        tpm_matrix,
        processed_dd_data,
        args.output_dir,
        logger,
        blueberry_ex_column,
        show=False,
    )
