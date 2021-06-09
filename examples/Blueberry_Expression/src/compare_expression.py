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

from examples.Blueberry_Expression.src.import_blueberry_gene_anno import import_genes
from transposon.density_data import DensityData
from transposon.gene_data import GeneData


def supply_density_data_files(path_to_folder):
    """
    Iterate over a folder containing the H5 files of TE Density output and
    return a list of absolute file paths.

    Args:
        path_to_folder (str): path to the folder containing multiple h5 files
        of density data

    Returns:
        raw_file_list (list of str): A list containing the absolute paths to
        each relevant H5 file of density data

        processed_files (list of str): A list containing absolute paths to each
        previously processed DensityData file
    """
    raw_file_list = []  # init empty list to store filenames
    processed_file_list = []  # init empty list to store filenames
    for root, dirs, files in os.walk(path_to_folder):
        for a_file_object in files:
            # N.B very particular usage of abspath and join.
            a_file_object = os.path.abspath(os.path.join(root, a_file_object))
            if a_file_object.endswith(".h5"):  # MAGIC
                raw_file_list.append(a_file_object)
            if a_file_object.endswith("SenseSwapped.HDF5"):  # MAGIC
                processed_file_list.append(a_file_object)

    return (raw_file_list, processed_file_list)


def verify_sense_swapped_files(
    density_files, processed_density_files, gene_data_list, logger
):
    """
    Return a list of DensityData objects for a given set of TE density (raw
    HDF5) files. Initializes multiple objects

    Args:
        density_files (list of str): A list of str paths to the raw h5 output
        from the TE density tool

        processed_density_files (list of str): A list of str paths to the
        previously processed density data

        gene_data_list (list of str): A list of GeneData instances

        logger (logging.Logger)
    """
    density_data = []
    if len(density_files) == len(processed_density_files):
        for filename in processed_density_files:
            current_hdf5_file_chromosome = re.search(
                "Vacc_Cory_(.*?)_SenseSwapped.HDF5", filename
            ).group(
                1
            )  # MAGIC
            for gene_data_obj in gene_data_list:
                if gene_data_obj.chromosome_unique_id == current_hdf5_file_chromosome:
                    logger.info(
                        "Initializing DensityData from previously processed file %s"
                        % filename
                    )
                    prev_processed_density_data = DensityData(
                        filename, gene_data_obj, logger, sense_swap=False
                    )
                    density_data.append(prev_processed_density_data)
    else:
        for filename in density_files:
            current_hdf5_file_chromosome = re.search(
                "Vacc_Cory_(.*?).h5", filename
            ).group(1)
            # MAGIC
            for gene_data_obj in gene_data_list:
                if gene_data_obj.chromosome_unique_id == current_hdf5_file_chromosome:
                    logger.info("Initializing DensityData from raw file %s" % filename)
                    prev_processed_density_data = DensityData(
                        filename, gene_data_obj, logger, sense_swap=False
                    )
                    density_data.append(prev_processed_density_data)
    return density_data


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
    matrix = pd.read_csv(tpm_matrix_file, header="infer", sep="\t")
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

    ########################
    # TODO could probably refactor the above and below sections into their own
    # places so that it is more coherent. The above portion has to do with
    # getting a fully complete dataframe and the below portion has to do with
    # actually plotting it all.
    ########################
    step = 0.1
    cut_bins = np.arange(0, 1.1, step)  # MAGIC get bins for plotting
    cut_labels = [
        "(" + str(round(x, 2)) + " " + str(round(x + step, 2)) + "]" for x in cut_bins
    ]  # this is ordered from least to greatest
    cut_labels.pop()
    cut_labels[0] = cut_labels[0].replace("(", "[")  # MAGIC to get the
    # correct label for the first bin
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
            color="skyblue",
        )

        plt.legend(title=("Total Genes: " + str(total_genes)), loc="upper right")
        plt.ylabel(str(blueberry_ex_column + " log(TPM)"))
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

        title_info = data_column.split("_")
        # NB type, window, direction
        title_string = str(
            title_info[0] + " " + str(title_info[1]) + " BP " + title_info[2]
        )
        plt.title(title_string)

        figure = plt.gcf()
        figure.set_size_inches(10.5, 8)  # MAGIC get figure size to fit
        # everything nicely
        plt.xlabel("Density Bins")
        logger.info("Graph created for %s" % data_column)
        plt.savefig(
            os.path.join(args.output_dir, str(data_column + ".png")),
            bbox_inches="tight",
        )
        if show:
            plt.show()
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

    # MAGIC tuple for files
    density_files = supply_density_data_files(args.density_data_folder)[0]
    processed_density_files = supply_density_data_files(args.density_data_folder)[1]

    processed_dd_data = verify_sense_swapped_files(
        density_files, processed_density_files, gene_data_list, logger
    )

    # NOTE at this point I have a list of initialized DensityData (chromosome
    # level) and one large expression matrix (all chromosomes). But the
    # expression matrix doesn't actually have the chromosome information so I
    # can't easily split by that method. So I need to use the list of genes in
    # the density data to split the 'chromosomes' of the expression matrix.

    blueberry_ex_column = "Dra_4dpo_c1_S46_L005"  # MAGIC column for a single

    plot_expression_v_density_violin(
        tpm_matrix,
        processed_dd_data,
        args.output_dir,
        logger,
        blueberry_ex_column,
        show=False,
    )