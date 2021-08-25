#!/usr/bin/env/python

"""
Execute graphing commands
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import pandas as pd
import h5py
import re

from transposon.gene_data import GeneData
from transposon.density_data import DensityData
from transposon.import_filtered_genes import import_filtered_genes


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
    """
    raw_file_list = []  # init empty list to store filenames
    for root, dirs, files in os.walk(path_to_folder):
        for a_file_object in files:
            # N.B very particular usage of abspath and join.
            a_file_object = os.path.abspath(os.path.join(root, a_file_object))
            if a_file_object.endswith(".h5"):  # MAGIC
                raw_file_list.append(a_file_object)

    return raw_file_list


def return_list_density_data_objs(list_of_gene_data, HDF5_folder, file_substring):
    """
    Returns a list of DensityData instances from a list of GeneData and
    HDF5 files that are gathered from an input directory

    Args:
        list_of_gene_data ():
        HDF5_folder ():
        file_substring
    """
    processed_density_data = []
    for raw_hdf5_data_file in supply_density_data_files(HDF5_folder):
        current_hdf5_file_chromosome = re.search(
            file_substring, raw_hdf5_data_file
        ).group(1)
        for gene_data_obj in list_of_gene_data:
            if gene_data_obj.chromosome_unique_id == current_hdf5_file_chromosome:
                dd_data_obj = DensityData.verify_h5_cache(
                    raw_hdf5_data_file, gene_data_obj, logger
                )
                processed_density_data.append(dd_data_obj)
    return processed_density_data


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.abspath(os.path.join(dir_main, "../", "results/graphs"))
    parser = argparse.ArgumentParser(description="generate graphs")

    parser.add_argument(
        "o_sativa_gene_data",
        type=str,
        help="parent path to Arabidopsis' filtered gene data file",
    )

    parser.add_argument(
        "sativa_density_data_dir",
        type=str,
        help="Parent path of folders containing TE Density results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=output_default,
        help="parent directory to output results",
    )
    args = parser.parse_args()
    args.o_sativa_gene_data = os.path.abspath(args.o_sativa_gene_data)
    args.sativa_density_data_dir = os.path.abspath(args.sativa_density_data_dir)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Begin reading files:
    # Get the genes:
    one_gene_data = import_filtered_genes(args.o_sativa_gene_data, logger)

    chromosomes_o_sativa_panda_list = [
        dataframe for k, dataframe in one_gene_data.groupby("Chromosome")
    ]

    gene_data_o_sativa_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in chromosomes_o_sativa_panda_list
    ]

    processed_o_sativa_density_data = return_list_density_data_objs(
        gene_data_o_sativa_list,
        args.sativa_density_data_dir,
        "Sativa_(.*?).h5",
    )

    # TODO begin getting the percentiles for the WHOLE genome not just one
    # chromosome

    test_point = processed_o_sativa_density_data[1]
    test_point_order_left_1000 = test_point.left_orders[:, 1, :]

    for window_idx, window_val in enumerate(test_point.window_list):
        for te_idx, te_grouping in enumerate(test_point.order_list):
            all_upstream_1KB = [
                density_array.left_orders[te_idx, window_idx, :]  # MAGIC
                for density_array in processed_o_sativa_density_data
            ]

            major_concat = np.concatenate((all_upstream_1KB), axis=0)
            # print(major_concat.shape)
            all_chromosomes_cutoff_val = np.percentile(major_concat, 99)
            print(all_chromosomes_cutoff_val)
            # print(major_concat[np.where(major_concat > all_chromosomes_cutoff_val)].shape)

            list_of_genes_with_cutoff = []
            for density_array in processed_o_sativa_density_data:
                # print(
                # np.where(density_array.left_supers[8, 1, :] >= all_chromosomes_cutoff_val)
                # )
                y = np.where(
                    density_array.left_orders[te_idx, window_idx, :]
                    >= all_chromosomes_cutoff_val
                )
                print(density_array.gene_list[y].shape)
                list_of_genes_with_cutoff.append(density_array.gene_list[y])

            with open(
                os.path.join(
                    args.output_dir,
                    str(te_grouping + "_file_" + str(window_val) + ".tsv"),
                ),
                "w",
            ) as f_out:
                for list_of_lists in list_of_genes_with_cutoff:
                    for element in list_of_lists:
                        f_out.write(element + "\n")

    raise ValueError
