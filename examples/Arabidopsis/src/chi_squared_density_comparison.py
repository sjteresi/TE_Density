#!/usr/bin/env python

"""
Compare upstream and downstream TE Density values with a chi-squared test to
see if they are different.
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import numpy as np
import re
import pandas as pd

import scipy.stats

# import matplotlib.pyplot as plt
# import seaborn as sns

# from examples.Arabidopsis.src.import_Arabidopsis_gene_anno import import_genes
from transposon.density_data import DensityData
from transposon.gene_data import GeneData
from transposon.import_filtered_genes import import_filtered_genes


def get_means_of_dd(processed_dd_data):
    list_mean_ltr_up = [
        single_dd_obj.get_specific_slice("Order", "LTR", 1000, "Upstream").slice.mean()
        for single_dd_obj in processed_dd_data
    ]

    list_mean_ltr_down = [
        single_dd_obj.get_specific_slice(
            "Order", "LTR", 1000, "Downstream"
        ).slice.mean()
        for single_dd_obj in processed_dd_data
    ]

    list_mean_helitron_up = [
        single_dd_obj.get_specific_slice(
            "Order", "Helitron", 1000, "Upstream"
        ).slice.mean()
        for single_dd_obj in processed_dd_data
    ]

    list_mean_helitron_down = [
        single_dd_obj.get_specific_slice(
            "Order", "Helitron", 1000, "Downstream"
        ).slice.mean()
        for single_dd_obj in processed_dd_data
    ]

    print(list_mean_ltr_up)
    print(list_mean_ltr_down)
    mean_both_ltr = (np.array(list_mean_ltr_up) + np.array(list_mean_ltr_down)) / 2
    print(mean_both_ltr)

    print()
    print(
        scipy.stats.chisquare(np.array(list_mean_ltr_up), f_exp=np.array(mean_both_ltr))
    )
    print(
        scipy.stats.chisquare(
            np.array(list_mean_ltr_down), f_exp=np.array(mean_both_ltr)
        )
    )
    print()


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
    # TODO delete this
    print()

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

    # print(processed_dd_data[0].windows)
    # print(processed_dd_data[0].window_list)
    # print(processed_dd_data[0].window_index_dict)
    get_means_of_dd(processed_dd_data)

    # NOTE at this point I have a list of initialized DensityData (chromosome
    # level) and one large expression matrix (all chromosomes).
