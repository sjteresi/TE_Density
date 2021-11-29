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


def test_chi_square(upstream, downstream, midpoint):
    result = scipy.stats.chisquare([upstream, downstream], f_exp=[midpoint, midpoint])
    return result


def save_arrays_for_pat(processed_dd_data, output_dir):
    # TODO rename function and refactor heavily
    """ """
    ltr_upstream_array = {
        str(
            single_dd_obj.unique_chromosome_id + "_Up"
        ): single_dd_obj.get_specific_slice("Order", "LTR", 1000, "Upstream").slice
        for single_dd_obj in processed_dd_data
    }
    ltr_downstream_array = {
        str(
            single_dd_obj.unique_chromosome_id + "_Down"
        ): single_dd_obj.get_specific_slice("Order", "LTR", 1000, "Downstream").slice
        for single_dd_obj in processed_dd_data
    }

    for (chromosome_up, array_up), (chromosome_down, array_down) in zip(
        ltr_upstream_array.items(), ltr_downstream_array.items()
    ):
        if chromosome_up[:4] == chromosome_down[:4]:
            to_pd = {chromosome_up: array_up, chromosome_down: array_down}
            x = pd.DataFrame.from_dict(to_pd)
            x.to_csv(
                os.path.join(
                    output_dir, str(chromosome_up[:4] + "_LTR_1KB_Up_Down.tsv")
                ),
                sep="\t",
                header=True,
                index=False,
            )

    # NOTE do Helitron now
    helitron_upstream_array = {
        str(
            single_dd_obj.unique_chromosome_id + "_Up"
        ): single_dd_obj.get_specific_slice("Order", "Helitron", 1000, "Upstream").slice
        for single_dd_obj in processed_dd_data
    }
    helitron_downstream_array = {
        str(
            single_dd_obj.unique_chromosome_id + "_Down"
        ): single_dd_obj.get_specific_slice(
            "Order", "Helitron", 1000, "Downstream"
        ).slice
        for single_dd_obj in processed_dd_data
    }

    for (chromosome_up, array_up), (chromosome_down, array_down) in zip(
        helitron_upstream_array.items(), helitron_downstream_array.items()
    ):
        if chromosome_up[:4] == chromosome_down[:4]:
            to_pd = {chromosome_up: array_up, chromosome_down: array_down}
            x = pd.DataFrame.from_dict(to_pd)
            x.to_csv(
                os.path.join(
                    output_dir, str(chromosome_up[:4] + "_Helitron_1KB_Up_Down.tsv")
                ),
                sep="\t",
                header=True,
                index=False,
            )


def get_means_of_dd(processed_dd_data):
    # NOTE
    # This is really dumb and now how you should be doing chi-squared tests
    # which are observed frequencies of categorical data.

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
    midpoint_ltr = (np.array(list_mean_ltr_up) + np.array(list_mean_ltr_down)) / 2

    i = 1
    for upstream, downstream, midpoint in zip(
        list_mean_ltr_up, list_mean_ltr_down, midpoint_ltr
    ):
        print("chromosome #: %s" % i)
        print("AVG LTR upstream val: %s" % upstream)
        print("AVG LTR downstream val: %s" % downstream)
        print("Midpoint LTR val: %s" % midpoint)
        print(test_chi_square(upstream, downstream, midpoint))

        i += 1
        print()

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
    midpoint_helitron = (
        np.array(list_mean_helitron_up) + np.array(list_mean_helitron_down)
    ) / 2

    i = 1
    for upstream, downstream, midpoint in zip(
        list_mean_helitron_up, list_mean_helitron_down, midpoint_helitron
    ):
        print("chromosome #: %s" % i)
        print("AVG Helitron upstream val: %s" % upstream)
        print("AVG Helitron downstream val: %s" % downstream)
        print("Midpoint Helitron val: %s" % midpoint)
        print(test_chi_square(upstream, downstream, midpoint))
        i += 1
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
    # get_means_of_dd(processed_dd_data)
    save_arrays_for_pat(processed_dd_data, args.output_dir)

    # NOTE at this point I have a list of initialized DensityData (chromosome
    # level) and one large expression matrix (all chromosomes).
