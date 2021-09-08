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

    processed_o_sativa_density_data = DensityData.from_list_gene_data_and_hdf5_dir(
        gene_data_o_sativa_list, args.sativa_density_data_dir, "Sativa_(.*?).h5", logger
    )

    # TODO begin getting the percentiles for the WHOLE genome not just one
    # chromosome

    # TODO refactor because this is messy, tried to make as little hard-coded
    # as possible.
    # Essentially iterates over the densitydata, grabs all the arrays for the
    # 1KB upstream data for the Order set, and calculates the 99th percentile
    # for each TE type for all genes, then outputs a txt file with the genes
    # meeting that cutoff as individual rows in the output file

    # MAGIC, just get an iterable, not actually operating on a single
    # DensityData
    test_point = processed_o_sativa_density_data[1]
    test_point_order_left_1000 = test_point.left_orders[:, 1, :]

    for window_idx, window_val in enumerate(test_point.window_list):
        for te_idx, te_grouping in enumerate(test_point.order_list):
            all_upstream_1KB = [
                density_array.left_orders[te_idx, window_idx, :]  # MAGIC
                for density_array in processed_o_sativa_density_data
            ]

            major_concat = np.concatenate(all_upstream_1KB, axis=0)
            all_chromosomes_cutoff_val = np.percentile(major_concat, 99)
            # MAGIC 99th percentle cutoff

            list_of_genes_with_cutoff = []
            for density_array in processed_o_sativa_density_data:
                y = np.where(
                    density_array.left_orders[te_idx, window_idx, :]
                    >= all_chromosomes_cutoff_val
                )
                list_of_genes_with_cutoff.append(density_array.gene_list[y].tolist())

            filename_to_write = os.path.join(
                args.output_dir,
                str(te_grouping + "_file_" + str(window_val) + ".tsv"),
            )
            with open(
                filename_to_write,
                "w",
            ) as f_out:

                for list_of_lists in list_of_genes_with_cutoff:
                    for element in list_of_lists:
                        f_out.write(element + "\n")
                logger.info("Writing to: %s" % filename_to_write)
