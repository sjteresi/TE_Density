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
    cleaned_genes = import_filtered_genes(args.o_sativa_gene_data, logger)
    gene_dataframe_list = [
        dataframe for k, dataframe in cleaned_genes.groupby("Chromosome")
    ]
    # NB MAGIC get unique chromosome ID
    gene_data_list = [
        GeneData(dataframe, dataframe["Chromosome"].unique()[0])
        for dataframe in gene_dataframe_list
    ]
    processed_dd_data = DensityData.from_list_gene_data_and_hdf5_dir(
        gene_data_list, args.sativa_density_data_dir, "Sativa_(.*?).h5", logger
    )

    # NOTE
    cleaned_genes.reset_index(inplace=True)  # necessary
    # cleaned_genes = cleaned_genes.loc[cleaned_genes["Chromosome"] == "1"]
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

    upper_cutoff_val = np.percentile(gene_frame_w_ind_te_vals["LTR_1000_Upstream"], 99)
    lower_cutoff_val = np.percentile(gene_frame_w_ind_te_vals["LTR_1000_Upstream"], 1)

    genes_meeting_upper_cutoff = gene_frame_w_ind_te_vals.loc[
        gene_frame_w_ind_te_vals["LTR_1000_Upstream"] >= upper_cutoff_val
    ]["Gene_Name"].to_list()

    upper_cutoff_len = len(genes_meeting_upper_cutoff)

    print(f"Upper Cutoff Length: {upper_cutoff_len}")
    print(f"Upper Cutoff Val: {upper_cutoff_val}")

    print(lower_cutoff_val)
    genes_meeting_lower_cutoff = gene_frame_w_ind_te_vals.loc[
        gene_frame_w_ind_te_vals["LTR_1000_Upstream"] >= lower_cutoff_val
    ]["Gene_Name"].to_list()
    print(len(genes_meeting_lower_cutoff))  # NOTE this is 37741 genes when you
    # faithfully apply the cutoff

    # Take random sample equal to the length of the gene array of upper values.
    # Have to take a random sample because data not normally distributed and
    # the cutoff value for the 1st percentile is so low you would actually get
    # way too many genes if you actually applied it.
    genes_meeting_lower_cutoff = np.random.choice(
        gene_frame_w_ind_te_vals.loc[
            gene_frame_w_ind_te_vals["LTR_1000_Upstream"] >= lower_cutoff_val
        ]["Gene_Name"].to_list(),
        upper_cutoff_len,
        replace=False,
    )

    # NOTE begin writing all the values
    filename_to_write = os.path.join(
        args.output_dir,
        str("Upper_Sample_LTR_1000_Upstream" + ".tsv"),
    )
    with open(
        filename_to_write,
        "w",
    ) as f_out:

        for gene in genes_meeting_upper_cutoff:
            f_out.write(gene + "\n")
    logger.info("Writing to: %s" % filename_to_write)

    filename_to_write = os.path.join(
        args.output_dir,
        str("Lower_Sample_LTR_1000_Upstream" + ".tsv"),
    )
    with open(
        filename_to_write,
        "w",
    ) as f_out:

        for gene in genes_meeting_lower_cutoff:
            f_out.write(gene + "\n")
    logger.info("Writing to: %s" % filename_to_write)
