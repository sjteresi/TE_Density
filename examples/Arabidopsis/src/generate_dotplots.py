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
from collections import defaultdict

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

from transposon.gene_data import GeneData
from transposon.density_data import DensityData
from transposon.import_filtered_genes import import_filtered_genes


def plot_intra_density(dd_obj, order_or_super, output_dir, display=False):
    """
    Args:
        dd_obj (DensityData):

        order_or_super (str): String representing the category of TE, should be either
        'Order' or 'Superfamily', raises ValueError if it does not match.

        output_dir (str): Path representing output directory.

        display (str): Defaults to False, if True it will show the graph
        (usually opens another window).

    """
    # NOTE intra code has not been edited for multiple chromosomes yet
    # NOTE unused, example code
    plotting_dict = {}

    if order_or_super == "Order":
        te_index_dict = dd_obj.order_index_dict
        h5_frame = dd_obj.intra_orders
    elif order_or_super == "Superfamily":
        te_index_dict = dd_obj.super_index_dict
        h5_frame = dd_obj.intra_supers
    else:
        raise ValueError("Please provide Order or Superfamily")
    # NOTE
    # I could add a FOR loop here to go over each dd_obj and get the mean that
    # way, plotting that would be super easy and would only yield one line.
    for te_type, index_val in te_index_dict.items():
        plotting_dict[te_type] = np.mean(h5_frame[te_index_dict[te_type], :, :])

    for key, val in plotting_dict.items():
        plt.scatter(0, val, label=key)
        plt.legend()

    plt.yticks(np.arange(0, 0.8, 0.05))  # NOTE I want this to be consistent
    # but it is too high, can't do std dev? Needs to be consistent with the
    # multiplots
    plt.xticks([0])
    plt.xlabel("Window Position in BP")
    plt.ylabel("Intragenic TE Density")

    if order_or_super == "Order":
        plt.title("TE Orders")
    elif order_or_super == "Superfamily":
        plt.title("TE Superfamilies")
    else:
        raise ValueError("Please provide Order or Superfamily")

    plt.savefig(
        os.path.join(
            output_dir, str(order_or_super + "_" + dd_obj.genome_id + "_Intra_Plot.png")
        )
    )
    if display:
        plt.show()
    plt.close()


def plot_density_all(dd_obj, order_or_super, output_dir, logger, display=False):
    """
    Plot the mean density value of each TE category over each window for "left"
    and "right".

    Args:
        dd_obj (DensityData):

        order_or_super (str): String representing the category of TE, should be either
        'Order' or 'Superfamily', raises ValueError if it does not match.

        output_dir (str): Path representing output directory.

        logger (logging.Logger):

        display (str): Defaults to False, if True it will show the graph
        (usually opens another window).

    """
    if order_or_super == "Order":
        te_index_dict = dd_obj.order_index_dict
        left_h5_frame = dd_obj.left_orders
        intra_h5_frame = dd_obj.intra_orders
        right_h5_frame = dd_obj.right_orders
    elif order_or_super == "Superfamily":
        te_index_dict = dd_obj.super_index_dict
        left_h5_frame = dd_obj.left_supers
        intra_h5_frame = dd_obj.intra_supers
        right_h5_frame = dd_obj.right_supers
    else:
        raise ValueError("Please provide Order or Superfamily")

    # NOTE rename this to upstream and downstream
    data_left_dict = defaultdict(list)
    data_intra_dict = defaultdict(list)
    data_right_dict = defaultdict(list)
    for te_type, index_val in te_index_dict.items():

        if "Revision" in te_type:
            continue  # NB I don't want to graph the revision data, see
        # revision documentation as to why
        for window_idx in range(len(dd_obj.window_list)):

            data_left_dict[te_type].append(
                np.mean(left_h5_frame[te_index_dict[te_type], window_idx, :])
            )
            data_right_dict[te_type].append(
                np.mean(right_h5_frame[te_index_dict[te_type], window_idx, :])
            )
        data_intra_dict[te_type].append(
            np.mean(intra_h5_frame[te_index_dict[te_type], :, :])
        )

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey="col")
    fig.set_size_inches(16, 9.5)

    # define colors
    NUM_COLORS = sum(1 for te_type in te_index_dict.items())
    cm = plt.get_cmap("tab20")
    ax1.set_prop_cycle("color", [cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])
    ax2.set_prop_cycle("color", [cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])
    ax3.set_prop_cycle("color", [cm(1.0 * i / NUM_COLORS) for i in range(NUM_COLORS)])

    for key, val in data_left_dict.items():
        ax1.plot(
            dd_obj.window_list,
            val,
            label=key,
            linestyle=(0, (3, 1, 1, 1)),
            marker="o",
        )
    ax1.set(
        xlabel="BP Upstream",
        ylabel="TE Density",
        xlim=[max(dd_obj.window_list), min(dd_obj.window_list)],
        # yticks=np.arange(0, 0.51, 0.025), # MAGIC
        xticks=range(min(dd_obj.window_list), (max(dd_obj.window_list) + 1), 1000),
        ylim=[0.0, 0.08],  # MAGIC
    )
    # MAGIC, add panel 'A' label to the left side graph
    ax1.text(
        -0.01,
        1.05,
        "A",
        transform=ax1.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="right",
    )

    for key, val in data_intra_dict.items():
        ax2.scatter(
            0,
            val,
            label=key,
        )

    ax2.set(xlabel="Intragenic TEs", xticks=[])
    ax2.legend(loc="upper right", bbox_to_anchor=(0.76, 0.8))  # MAGIC
    # Heavily tailored location of legend in the intragenic plot, could just use
    # loc='center' but this occasionally obscured a point

    # MAGIC, add panel 'B' label to the center graph
    ax2.text(
        -0.01,
        1.05,
        "B",
        transform=ax2.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="right",
    )

    for key, val in data_right_dict.items():
        ax3.plot(
            dd_obj.window_list,
            val,
            label=key,
            linestyle=(0, (3, 1, 1, 1)),
            marker="o",
        )
    ax3.set(
        xlabel="BP Downstream",
        # yticks=np.arange(0, 0.51, 0.025),  # MAGIC
        xlim=[min(dd_obj.window_list), max(dd_obj.window_list)],
        xticks=range(min(dd_obj.window_list), (max(dd_obj.window_list) + 1), 1000),
        ylim=[0.0, 0.08],  # MAGIC
    )
    ax3.yaxis.tick_right()

    # MAGIC, add panel 'C' label to the right side graph
    ax3.text(
        -0.01,
        1.05,
        "C",
        transform=ax3.transAxes,
        fontsize=16,
        fontweight="bold",
        va="top",
        ha="right",
    )

    filename_to_save = os.path.join(
        output_dir,
        str(order_or_super + "_" + dd_obj.genome_id + "_Combined_Density_Plot.png"),
    )

    logger.info("Saving graphic to: %s" % filename_to_save)

    fig.suptitle(
        "Average TE Density of All Genes as a Function of Window Size and Location"
    )
    plt.savefig(
        filename_to_save,
        bbox_inches="tight",
    )
    if display:
        plt.show()


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.abspath(os.path.join(dir_main, "../", "results/graphs"))
    parser = argparse.ArgumentParser(description="generate graphs")

    parser.add_argument(
        "one_chromosome_hdf5_file",
        type=str,
        help="parent path to an HDF5 file of TE density data",
    )

    parser.add_argument(
        "arabidopsis_gene_data",
        type=str,
        help="parent path to Arabidopsis' filtered gene data file",
    )
    parser.add_argument(
        "chromosome_id",
        type=int,
        help="ID of the chromosome of the HDF5 data, so that we may appropriately subset the GeneData file, which is all chromosomes",
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
    args.one_chromosome_hdf5_file = os.path.abspath(args.one_chromosome_hdf5_file)
    args.arabidopsis_gene_data = os.path.abspath(args.arabidopsis_gene_data)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Begin reading files:
    # Get the genes:
    gene_data = import_filtered_genes(args.arabidopsis_gene_data, logger)
    # Subset the pandas file to get the appropriate chromosome only
    gene_data = gene_data.loc[gene_data["Chromosome"] == args.chromosome_id]
    logger.info("Initializing GeneData from cleaned annotation file")
    # Wrap as GeneData
    gene_data = GeneData(gene_data, str("Arabidopsis_" + str(args.chromosome_id)))

    logger.info("Initializing DensityData from HDF5 file")
    processed_arabidopsis_density_data = DensityData.verify_h5_cache(
        args.one_chromosome_hdf5_file, gene_data, logger
    )

    plot_density_all(
        processed_arabidopsis_density_data,
        "Order",
        args.output_dir,
        logger,
        display=False,
    )

    plot_density_all(
        processed_arabidopsis_density_data,
        "Superfamily",
        args.output_dir,
        logger,
        display=False,
    )
