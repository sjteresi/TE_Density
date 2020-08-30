#!/usr/bin/env python3

"""
Calculate transposable element density.
"""

__author__ = "Scott Teresi, Michael Teresi"

import argparse
import os

import logging
import coloredlogs
import numpy as np
import pandas as pd
from tqdm import tqdm
from configparser import ConfigParser
import sys
import time

from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.overlap import OverlapWorker


def validate_args(args, logger):
    """Raise if an input argument is invalid."""

    if not os.path.isfile(args.genes_input_file):
        logger.critical("argument 'genes_input_dir' is not a file")
        raise ValueError("%s is not a directory" % (args.genes_input_file))
    if not os.path.isfile(args.tes_input_file):
        logger.critical("argument 'tes_input_dir' is not a file")
        raise ValueError("%s is not a directory" % (args.tes_input_file))
    if not os.path.isdir(args.overlap_dir):
        logger.critical("argument 'overlap_dir' is not a directory")
        raise ValueError("%s is not a directory" % (args.overlap_dir))
    if not os.path.isdir(args.output_dir):
        logger.critical("argument 'output_dir' is not a directory")
        raise ValueError("%s is not a directory" % (args.output_dir))


if __name__ == "__main__":
    """Command line interface to calculate density."""

    parser = argparse.ArgumentParser(description="calculate TE density")
    path_main = os.path.abspath(__file__)
    parser.add_argument("genes_input_file", type=str, help="parent path of gene file")
    parser.add_argument(
        "tes_input_file", type=str, help="parent path of transposon file"
    )
    parser.add_argument(
        "genome_id", type=str, help="string of the genome to be run, for clarity"
    )
    parser.add_argument(
        "--config_file",
        "-c",
        type=str,
        default=os.path.join(path_main, "../../", "config/test_run_config.ini"),
        help="parent path of config file",
    )
    parser.add_argument(
        "--overlap_dir",
        "-s",
        type=str,
        default=os.path.abspath("/tmp"),
        help="parent directory to output overlap data",
    )
    parser.add_argument(
        "--filtered_input_data",
        "-f",
        type=str,
        default=os.path.join(path_main, "../..", "filtered_input_data"),
        help="parent directory for cached input data",
    )
    # NOTE
    # TODO if the user sets their own filtered_input_data location, this
    # h5_cache_loc location will not exactly follow their filtered_input_data
    # convention
    # ANSWER then don't allow them to do so
    # also what does 'not exactly follow' mean and what is the result?
    parser.add_argument(
        "--input_h5_cache_loc",
        "-h5",
        type=str,
        default=os.path.join(
            path_main, "../..", "filtered_input_data/chromosome_h5_cache"
        ),
        help="parent directory for h5 cached chromosome input data",
    )
    parser.add_argument("--reset_h5", action="store_true")
    parser.add_argument("--contig_del", action="store_false")
    parser.add_argument("--revise_anno", action="store_true")

    parser.add_argument(
        "--revised_input_data",
        "-r",
        type=str,
        default=os.path.join(
            path_main, "../..", "filtered_input_data/revised_input_data"
        ),
        help="parent directory for cached revised input data",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=os.path.join(path_main, "../..", "results"),
        help="parent directory to output results",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.genes_input_file = os.path.abspath(args.genes_input_file)
    args.tes_input_file = os.path.abspath(args.tes_input_file)
    args.config_file = os.path.abspath(args.config_file)
    args.overlap_dir = os.path.abspath(args.overlap_dir)
    args.filtered_input_data = os.path.abspath(args.filtered_input_data)
    args.input_h5_cache_loc = os.path.abspath(args.input_h5_cache_loc)
    args.revised_input_data = os.path.abspath(
        os.path.join(args.revised_input_data, args.genome_id)
    )

    # Create directories
    if not os.path.exists(args.filtered_input_data):
        os.makedirs(args.filtered_input_data)
    if not os.path.exists(args.input_h5_cache_loc):
        os.makedirs(args.input_h5_cache_loc)
    if not os.path.exists(args.revised_input_data):
        os.makedirs(args.revised_input_data)

    args.output_dir = os.path.abspath(args.output_dir)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Want a higher recursion limit for the code
    sys.setrecursionlimit(11 ** 6)

    # logger.info("Start processing directory '%s'"%(args.input_dir))
    for argname, argval in vars(args).items():
        logger.debug("%-12s: %s" % (argname, argval))
    validate_args(args, logger)


    # TODO move all preprocessing out of this input file
    raise NotImplementedError()

    logger.info("Reading config file and making parameter dictionary...")
    parser = ConfigParser()
    parser.read(args.config_file)
    first_window_size = parser.getint("density_parameters", "first_window_size")
    window_delta = parser.getint("density_parameters", "window_delta")
    last_window_size = parser.getint("density_parameters", "last_window_size")
    alg_parameters = {
        first_window_size: first_window_size,
        window_delta: window_delta,
        last_window_size: last_window_size,
    }

    # Process data
    logger.info("Process data...")
    process(
        alg_parameters,
        gene_data_unwrapped,
        te_data_unwrapped,
        args.overlap_dir,
        args.genome_id,
        args.filtered_input_data,
        args.reset_h5,
        args.input_h5_cache_loc,
        args.genes_input_file,
        args.tes_input_file,
    )
