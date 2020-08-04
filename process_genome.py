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
from transposon.replace_names import te_annot_renamer
from transposon.verify_cache import (verify_chromosome_h5_cache,
    verify_TE_cache, verify_gene_cache, revise_annotation)
from transposon.revise_annotation import Revise_Anno


def validate_args(args, logger):
    """Raise if an input argument is invalid."""

    if not os.path.isfile(args.genes_input_file):
        logger.critical("argument 'genes_input_dir' is not a file")
        raise ValueError("%s is not a file" % (args.genes_input_file))
    if not os.path.isfile(args.tes_input_file):
        logger.critical("argument 'tes_input_dir' is not a file")
        raise ValueError("%s is not a file" % (args.tes_input_file))
    if not os.path.isdir(args.overlap_dir):
        logger.critical("argument 'overlap_dir' is not a directory")
        raise ValueError("%s is not a directory" % (args.overlap_dir))
    if not os.path.isdir(args.output_dir):
        logger.critical("argument 'output_dir' is not a directory")
        raise ValueError("%s is not a directory" % (args.output_dir))


def process(alg_parameters, gene_data_unwrapped, te_data_unwrapped, overlap_dir,
            genome_id, filtered_input_data, reset_h5, h5_cache_location,
            genes_input_file, tes_input_file):
    """
    Run the algorithm

    Args:
        alg_parameters (dict): A dictionary containing the parameters used to
    tune the windowing function. This dictionary is created from the config
    file in main.

        overlap_dir (string): A string path of the directory path to output the
        overlap files. This comes from the ArgumentParser obj and defaults to
        /tmp. You can edit the location of the directory with the -s flag when
        calling density.py. overlap_dir is used when calling OverlapWorker.

        reset_h5 (bool): True/False whether or not to overwrite the H5 cache of
        GeneData and TransposonData if it currently exists
    """
    grouped_genes = split(gene_data_unwrapped, 'Chromosome')
    grouped_TEs = split(te_data_unwrapped, 'Chromosome')
    # check docstring for my split func
    check_groupings(grouped_genes, grouped_TEs, logger, genome_id)

    gene_progress = tqdm(
        total=len(grouped_genes), desc="chromosome  ", position=0, ncols=80)
    _temp_count = 0

    for chromosome_of_gene_data, chromosome_of_te_data in zip(grouped_genes,
                                                              grouped_TEs):
        wrapped_genedata_by_chrom = GeneData(chromosome_of_gene_data.copy(deep=True))
        wrapped_genedata_by_chrom.add_genome_id(genome_id)
        wrapped_tedata_by_chrom = TransposonData(chromosome_of_te_data.copy(deep=True))
        wrapped_tedata_by_chrom.add_genome_id(genome_id)
        chrom_id = wrapped_genedata_by_chrom.chromosome_unique_id

        # NOTE
        # MICHAEL, these are how the H5 files are saved
        h5_g_filename = os.path.join(h5_cache_location, str(chrom_id +
                                     '_GeneData.h5'))
        h5_t_filename = os.path.join(h5_cache_location, str(chrom_id +
                                     '_TEData.h5'))

        verify_chromosome_h5_cache(wrapped_genedata_by_chrom, wrapped_tedata_by_chrom,
                                   h5_g_filename, h5_t_filename, reset_h5, h5_cache_location,
                                   genes_input_file, tes_input_file, chrom_id, logger)

        def window_it(temp_param):
            return range(temp_param[first_window_size],
                         temp_param[last_window_size],
                         temp_param[window_delta])

        n_genes = sum(1 for g in wrapped_genedata_by_chrom.names)
        # TODO create status bar above, reuse and reset here
        # 'total' is a public member
        sub_progress = tqdm(total=n_genes, desc="  genes     ", position=1, ncols=80)
        overlap = OverlapWorker(overlap_dir)

        def progress():
            sub_progress.update(1)
            gene_progress.refresh()
        overlap.calculate(
            wrapped_genedata_by_chrom, wrapped_tedata_by_chrom,
            window_it(alg_parameters), wrapped_genedata_by_chrom.names, progress)

        _temp_count += 1
        gene_progress.update(1)


if __name__ == '__main__':
    """Command line interface to calculate density."""

    parser = argparse.ArgumentParser(description="calculate TE density")
    path_main = os.path.abspath(__file__)
    parser.add_argument('genes_input_file', type=str,
                        help='parent path of gene file')
    parser.add_argument('tes_input_file', type=str,
                        help='parent path of transposon file')
    parser.add_argument('genome_id', type=str,
                        help='string of the genome to be run, for clarity')
    parser.add_argument('--config_file', '-c', type=str,
                        default=os.path.join(path_main, '../',
                                             'config/test_run_config.ini'),
                        help='parent path of config file')
    parser.add_argument('--overlap_dir', '-s', type=str,
                        default=os.path.abspath('/tmp'),
                        help='parent directory to output overlap data')
    parser.add_argument('--filtered_input_data', '-f', type=str,
                        default=os.path.join(path_main, '../',
                                             'filtered_input_data'),
                        help='parent directory for cached input data')
    # TODO if the user sets their own filtered_input_data location, this
    # h5_cache_loc location will not exactly follow their filtered_input_data
    # convention
    parser.add_argument('--input_h5_cache_loc', '-h5', type=str,
                        default=os.path.join(path_main, '../',
                                             'filtered_input_data/input_h5_cache'),
                        help='parent directory for h5 cached input data')
    parser.add_argument('--reset_h5', action='store_true')
    parser.add_argument('--contig_del', action='store_false')
    parser.add_argument('--revise_anno', action='store_true')

    parser.add_argument('--revised_input_data', '-r', type=str,
                        default=os.path.join(path_main, '../',
                                             'filtered_input_data/revised_input_data'),
                        help='parent directory for cached revised input data')
    parser.add_argument('--output_dir', '-o', type=str,
                        default=os.path.join(path_main, '../', 'results'),
                        help='parent directory to output results')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='set debugging level to DEBUG')

    args = parser.parse_args()
    args.genes_input_file = os.path.abspath(args.genes_input_file)
    args.tes_input_file = os.path.abspath(args.tes_input_file)
    args.config_file = os.path.abspath(args.config_file)
    args.overlap_dir = os.path.abspath(args.overlap_dir)
    args.filtered_input_data = os.path.abspath(args.filtered_input_data)
    args.input_h5_cache_loc = os.path.abspath(args.input_h5_cache_loc)
    args.revised_input_data = os.path.abspath(args.revised_input_data)
    args.output_dir = os.path.abspath(args.output_dir)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Want a higher recursion limit for the code
    sys.setrecursionlimit(11**6)

    # logger.info("Start processing directory '%s'"%(args.input_dir))
    for argname, argval in vars(args).items():
        logger.debug("%-12s: %s" % (argname, argval))
    validate_args(args, logger)

    # NOTE Imports
    # FUTURE move this preprocessing to it's object

    logger.info("Checking disk for previously filtered data...")
    # NOTE
    # MAGIC NUMBER for g_fname, and t_fname trying to get filename without
    # extension, will produce an unexpected, but functional filename if their
    # input filename has multiple . in it
    g_fname = os.path.basename(os.path.splitext(args.genes_input_file)[0])
    t_fname = os.path.basename(os.path.splitext(args.tes_input_file)[0])
    # the genome input file, maybe make it user provided
    cleaned_genes = os.path.join(args.filtered_input_data, str('Cleaned_' +
                                                               g_fname +
                                                               '.tsv'))
    cleaned_transposons = os.path.join(args.filtered_input_data, str('Cleaned_' +
                                                                     t_fname +
                                                                     '.tsv'))
    revised_transposons = os.path.join(args.revised_input_data, str('Revised_'
                                                                    + t_fname
                                                                     + '.tsv'))
    gene_data_unwrapped = verify_gene_cache(args.genes_input_file,
                                            cleaned_genes, args.contig_del,
                                            logger)
    te_data_unwrapped = verify_TE_cache(args.tes_input_file, cleaned_transposons,
                              te_annot_renamer, args.contig_del,
                              logger)

    # Revise the TE annotation
    te_data_unwrapped = revise_annotation(te_data_unwrapped, args.revise_anno,
                                          revised_transposons,
                                          args.revised_input_data,
                                          logger, args.genome_id)

    # NOTE neither Gene_Data or TE_Data are wrapped yet

    logger.info("Reading config file and making parameter dictionary...")
    parser = ConfigParser()
    parser.read(args.config_file)
    first_window_size = parser.getint('density_parameters', 'first_window_size')
    window_delta = parser.getint('density_parameters', 'window_delta')
    last_window_size = parser.getint('density_parameters', 'last_window_size')
    alg_parameters = {first_window_size: first_window_size,
                      window_delta: window_delta,
                      last_window_size: last_window_size}

    # Process data
    logger.info("Process data...")
    process(alg_parameters, gene_data_unwrapped, te_data_unwrapped,
            args.overlap_dir, args.genome_id, args.filtered_input_data,
            args.reset_h5, args.input_h5_cache_loc, args.genes_input_file,
            args.tes_input_file)
