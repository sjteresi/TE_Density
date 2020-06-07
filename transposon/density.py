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
import time

from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.overlap import OverlapWorker
from transposon.replace_names import te_annot_renamer
from transposon.verify_cache import (verify_chromosome_h5_cache,
    verify_TE_cache, verify_gene_cache)


def get_nulls(my_df):
    """
    Print out a count of null values per column.
    Print out the row IDs where the null values exist

    Args:
        my_df (Pandaframes): Pandaframe to check null values in
    """
    null_columns = my_df.columns[my_df.isnull().any()]
    count_of_null = my_df[null_columns].isnull().sum()
    print('Counts of null values per column: ' '\n', count_of_null, '\n')
    rows_where_null = my_df[my_df.isnull().any(axis=1)][null_columns].head()
    print('Rows where null exist: ', '\n', rows_where_null, '\n')


def drop_nulls(my_df, status=False):
    """
    Drop null values inside a Pandaframe

    Args:
        my_df (Pandaframes): Pandaframe to drop null values
    """
    if status:
        print('DROPPING ROWS WITH AT LEAST ONE NULL VALUE!!!')
    my_df = my_df.dropna(axis=0, how='any')
    return my_df


def swap_columns(dataframe, col_condition, col_1, col_2):
    """
    Swap the values of two designated columns for a row based on a column
    condition, return the dataframe

    Args:
        my_df (Pandaframes): Pandaframe to swap columns in.
    """
    dataframe.loc[col_condition, [col_1, col_2]] = \
        dataframe.loc[col_condition, [col_2, col_1]].values
    return dataframe


def split(dataframe, group):
    """Return list of the dataframe with each element being a subset of the df.

    I use this function to split by chromosome so that we may later do
    chromosome element-wise operations.
    """

    grouped_df = dataframe.groupby(group)
    return [grouped_df.get_group(x) for x in grouped_df.groups]


def check_density_shape(densities, transposon_data):
    """Checks to make sure the density output is of the same dimension as the
    transposon_data input.

    Args:
        densities (numpy.ndarray):
            density calculations, returned from rho functions.
        transposon_data (TransposonData):
            transposon container
    """
    if densities.shape != transposon_data.starts.shape:
        msg = ("Density dataframe shape not the same size as the TE dataframe")
        logger.critical(msg)
        raise ValueError(msg)


def validate_window(left_window_start, left_window_stop, window_length):
    """
    Validate the window specifically to make sure it doesn't go into the
    negative values. Function invoked to validate the left-hand window.

    Args:
        window_start (int): integer value for the left window start value, we need
        to make sure it isn't negative.

        TODO add description

        window_length (int)
    """
    if left_window_start < 0:
        msg = ("window_start is not 0 or a positive value")
        logger.critical(msg)
        raise ValueError(msg)
    if left_window_start == 0:
        window_length = left_window_stop - left_window_start + 1
    return window_length


def rho_left_window(gene_data, gene_name, transposon_data, window):
    """Density to the left (downstream) of a gene.
    When TE is between gene and window

    Relevant things: No. TE bases / window
    GeneStart so that we may calculate the WindowStart
    WindowStart is unique to each gene. Calculated via Gstart - Window
    WindowStart is the leftmost window

    Args:
        window (int): integer value of the current window
        gene_data (transponson.data.GeneData): gene container
        gene_name (hashable): name of gene to use
        transposon_data (transponson.data.TransposonData): transposon container
    """
    transposon_data.check_shape()
    gene_datum = gene_data.get_gene(gene_name)
    g_start, g_stop, g_length = gene_datum.start_stop_len

    # Define windows
    win_length = gene_datum.win_length(window)
    win_start = gene_datum.left_win_start(window)
    win_stop = gene_datum.left_win_stop
    win_length = validate_window(win_start, g_start, win_length)

    # Set bounds and perform density calculation
    lower_bound = np.maximum(win_start, transposon_data.starts)
    upper_bound = np.minimum(win_stop, transposon_data.stops)
    te_overlaps = np.maximum(0, (upper_bound - lower_bound + 1))
    densities = np.divide(
        te_overlaps,
        win_length,
        out=np.zeros_like(te_overlaps, dtype='float')
    )
    check_density_shape(densities, transposon_data)
    return densities


def rho_intra(gene_data, gene_name, transposon_data):
    """Intra density for one gene wrt transposable elements.

    The relevant ares is the gene for intra density.

    Args:
        gene_data (transponson.data.GeneData): gene container
        gene_name (hashable): name of gene to use
        transposon_data (transponson.data.TransposonData): transposon container
    """
    transposon_data.check_shape()
    gene_datum = gene_data.get_gene(gene_name)
    g_start, g_stop, g_length = gene_datum.start_stop_len
    lower = np.minimum(g_stop, transposon_data.stops)
    upper = np.maximum(g_start, transposon_data.starts)
    te_overlaps = np.maximum(0, (lower - upper + 1))

    densities = np.divide(
        te_overlaps,
        g_length,
        out=np.zeros_like(te_overlaps, dtype='float'),
        where=g_length != 0
    )

    check_density_shape(densities, transposon_data)
    return densities


def rho_right_window(gene_data, gene_name, transposon_data, window):
    """Density to the right (upstream) of a gene.
    When TE is between gene and window

    Relevant things: No. TE bases / window
    GeneStop so that we may calculate the WindowStop
    WindowStop is unique to each gene. Calculated via Gstop + Window
    WindowStop is the rightmost window

    Args:
        window (int): integer value of the current window
        gene_data (transponson.data.GeneData): gene container
        gene_name (hashable): name of gene to use
        transposon_data (transponson.data.TransposonData): transposon container
    """
    transposon_data.check_shape()
    gene_datum = gene_data.get_gene(gene_name)
    g_start, g_stop, g_length = gene_datum.start_stop_len

    # Define windows
    win_length = gene_datum.win_length(window)
    win_start = gene_datum.right_win_start
    win_stop = gene_datum.right_win_stop(window)

    # Set bounds and perform density calculation
    lower_bound = np.maximum(win_start, transposon_data.starts)
    upper_bound = np.minimum(win_stop, transposon_data.stops)
    te_overlaps = np.maximum(0, (upper_bound - lower_bound + 1))
    densities = np.divide(
        te_overlaps,
        win_length,
        out=np.zeros_like(te_overlaps, dtype='float')
    )
    check_density_shape(densities, transposon_data)
    return densities


def init_empty_densities(my_genes, my_tes, window):
    """Initializes all of the empty columns we need in the gene file. """

    # NOTE This function is a candidate for deletion
    order_list = my_tes.Order.unique()
    superfamily_list = my_tes.SuperFamily.unique()
    directions = ['_downstream', '_intra', '_upstream']
    # left, center, right

    for an_order in order_list:
        for direction in directions:
            col_name = (str(window) + '_' + an_order + direction)
            my_genes[col_name] = np.nan

    for a_superfamily in superfamily_list:
        for direction in directions:
            col_name = (str(window) + '_' + a_superfamily + direction)
            my_genes[col_name] = np.nan
    my_genes['TEs_inside'] = np.nan
    return my_genes


def check_groupings(grouped_genes, grouped_TEs, logger, genome_id):
    """Validates the gene / TE pairs.

    This is just to make sure that each pair of chromosomes are right.
    Correct subsetting would be managed by the custom split command.

    Args:
        grouped_genes (list of pandaframes): Gene dataframes separated by chromosome
        grouped_TEs (list of pandaframes): TE dataframes separated by chromosome
        genome_id (str) a string of the genome name.
    """
    for g_element, t_element in zip(grouped_genes, grouped_TEs):
        # print(g_element.Chromosome.iloc[:].values[0])
        if g_element.Chromosome.iloc[:].values[0] != t_element.Chromosome.iloc[:].values[0]:
            msg = 'Chromosomes do not match for the grouped_genes or grouped_TEs'
            logger.critical(msg)
            raise ValueError(msg)
        try:
            sub_gene = GeneData(g_element, genome_id)
            subgene_uid = sub_gene.chromosome_unique_id
        except RuntimeError as r_err:
            logging.critical("sub gene grouping is not unique: {}".format(sub_gene))
            raise r_err


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


def process(alg_parameters, Gene_Data, TE_Data, overlap_dir, genome_id,
            filtered_input_data, reset_h5, h5_cache_location, genes_input_file,
            tes_input_file):
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




    grouped_genes = split(Gene_Data, 'Chromosome')  # check docstring for my split func
    grouped_TEs = split(TE_Data, 'Chromosome')  # check docstring for my split func
    check_groupings(grouped_genes, grouped_TEs, logger, genome_id)
    # Think of the 7 main "chromosomes" as "meta-chromosomes" in reality there
    # are 4 actual chromosomes per "meta-chromosome" label. So Fvb1 is
    # meta-chromosome 1, and within that Fvb1-1 of genes should only be
    # matching with Fvb1-1 of TEs, not Fvb1-2. The first number, what I am
    # calling the "meta-chromosome" is just denoting that it is the first
    # chromosome, where the second number is the actual physical chromosome,
    # and we use the number to denote which subgenome it is assigned to.

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

def str2bool(v):
    """
    Make sure the argparse boolean arguments work correctly.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('True', 'true'):
        return True
    elif v.lower() in ('False', 'false'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected')

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
    parser.add_argument('--contig_del', nargs='?', type=str2bool, const=True,
                        default=True)
    parser.add_argument('--config_file', '-c', type=str,
                        default=os.path.join(path_main, '../../',
                                             'config/test_run_config.ini'),
                        help='parent path of config file')
    parser.add_argument('--overlap_dir', '-s', type=str,
                        default=os.path.abspath('/tmp'),
                        help='parent directory to output overlap data')
    parser.add_argument('--filtered_input_data', '-f', type=str,
                        default=os.path.join(path_main, '../..',
                                             'filtered_input_data'),
                        help='parent directory for cached input data')
    parser.add_argument('--h5_cache_loc', '-h5', type=str,
                        default=os.path.join(path_main, '../..',
                                             'filtered_input_data/h5_cache'),
                        help='parent directory for h5 cached input data')
    parser.add_argument('--reset_h5', nargs='?', type=str2bool, const=True, default=False)
    parser.add_argument('--output_dir', '-o', type=str,
                        default=os.path.join(path_main, '../..', 'results'),
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
    args.h5_cache_loc = os.path.abspath(args.h5_cache_loc)
    args.output_dir = os.path.abspath(args.output_dir)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)


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
    Gene_Data = verify_gene_cache(args.genes_input_file, cleaned_genes, args.contig_del, logger)
    TE_Data = verify_TE_cache(args.tes_input_file, cleaned_transposons,
                              te_annot_renamer, args.contig_del, logger)
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
    process(alg_parameters, Gene_Data, TE_Data, args.overlap_dir,
            args.genome_id, args.filtered_input_data, args.reset_h5,
            args.h5_cache_loc, args.genes_input_file, args.tes_input_file)
