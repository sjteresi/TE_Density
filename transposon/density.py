#!/usr/bin/env python3

"""
Calculate transposable element density.
"""

__author__ = "Scott Teresi, Michael Teresi"

import os
import time
import argparse
import coloredlogs
import logging

from tqdm import tqdm
import numpy as np

from transposon.data import GeneData, TransposonData
from transposon.import_genes import import_genes
from transposon.import_transposons import import_transposons
from transposon.overlap import OverlapData
from transposon.replace_names import te_annot_renamer


def get_nulls(my_df):
    null_columns = my_df.columns[my_df.isnull().any()]
    count_of_null = my_df[null_columns].isnull().sum()
    print('Counts of null values per column: ' '\n', count_of_null, '\n')
    rows_where_null = my_df[my_df.isnull().any(axis=1)][null_columns].head()
    print('Rows where null exist: ', '\n', rows_where_null, '\n')


def drop_nulls(my_df, status=False):
    if status:
        print('DROPPING ROWS WITH AT LEAST ONE NULL VALUE!!!')
    my_df = my_df.dropna(axis=0, how='any')
    return my_df


def get_unique(my_df_col):
    return my_df_col.unique()


def swap_columns(df, col_condition, c1, c2):
    df.loc[col_condition, [c1, c2]] = \
        df.loc[col_condition, [c2, c1]].values
    return df


def split(df, group):
    """Return list of the dataframe with each element being a subset of the df.

    I use this function to split by chromosome so that we may later do
    chromosome element-wise operations.
    """

    gb = df.groupby(group)
    return [gb.get_group(x) for x in gb.groups]


def check_density_shape(densities, transposon_data):
    """Checks to make sure the density output is of the same dimension as the
    transposon_data input.

    Args:
        densities (numpy.ndarray):
            density calculations, returned from rho functions.
        transposon_data (transposon.data.TransposonData):
            transposon container
    """
    if densities.shape != transposon_data.starts.shape:
        msg = ("Density dataframe shape not the same size as the TE dataframe")
        logger.critical(msg)
        raise ValueError(msg)


def validate_window(window_start, g_start, window_length):
    if window_start < 0:
        msg = ("window_start is not 0 or a positive value")
        logger.critical(msg)
        raise ValueError(msg)
    if window_start == 0:
        window_length = g_start - window_start + 1
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
    win_start = gene_datum.left_win_start(win_length)
    win_stop = gene_datum.left_win_stop()
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
    win_start = gene_datum.right_win_start()
    win_stop = gene_datum.right_win_stop(win_length)

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


def density_algorithm(genes, tes, window, increment, max_window):
    """
    te data frame has columns: SEE import_transposons
    te data frame has rows: each row is a temp

    Args:


    """
    # NOTE create 2 structs to hold files / param & result

    try:
        get_unique(genes.Chromosome) == get_unique(tes.Chromosome)
    except:
        raise ValueError("You do not have the same chromosomes in your files")

    windows = list(range(window, max_window, increment))
    logging.info(" windows are {}:{}:{}  -->  {}"
                 .format(window, increment, max_window, windows))

    # Use the subsets in main?
    while window <= max_window:
        logging.debug(" Gene df shape:  {}".format(genes.values.shape))
        logging.debug(" TE df shape:  {}".format(tes.values.shape))
        # Perform the windowing operations
        # Multiple tests need to be done here for each window
        # All the tests must be run for that window and then rerun for the
        # next window
        # I foresee that for each gene file we will have an output wiht all
        # of the original gene data and the output will be 500_LTR_Upstream
        # The TE types (LTR) are given by the TE_Dataframe.Order and
        # TE_Dataframe.SuperFamily columns

        # The init_empty_densities function was to add the appropriate
        # columns, we may not need to worry about that for now
        # All the commented code below are my attempts to do the work
        # -----------------------------
        window += increment


def init_empty_densities(my_genes, my_tes, window):
    """Initializes all of the empty columns we need in the gene file. """

    # NOTE This function is a candidate for deletion
    Order_List = my_tes.Order.unique()
    SuperFamily_List = my_tes.SuperFamily.unique()
    Directions = ['_downstream', '_intra', '_upstream']
    # left, center, right

    for an_order in Order_List:
        for direction in Directions:
            col_name = (str(window) + '_' + an_order + direction)
            my_genes[col_name] = np.nan

    for a_superfamily in SuperFamily_List:
        for direction in Directions:
            col_name = (str(window) + '_' + a_superfamily + direction)
            my_genes[col_name] = np.nan
    my_genes['TEs_inside'] = np.nan
    return my_genes


def check_groupings(grouped_genes, grouped_TEs, logger):
    """Validates the gene / TE pairs.

    This is just to make sure that each pair of chromosomes are right.
    Correct subsetting would be managed by the custom split command.

    Args:
        grouped_genes (list of pandaframes): Gene dataframes separated by chromosome
        grouped_TEs (list of pandaframes): TE dataframes separated by chromosome
    """
    for g_element, t_element in zip(grouped_genes, grouped_TEs):
        # print(g_element.Chromosome.iloc[:].values[0])
        if g_element.Chromosome.iloc[:].values[0] != t_element.Chromosome.iloc[:].values[0]:
            msg = 'Chromosomes do not match for the grouped_genes or grouped_TEs'
            logger.critical(msg)
            raise ValueError(msg)


def validate_args(args, logger):
    """Raise if an input argument is invalid."""

    if not os.path.isfile(args.genes_input_file):
        logger.critical("argument 'genes_input_dir' is not a file")
        raise ValueError("%s is not a directory" % (abs_path))
    if not os.path.isfile(args.tes_input_file):
        logger.critical("argument 'tes_input_dir' is not a file")
        raise ValueError("%s is not a directory" % (abs_path))
    if not os.path.isdir(args.output_dir):
        logger.critical("argument 'output_dir' is not a directory")
        raise ValueError("%s is not a directory" % (abs_path))


def process():
    grouped_genes = split(Gene_Data, 'Chromosome')  # check docstring for my split func
    grouped_TEs = split(TE_Data, 'Chromosome')  # check docstring for my split func
    check_groupings(grouped_genes, grouped_TEs, logger)
    # Think of the 7 main "chromosomes" as "meta-chromosomes" in reality there
    # are 4 actual chromosomes per "meta-chromosome" label. So Fvb1 is
    # meta-chromosome 1, and within that Fvb1-1 of genes should only be
    # matching with Fvb1-1 of TEs, not Fvb1-2. The first number, what I am
    # calling the "meta-chromosome" is just denoting that it is the first
    # chromosome, where the second number is the actual physical chromosome,
    # and we use the number to denote which subgenome it is assigned to.

    gene_progress = tqdm(total=len(grouped_genes), desc="metachrome  ", position=0, ncols=80)
    for sub_gene, sub_te in zip(grouped_genes, grouped_TEs):
        gene_data = GeneData(sub_gene)
        te_data = TransposonData(sub_te)
        # TODO validate the gene / te pair

        # NOTE Everything seems good to me, the method for validating (my split
        # function and the check_groupings) could probably be done better but I
        # don't know, what do you think Michael? It is nice to be able
        # to put it through a for loop for our workers's purposes, see
        # the above for loop. -S

        # FUTURE multiprocess starting here
        # create workers
        # create accumulators

        window_it = lambda: range(100, 1000, 100)  # TODO remove magic numbers, parametrize
        n_genes = sum(1 for g in gene_data.names)
        sub_progress = tqdm(total=n_genes, desc="  genes     ", position=1, ncols=80)
        overlap = OverlapData(gene_data, te_data)
        def update_win_prog():
            sub_progress.update(1)
            gene_progress.refresh()
        overlap.calculate(window_it(), gene_data.names, update_win_prog)
        # for window in window_it():
        #     # create density request, push
        #
        #
        #     # density_algorithm(
        #     #                 Gene_Data,
        #     #                 TE_Data,
        #     #                 window=1000,
        #     #                 increment=500,
        #     #                 max_window=10000
        #     #                 )
        #
        #     time.sleep(0.075)
        #     window_progress.update(1)
        #     gene_progress.refresh()
        # collapse accumulated results (i.e. do the division)
        # combine all the results
        # write to disk
        gene_progress.update(1)


if __name__ == '__main__':
    """Command line interface to calculate density."""

    parser = argparse.ArgumentParser(description="calculate TE density")
    path_main = os.path.abspath(__file__)
    parser.add_argument('genes_input_file', type=str,
                        help='parent path of gene file')
    parser.add_argument('tes_input_file', type=str,
                        help='parent path of transposon file')
    parser.add_argument('--output_dir', '-o', type=str,
                        default=os.path.join(path_main, '../..', 'results'),
                        help='parent directory to output results')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='set debugging level to DEBUG')
    args = parser.parse_args()
    args.output_dir = os.path.abspath(args.output_dir)
    args.genes_input_file = os.path.abspath(args.genes_input_file)
    args.tes_input_file = os.path.abspath(args.tes_input_file)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # logger.info("Start processing directory '%s'"%(args.input_dir))
    for argname, argval in vars(args).items():
        logger.debug("%-12s: %s" % (argname, argval))
    validate_args(args, logger)

    # FUTURE move this preprocessing to it's object
    logger.info("Importing genes, this may take a moment...")
    Gene_Data = import_genes(args.genes_input_file)
    logger.info("Importing transposons, this may take a moment...")
    TE_Data = import_transposons(args.tes_input_file, te_annot_renamer)
    
    print(TE_Data.head())

    process()


    # Scott Test
    # genes = GeneData(Gene_Data)
    # print(genes.unique_genes)
    #print(genes.get_gene('maker-Fvb1-1-snap-gene-0.15'))
    #TEs = TransposonData(TE_Data)
    #print(TEs)
    # print(genes.get_gene('maker-Fvb1-1-snap-gene-0.15').start_stop_len)
    # genes.get_gene('maker-Fvb1-1-snap-gene-0.15', 500)
    # print(genes.get_gene('maker-Fvb1-1-snap-gene-0.15').left_win_start)
    # print(genes.get_gene('maker-Fvb1-1-snap-gene-0.15').calc_win_length(500))
    # print(genes.get_gene('maker-Fvb1-1-snap-gene-0.15').start)
    # print(genes.get_gene('maker-Fvb1-1-snap-gene-0.15').stop)
    # print(genes.get_gene('maker-Fvb1-1-snap-gene-0.15').left_win_start(50))
    # print(genes.get_gene('maker-Fvb1-1-snap-gene-0.15').left_win_stop())
