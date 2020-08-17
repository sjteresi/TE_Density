#!/usr/bin/env python3

"""
Preprocessing for input data.
"""

__author__ = "Scott Teresi, Michael Teresi"

import logging
import coloredlogs
import numpy as np
import pandas as pd
from tqdm import tqdm

from transposon.gene_data import GeneData


def split(dataframe, group):
    """Return list of dataframes with each element being a subset of the df.

    I use this function to split by chromosome so that we may later do
    chromosome element-wise operations.

    This function is also used in revise_annotation.py to split on transposon
    identities.
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
    # NOTE: Candidate for deletion but this is used in the rho functions (also
    # candidates for deletion).
    if densities.shape != transposon_data.starts.shape:
        msg = "Density dataframe shape not the same size as the TE dataframe"
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
        msg = "window_start is not 0 or a positive value"
        logger.critical(msg)
        raise ValueError(msg)
    if left_window_start == 0:
        window_length = left_window_stop - left_window_start + 1
    return window_length


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
        if (
            g_element.Chromosome.iloc[:].values[0]
            != t_element.Chromosome.iloc[:].values[0]
        ):
            msg = "Chromosomes do not match for the grouped_genes or grouped_TEs"
            logger.critical(msg)
            raise ValueError(msg)
        try:
            sub_gene = GeneData(g_element, genome_id)
            subgene_uid = sub_gene.chromosome_unique_id
        except RuntimeError as r_err:
            logging.critical("sub gene grouping is not unique: {}".format(sub_gene))
            raise r_err
