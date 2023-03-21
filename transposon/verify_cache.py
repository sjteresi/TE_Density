#!/usr/bin/env python3

"""
Contains methods for verifying various cached input data.
"""

__author__ = "Scott Teresi"

import os
import sys
import pandas as pd
from transposon.revise_annotation import ReviseAnno
from transposon.import_filtered_genes import import_filtered_genes
from transposon.import_filtered_TEs import import_filtered_TEs


def verify_chromosome_h5_cache(
    gene_data_obj,
    te_data_obj,
    g_filepath,
    t_filepath,
    reset_h5,  # TODO edit this later
    cache_location,
    genes_input_file,
    tes_input_file,
    chrom_id,
    logger,
):
    """Determine whether or not previously saved gene_data and TransposonData
    exist in .tsv format. Each file represents either gene_data or
    TransposonData for one chromosome at a time. Save/update files as
    necessary.

    When are the tsv cache files written?
    1. Files will be written if there are no current corresponding files
    saved on disk.
    2. If a command-line option is passed to density.py to re-write the
    files, this option defaults to NOT re-write.
    # TODO check
    3. If enough time has passed between the creation of the H5 file and the
    current run-time of the program. TODO talk to Michael more about this.

    Args:
        gene_data_obj (gene_data): Instance of gene_data
        te_data_obj (TransposonData): Instance of TransposonData
        g_filepath (str): The string of the filename in which to save the
            gene_data as an H5 file.
        t_filepath (str): The string of the filename in which to save the
            TransposonData as an H5 file.
        reset_h5 (bool): Boolean, whether or not to completely rewrite the
            cache of all H5 files. True means that we will rewrite.
        cache_location (str): The location (file path) in which to store the
            cache gene and TE data files.
        genes_input_file (str): The file path of the file that was used to
            generate the gene_data instance.
        tes_input_file (str): The file path of the file that was used to
            generate the TransposonData instance.
        chrom_id (str): A string representation of the current chromosome. Used
        to name each H5 file.
    """
    if reset_h5:
        logger.info("overwrite: %s" % g_filepath)
        logger.info("overwrite: %s" % t_filepath)
        gene_data_obj.write(g_filepath)
        te_data_obj.write(t_filepath)

    if os.path.exists(g_filepath) and os.path.exists(t_filepath):
        gene_annot_time = os.path.getmtime(genes_input_file)
        te_annot_time = os.path.getmtime(tes_input_file)
        gene_h5_time = os.path.getmtime(g_filepath)
        te_h5_time = os.path.getmtime(t_filepath)

        if (gene_annot_time > gene_h5_time) and (te_annot_time > te_h5_time):
            logger.info("cache is too old for chromosome '%s'" % chrom_id)
            logger.info("write: %s" % g_filepath)
            logger.info("write: %s" % t_filepath)
            gene_data_obj.write(g_filepath)
            te_data_obj.write(t_filepath)

        elif (gene_annot_time < gene_h5_time) and (te_annot_time < te_h5_time):
            # No need to re-write a current cache
            return

    elif reset_h5 or (not (os.path.exists(g_filepath) and os.path.exists(t_filepath))):
        gene_data_obj.write(g_filepath)
        te_data_obj.write(t_filepath)
    else:
        logger.critical(
            """During the verification of the H5 cache nothing was
                       saved because 0 conditions were met."""
        )


def verify_TE_cache(tes_input_file, logger):
    """Read a preprocessed/filtered TE annotation file from disk; return a
    pandaframe of the file, no modifications are made to the data.

    Args:
        tes_input_file (str): A command line argument, this is the location
            of the processed TE annotation file.

    Returns:
        te_data (pandas.DataFrame): A pandas dataframe of the TE data
    """
    logger.info("Reading pre-filtered TE annotation file: %s" % tes_input_file)
    te_data = import_filtered_TEs(tes_input_file, logger)
    return te_data


def verify_gene_cache(genes_input_file, logger):
    """Read a preprocessed/filtered gene annotation file from disk; return a
    pandaframe of the file, no modifications are made to the data.

    Args:
        genes_input_file (str): A command line argument, this is the location
            of the processed gene annotation file.

    Returns:
        gene_data (pandas.DataFrame): the gene data container
    """
    logger.info("Reading pre-filtered gene annotation file: %s" % genes_input_file)
    gene_data = import_filtered_genes(genes_input_file, logger)
    return gene_data


def revise_annotation(
    te_data, revise_anno, revised_transposons_loc, revised_cache_loc, logger, genome_id
):
    """Remove overlapping elements of the same type.

    Revises the annotation so that elements of the same type do not overlap at
    all. Will essentially merge elements together, elongating them. This is
    done so that the mathematics of density make sense. You can elect to not
    use the revised annotation through a command-line argument to density.py,
    however given that TEs often overlap with one another in annotatios (not
    just being nested in one another) it can lead to some difficulties in
    accurately assessing density and obfuscate the density results.

    Args:
        te_data (pandas.core.DataFrame): A PandaFrame of the TE data,
        previously imported from raw and filtered or imported from a previously
        filtered data file that was saved to disk.

        revise_anno (bool): A boolean of whether or not to use/create a revised
        annotation

        revised_transposons (str): A string representing the path of a
        previously filtered (cleaned) and revised TE annotation.

        revised_cache_loc (): Directory for output files

        logger ():

        genome_id (str): String of the genome ID

    Returns:
        te_data (pandaframe): A pandas dataframe of the revised TE data
    """

    if os.path.exists(revised_transposons_loc) and not revise_anno:
        logger.info("load revised TE: %s" % revised_transposons_loc)
        te_data = import_filtered_TEs(revised_transposons_loc, logger)
    else:
        logger.info("creating revised TE dataset...")
        logger.info("revising the TE dataset will take a long time!")
        # N.B we want higher recursion limit for the code
        sys.setrecursionlimit(11 ** 6)
        revised_te_data = ReviseAnno(
            te_data, revised_transposons_loc, revised_cache_loc, genome_id
        )
        revised_te_data.create_superfam()
        revised_te_data.create_order()
        revised_te_data.create_nameless()
        revised_te_data.verify_files()
        te_data = revised_te_data.whole_te_annotation
    return te_data
