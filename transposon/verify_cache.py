#!/usr/bin/env python3

"""
Contains methods for verifying various cached input data.
"""

__author__ = "Scott Teresi"

import os
import sys
import pandas as pd
from transposon.revise_annotation import Revise_Anno
from transposon.import_filtered_genes import import_filtered_genes
from transposon.import_filtered_TEs import import_filtered_TEs


def verify_chromosome_h5_cache(
        gene_data_obj,
        te_data_obj,
        h5_g_filename,
        h5_t_filename,
        reset_h5,
        h5_cache_location,
        genes_input_file,
        tes_input_file,
        chrom_id,
        logger,
):
    """Determine whether or not previously saved gene_data and TransposonData
    exist in H5 format. Each h5 file represents either gene_data or
    TransposonData for one chromosome at a time. Save/Update H5 files as
    necessary.

    When are H5 files written?
    1. H5 files will be written if there are no current corresponding H5 files
    saved on disk.
    2. If a command-line option is passed to density.py to re-write the H5
    files, this option defaults to not re-write.
    3. If enough time has passed between the creation of the H5 file and the
    current run-time of the program. TODO talk to Michael more about this.

    Args:
        gene_data_obj (gene_data): Instance of gene_data
        te_data_obj (TransposonData): Instance of TransposonData
        h5_g_filename (str): The string of the filename in which to save the
            gene_data as an H5 file.
        h5_t_filename (str): The string of the filename in which to save the
            TransposonData as an H5 file.
        reset_h5 (bool): Boolean, whether or not to completely rewrite the
            cache of all H5 files. True means that we will rewrite.
        h5_cache_location (str): The location (file path) in which to store the
            h5 files. Defaults to /filtered_input_data/h5_cache
        genes_input_file (str): The file path of the file that was used to
            generate the gene_data instance.
        tes_input_file (str): The file path of the file that was used to
            generate the TransposonData instance.
        chrom_id (str): A string representation of the current chromosome. Used
        to name each H5 file.
    """
    if reset_h5:
        logger.info("overwrite: %s" % h5_g_filename)
        logger.info("overwrite: %s" % h5_t_filename)
        gene_data_obj.write(h5_g_filename)
        te_data_obj.write(h5_t_filename)

    if os.path.exists(h5_g_filename) and os.path.exists(h5_t_filename):
        gene_annot_time = os.path.getmtime(genes_input_file)
        te_annot_time = os.path.getmtime(tes_input_file)
        gene_h5_time = os.path.getmtime(h5_g_filename)
        te_h5_time = os.path.getmtime(h5_t_filename)

        if (gene_annot_time > gene_h5_time) and (te_annot_time > te_h5_time):
            logger.info("cache is too old for chromosome '%s'" % chrom_id)
            logger.info("write: %s" % h5_g_filename)
            logger.info("write: %s" % h5_t_filename)
            gene_data_obj.write(h5_g_filename)
            te_data_obj.write(h5_t_filename)

        elif (gene_annot_time < gene_h5_time) and (te_annot_time < te_h5_time):
            # No need to re-write a current cache
            return

    elif reset_h5 or (
        not (os.path.exists(h5_g_filename) and os.path.exists(h5_t_filename))
    ):
        gene_data_obj.write(h5_g_filename)
        te_data_obj.write(h5_t_filename)
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

        revised_cache_loc ():

        logger ():

        genome_id (str): String of the genome ID

    Returns:
        te_data (pandaframe): A pandas dataframe of the TE data
    """

    if os.path.exists(revised_transposons_loc) and not revise_anno:
        logger.info("load revised TE: %s" % revised_transposons_loc)
        te_data = import_filtered_TEs(revised_transposons_loc, logger)
    else:
        logger.info("creating revised TE dataset...")
        logger.info("revising the TE dataset will take a long time!")
        # N.B we want higher recursion limit for the code
        sys.setrecursionlimit(11 ** 6)
        revised_te_data = Revise_Anno(te_data, revised_cache_loc, genome_id)
        revised_te_data.create_superfam()
        revised_te_data.create_order()
        revised_te_data.create_nameless()
        logger.info("write revised TE: %s" % revised_transposons_loc)
        revised_te_data.save_updated_te_annotation(revised_transposons_loc)
        te_data = revised_te_data.whole_te_annotation
    return te_data
