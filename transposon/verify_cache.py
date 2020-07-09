#!/usr/bin/env python3

"""
Contains methods for verifying various cached input data.
"""

__author__ = "Scott Teresi"

import os
import pandas as pd
from transposon.import_genes import import_genes
from transposon.import_transposons import import_transposons


def verify_chromosome_h5_cache(GeneData_obj, TE_Data_obj, h5_g_filename,
                               h5_t_filename, reset_h5, h5_cache_location,
                               genes_input_file, tes_input_file, chrom_id, logger):
    """Determine whether or not previously saved GeneData and TransposonData
    exist in H5 format. Each h5 file represents either GeneData or
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
        GeneData_obj (GeneData): Instance of GeneData
        TE_Data_obj (TransposonData): Instance of TransposonData
        h5_g_filename (str): The string of the filename in which to save the
            GeneData as an H5 file.
        h5_t_filename (str): The string of the filename in which to save the
            TransposonData as an H5 file.
        reset_h5 (bool): Boolean, whether or not to completely rewrite the
            cache of all H5 files. True means that we will rewrite.
        h5_cache_location (str): The location (file path) in which to store the
            h5 files. Defaults to /filtered_input_data/h5_cache
        genes_input_file (str): The file path of the file that was used to
            generate the GeneData instance.
        tes_input_file (str): The file path of the file that was used to
            generate the TransposonData instance.
        chrom_id (str): A string representation of the current chromosome. Used
        to name each H5 file.
    """
    logger.info(reset_h5)
    if reset_h5 == True:
        logger.info("Writing h5 cache anew because of command-line arg")
        GeneData_obj.write(h5_g_filename)
        TE_Data_obj.write(h5_t_filename)

    if os.path.exists(h5_g_filename) and os.path.exists(h5_t_filename):
        gene_annot_time = os.path.getmtime(genes_input_file)
        te_annot_time = os.path.getmtime(tes_input_file)
        gene_h5_time = os.path.getmtime(h5_g_filename)
        te_h5_time = os.path.getmtime(h5_t_filename)

        if (gene_annot_time > gene_h5_time) and (te_annot_time > te_h5_time):
            logger.info("""Chromosome: %s: Writing Gene_Data and TE_Data to disk
                        in H5 format because annotation file is newer than h5
                        cache.""", chrom_id)
            GeneData_obj.write(h5_g_filename)
            TE_Data_obj.write(h5_t_filename)

        elif (gene_annot_time < gene_h5_time) and (te_annot_time < te_h5_time):
            # No need to re-write a current cache
            return

    elif reset_h5 or (not(os.path.exists(h5_g_filename) and
                          os.path.exists(h5_t_filename))):
        logger.info("""Chromosome: %s: Writing Gene_Data and TE_Data
                    to disk in H5 format""", chrom_id)
        GeneData_obj.write(h5_g_filename)
        TE_Data_obj.write(h5_t_filename)
    else:
        logger.warning('''During the verification of the H5 cache nothing was
                       saved because 0 conditions were met.''')


def verify_TE_cache(tes_input_file, cleaned_transposons, te_annot_renamer,
                    contig_del, logger):
    """Determine whether or not previously filtered TE data exists, if it
    does, read it from disk. If it does not, read the raw annotation file and
    make a filtered dataset for future import.

    Args:
        tes_input_file (str): A command line argument, this is the location
            of the TE annotation file.

        cleaned_tranposons (str): A string representing the path of a previously
            filtered TE file via import_transposons().

        te_annot_renamer (function containing a dictionary and other methods):
            imported from separate file within the repository. This file
            performs the more specific filtering steps on the TEs such as
            changing the annotation details for specific TE types.

        contig_del (bool): A boolean of whether to remove contigs on import

    Returns:
        TE_Data (pandaframe): A pandas dataframe of the TE data
    """
    logger.info("Verifying TransposonData cache...")
    if os.path.exists(cleaned_transposons):
        te_annot_time = os.path.getmtime(tes_input_file)
        cleaned_te_time = os.path.getmtime(cleaned_transposons)
        if te_annot_time > cleaned_te_time:
            logger.info("""TE annotation file is newer than the previously
                        filtered data set. Importing TE data from the
                        annotation file and re-writing the filtered input
                        data""")
            TE_Data = import_transposons(tes_input_file, te_annot_renamer,
                                         contig_del)
            TE_Data.to_csv(cleaned_transposons, sep='\t', header=True, index=False)
        else:
            logger.info("Importing filtered transposons from disk...")
            TE_Data = pd.read_csv(cleaned_transposons, header='infer',
                                  dtype={'Start': 'float32', 'Stop': 'float32',
                                         'Length': 'float32'}, sep='\t')
    else:
        logger.info("Previously filtered TE dataset does not exist...")
        logger.info("Importing unfiltered TE dataset from annotation file...")
        TE_Data = import_transposons(tes_input_file, te_annot_renamer,
                                     contig_del)
        TE_Data.to_csv(cleaned_transposons, sep='\t', header=True, index=False)
    return TE_Data


def verify_gene_cache(genes_input_file, cleaned_genes, contig_del, logger):
    """Determine whether or not previously filtered gene data exists, if it
    does, read it from disk. If it does not, read the raw annotation file and
    make a filtered dataset for future import.

    Args:
        genes_input_file (str): A command line argument, this is the location
            of the gene annotation file.

        cleaned_genes (str): A string representing the path of a previously
            filtered gene file via import_genes().

        contig_del (bool): A boolean of whether to remove contigs on import

    Returns:
        Gene_Data (pandaframe): A pandas dataframe of the Gene data
    """
    logger.info("Verifying GeneData cache...")
    if os.path.exists(cleaned_genes):
        gene_annot_time = os.path.getmtime(genes_input_file)
        cleaned_gene_time = os.path.getmtime(cleaned_genes)
        if gene_annot_time > cleaned_gene_time:
            logger.info("""Gene annotation file is newer than the previously
                        filtered data set. Importing gene data from the
                        annotation file and re-writing the filtered input
                        data""")
            Gene_Data = import_genes(genes_input_file, contig_del)
            Gene_Data.to_csv(cleaned_genes, sep='\t', header=True, index=True)

        else:
            logger.info("Importing filtered gene dataset from disk...")
            Gene_Data = pd.read_csv(cleaned_genes, header='infer', sep='\t',
                                    dtype={'Start': 'float32', 'Stop': 'float32',
                                           'Length': 'float32'}, index_col='Gene_Name')
    else:
        logger.info("Previously filtered gene dataset does not exist...")
        logger.info("Importing unfiltered gene dataset from annotation file...")
        Gene_Data = import_genes(genes_input_file, contig_del)
        Gene_Data.to_csv(cleaned_genes, sep='\t', header=True, index=True)
    return Gene_Data
