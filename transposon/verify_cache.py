#!/usr/bin/env python3

"""
Contains methods for verifying various cached input data.
"""

__author__ = "Scott Teresi"

import os
import pandas as pd
from transposon.import_genes import import_genes
from transposon.import_transposons import import_transposons
from transposon.revise_annotation import Revise_Anno


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
    2. If a command-line option is passed to process.py to re-write the H5
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
        logger.info("Writing h5 cache anew because of command-line arg")
        gene_data_obj.write(h5_g_filename)
        te_data_obj.write(h5_t_filename)

    if os.path.exists(h5_g_filename) and os.path.exists(h5_t_filename):
        gene_annot_time = os.path.getmtime(genes_input_file)
        te_annot_time = os.path.getmtime(tes_input_file)
        gene_h5_time = os.path.getmtime(h5_g_filename)
        te_h5_time = os.path.getmtime(h5_t_filename)

        if (gene_annot_time > gene_h5_time) and (te_annot_time > te_h5_time):
            logger.info(
                """Chromosome: %s: Writing gene_data and te_data to disk
                        in H5 format because annotation file is newer than h5
                        cache.""",
                chrom_id,
            )
            gene_data_obj.write(h5_g_filename)
            te_data_obj.write(h5_t_filename)

        elif (gene_annot_time < gene_h5_time) and (te_annot_time < te_h5_time):
            # No need to re-write a current cache
            return

    elif reset_h5 or (
        not (os.path.exists(h5_g_filename) and os.path.exists(h5_t_filename))
    ):
        logger.info(
            """Chromosome: %s: Writing gene_data and te_data
                    to disk in H5 format""",
            chrom_id,
        )
        gene_data_obj.write(h5_g_filename)
        te_data_obj.write(h5_t_filename)
    else:
        logger.warning(
            """During the verification of the H5 cache nothing was
                       saved because 0 conditions were met."""
        )


def verify_TE_cache(
    tes_input_file, cleaned_transposons, te_annot_renamer, contig_del, logger
):
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
        te_data (pandaframe): A pandas dataframe of the TE data
    """
    logger.info("Verifying TransposonData cache...")
    if os.path.exists(cleaned_transposons):
        te_annot_time = os.path.getmtime(tes_input_file)
        cleaned_te_time = os.path.getmtime(cleaned_transposons)
        if te_annot_time > cleaned_te_time:
            logger.info(
                """TE annotation file is newer than the previously
                        filtered data set. Importing TE data from the
                        annotation file and re-writing the filtered input
                        data"""
            )
            te_data = import_transposons(tes_input_file, te_annot_renamer, contig_del)
            te_data.sort_values(by=["Chromosome", "Start"], inplace=True)
            te_data.to_csv(cleaned_transposons, sep="\t", header=True, index=False)
        else:
            logger.info("Importing filtered transposons from disk...")
            te_data = pd.read_csv(
                cleaned_transposons,
                header="infer",
                dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
                sep="\t",
            )
    else:
        logger.info("Previously filtered TE dataset does not exist...")
        logger.info("Importing unfiltered TE dataset from annotation file...")
        te_data = import_transposons(tes_input_file, te_annot_renamer, contig_del)
        te_data.sort_values(by=["Chromosome", "Start"], inplace=True)
        te_data.to_csv(cleaned_transposons, sep="\t", header=True, index=False)
    return te_data


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
        gene_data (pandaframe): A pandas dataframe of the Gene data
    """
    logger.info("Verifying gene_data cache...")
    if os.path.exists(cleaned_genes):
        gene_annot_time = os.path.getmtime(genes_input_file)
        cleaned_gene_time = os.path.getmtime(cleaned_genes)
        if gene_annot_time > cleaned_gene_time:
            logger.info(
                """Gene annotation file is newer than the previously
                        filtered data set. Importing gene data from the
                        annotation file and re-writing the filtered input
                        data"""
            )
            gene_data = import_genes(genes_input_file, contig_del)
            gene_data.sort_values(by=["Chromosome", "Start"], inplace=True)
            gene_data.to_csv(cleaned_genes, sep="\t", header=True, index=True)

        else:
            logger.info("Importing filtered gene dataset from disk...")
            gene_data = pd.read_csv(
                cleaned_genes,
                header="infer",
                sep="\t",
                dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
                index_col="Gene_Name",
            )
    else:
        logger.info("Previously filtered gene dataset does not exist...")
        logger.info("Importing unfiltered gene dataset from annotation file...")
        gene_data = import_genes(genes_input_file, contig_del)
        gene_data.sort_values(by=["Chromosome", "Start"], inplace=True)
        gene_data.to_csv(cleaned_genes, sep="\t", header=True, index=True)
    return gene_data


def revise_annotation(
    TE_Data,
    create_revise_anno,
    dont_use_revise_anno,
    revised_transposons_loc,
    revised_cache_loc,
    logger,
    genome_id,
):
    """
    Revises the annotation so that elements of the same type do not overlap at
    all. Will essentially merge elements together, elongating them. This is
    done so that the mathematics of density make sense. You can elect to not
    use the revised annotation through a command-line argument to process.py,
    however given that TEs often overlap with one another in annotatios (not
    just being nested in one another) it can lead to some difficulties in
    accurately assessing density and obfuscate the density results.

    Args:
        TE_Data (pandas.core.DataFrame): A PandaFrame of the TE data,
        previously imported from raw and filtered or imported from a previously
        filtered data file that was saved to disk.

        create_revise_anno (bool): A boolean of whether or not to create a revised
        annotation

        dont_use_revise_anno (bool): A boolean of whether or not to use a revised
        annotation, if False, we just return the input TE annotation

        revised_transposons (str): A string representing the path of a
        previously filtered (cleaned) and revised TE annotation.

        revised_cache_loc ():

        logger ():

        genome_id (str): String of the genome ID

    Returns:
        TE_Data (pandaframe): A pandas dataframe of the TE data
    """
    if dont_use_revise_anno:
        logger.info(
            """User provided flag to NOT use a revised annotation, please
                    refer to the README for documentation on why this may not
                    be a great idea."""
        )
        return TE_Data

    else:
        logger.info("Checking to see if revised TE dataset exists...")
        if create_revise_anno:
            logger.info("Flag provided, forcing creation of revised TE dataset...")
            revised_TE_Data = Revise_Anno(TE_Data, revised_cache_loc, logger, genome_id)
            revised_TE_Data.create_superfam()
            revised_TE_Data.create_order()
            revised_TE_Data.create_nameless()
            logger.info("Saving revised TE dataset...")
            revised_TE_Data.save_whole_te_annotation(revised_transposons_loc)
            TE_Data = revised_TE_Data.whole_te_annotation

        elif os.path.exists(revised_transposons_loc):
            logger.info("Exists. Importing revised TE dataset from disk...")
            TE_Data = pd.read_csv(
                revised_transposons_loc,
                header="infer",
                dtype={"Start": "float32", "Stop": "float32", "Length": "float32"},
                sep="\t",
            )
        else:
            logger.info("Previously revised TE dataset does not exist...")
            logger.info("Creating revised TE dataset...")
            revised_TE_Data = Revise_Anno(TE_Data, revised_cache_loc, logger, genome_id)
            revised_TE_Data.create_superfam()
            revised_TE_Data.create_order()
            revised_TE_Data.create_nameless()
            logger.info("Saving revised TE dataset...")
            revised_TE_Data.save_updated_te_annotation(revised_transposons_loc)
            TE_Data = revised_TE_Data.whole_te_annotation
        return TE_Data