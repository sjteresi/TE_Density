#!/usr/bin/env python3

"""
Preprocess input data.
"""

from transposon.replace_names import te_annot_renamer
from transposon.verify_cache import (
    verify_chromosome_h5_cache,
    verify_TE_cache,
    verify_gene_cache,
    revise_annotation,
)
from transposon.revise_annotation import Revise_Anno

# SCOTT you should probably use class here b/c there will be a set of functions w/ a
# shared state, e.g: filepaths, input files, order that they are genes_processed, etc
# you can move the functions below into bound methods
class Preprocessor:

    def __init__(self):

        pass

            # NOTE
            # MAGIC NUMBER for g_fname, and t_fname trying to get filename without
            # extension, will produce an unexpected, but functional filename if their
            # input filename has multiple . in it
            # g_fname = os.path.basename(os.path.splitext(args.genes_input_file)[0])
            # t_fname = os.path.basename(os.path.splitext(args.tes_input_file)[0])
            # # the genome input file, maybe make it user provided
            # cleaned_genes = os.path.join(
            #     args.filtered_input_data, str("Cleaned_" + g_fname + ".tsv")
            # )
            # cleaned_transposons = os.path.join(
            #     args.filtered_input_data, str("Cleaned_" + t_fname + ".tsv")
            # )
            # revised_transposons = os.path.join(
            #     args.revised_input_data, str("Revised_" + t_fname + ".tsv")
            # )
            # logging.debug("revised_transposons:  %s" % revised_transposons)
            # gene_data_unwrapped = verify_gene_cache(
            #     args.genes_input_file, cleaned_genes, args.contig_del, logger
            # )


def cache_genome_files(gene_filepath, transposon_filepath, force_refresh=False):
    """Helper to cache both gene / transposon data files."""

    raise NotImplementedError()  # TODO
    # call cache_gene_data
    # call cache_transposon_data

def cache_gene_data(gene_filepath, force_refresh=False):
    """Caches GeneData files provided the raw gene file (tsv).

    Returns:
        list(str): filepaths to the GeneData files.
    """

    raise NotImplementedError()  # TODO

        # # Revise the TE annotation
        # te_data_unwrapped = revise_annotation(
        #     te_data_unwrapped,
        #     args.revise_anno,
        #     revised_transposons,
        #     args.revised_input_data,
        #     logger,
        #     args.genome_id,
        # )

def cache_transposon_data(gene_filepath, force_refresh=False):
    """Caches TransposonData files provided the raw gene file (tsv).

    Returns:
        list(str): filepaths to the GeneData files.
    """

    raise NotImplementedError()  # TODO

        # te_data_unwrapped = verify_TE_cache(
        #     args.tes_input_file,
        #     cleaned_transposons,
        #     te_annot_renamer,
        #     args.contig_del,
        #     logger,
        # )
