#!/usr/bin/env python3

"""
Preprocess input data.
"""

import errno
import os

from transposon import raise_if_no_file, raise_if_no_dir

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
    """Prepares input data for downstream processing.

    Notable tasks include removing unwanted fields and writing into new formats.
    """

    CLEAN_PREFIX = "Cleaned_"
    REVISED_PREFIX = "Revised_"
    EXT = "tsv"

    # TODO SCOTT what is 'contig_del'? what does it mean to remove contigs? what's a contig?
    # TODO SCOTT can you reduce the number of inputs here?
    # does filtered / revised have to be different?
    def __init__(self,
    gene_file,
    transposon_file,
    filtered_dir,
    revised_dir,
    genome_id,
    revise_transposons=False
    contig_del=False,
    logger=None):

        self._logger = logger or logging.getLogger(self.__class__.__name__)
        raise_if_no_file(gene_file)
        raise_if_no_file(transposon_file)
        os.makedirs(args.filtered_input_data, exist_ok=True)
        os.makedirs(args.revised_input_data, exist_ok=True)

        self.filtered_dir = filtered_dir
        self.revised_dir = revised_dir

        self.gene_in = str(gene_file)
        self.gene_out = self._processed_filename(self.gene_in, self.filtered_dir, self.CLEAN_PREFIX)
        self.te_in = str(transposon_file)
        self.te_out = self._processed_filename(self.te_in, self.filtered_dir, self.CLEAN_PREFIX)
        self.te_revised = self._processed_filename(self.te_in, self.revised_dir, self.REVISED_PREFIX)
        self.contiguous_delete = contig_del
        self.do_transposon_revisions = revise_transposons

        self.gene_frame = None
        self.transposon_frame = None
        self.gene_files = None
        self.transposon_files = None

    def process(self):
        """Helper to execute all preprocessing tasks.

        Mutates self.
        """

        self.gene_frame = self._filter_genes()
        transposon_frame_ = self._filter_tranposons()
        self.transposon_frame = self._revise_transposon_annnotations(transposon_frame)

        raise NotImplementedError()
        self.gene_files = list()
        self.transposon_files = list()

    @property
    def data_files(self):
        """Generate paths of GeneData, TransposonData pairs."""

        raise NotImplementedError()

    @classmethod
    def _processed_filename(cls, filepath, filtered_dir, prefix):
        """Preprocessed filepath given original.

        Args:
            filepath(str): path to input file (e.g. gene or transposon)
            filtered_dir(str): directory to output cleaned files
            prefix(str): prepend to filename
        Raises:
            ValueError: bad input filetype
        Returns:
            str: expected filepath
        """

        name, ext = os.path.splitext(gene_or_te_filepath)
        if ext != "." + cls.EXT:
            raise ValueError("unkown genome file ext '%s' at %s" % (ext, filepath))
        if not os.path.isdir(filtered_dir):
        newname = prefix + name + "." + cls.EXT
        return os.path.join(filtered_dir, newname)

    def _filter_genes(self):
        """Updates filtered gene file if necessary.

        Returns:
            pandas.DataFrame: preprocessed gene frame
        """

        # NB this does not create any GeneData
        gene_data_unwrapped = verify_gene_cache(
            self.gene_in, self.gene_out, self.contiguous_delete, self._logger
        )
        return gene_data_unwrapped

    def _filter_tranposons(self):
        """Updates filtered transposon file if necessary.

        Returns:
            pandas.DataFrame: preprocessed transposon frame
        """

        # NB this does not create any TransposonData
        te_data_unwrapped = verify_TE_cache(
            self.te_in,
            self.te_out,
            te_annot_renamer,
            self.contiguous_delete,
            self._logger,
        )
        return te_data_unwrapped

    def _revise_transposon_annnotations(self, transposon_frame):

        te_data_unwrapped = revise_annotation(
            transposon_frame,
            self.do_transposon_revisions,
            self.te_revised,
            self.revised_dir
            self._logger,
            args.genome_id,
        )
        return te_data_unwrapped


    def _cache_data_files(self):
        """Update GeneData and TransposonData files on disk."""


        raise NotImplementedError()

        # somethng like this?
            # grouped_genes = split(gene_data_unwrapped, "Chromosome")
            # grouped_TEs = split(te_data_unwrapped, "Chromosome")
            # # check docstring for my split func
            # check_groupings(grouped_genes, grouped_TEs, logger, genome_id)
            # for chromosome_of_gene_data, chromosome_of_te_data in zip(
            #     grouped_genes, grouped_TEs
            # ):
                # wrapped_genedata_by_chrom = GeneData(chromosome_of_gene_data.copy(deep=True))
                # wrapped_genedata_by_chrom.add_genome_id(genome_id)
                # wrapped_tedata_by_chrom = TransposonData(chromosome_of_te_data.copy(deep=True))
                # wrapped_tedata_by_chrom.add_genome_id(genome_id)
                # chrom_id = wrapped_genedata_by_chrom.chromosome_unique_id
                #
                # # NOTE
                # # MICHAEL, these are how the H5 files are saved
                # h5_g_filename = os.path.join(h5_cache_location, str(chrom_id + "_GeneData.h5"))
                # h5_t_filename = os.path.join(h5_cache_location, str(chrom_id + "_TEData.h5"))
                #
                # verify_chromosome_h5_cache(
                #     wrapped_genedata_by_chrom,
                #     wrapped_tedata_by_chrom,
                #     h5_g_filename,
                #     h5_t_filename,
                #     reset_h5,
                #     h5_cache_location,
                #     genes_input_file,
                #     tes_input_file,
                #     chrom_id,
                #     logger,
                # )
