#!/usr/bin/env/python

"""
Barebones initialization of DensityData class intended for demonstration
purposes
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs

from transposon.gene_data import GeneData
from transposon.density_data import DensityData
from transposon.import_filtered_genes import import_filtered_genes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="generate graphs")

    parser.add_argument(
        "dummy_chromosome_hdf5_data",
        type=str,
        help="parent path to an HDF5 file of TE density data",
    )

    parser.add_argument(
        "cleaned_gene_annotation",
        type=str,
        help="""parent path to the cleaned gene data file that was derived from
        preprocessing""",
    )
    parser.add_argument(
        "dummy_chromosome_id",
        help="""ID of the chromosome of the HDF5 data, so that we may
        appropriately subset the GeneData file, which is contains information
        for all chromosomes. Reference your cleaned gene data file to get the
        appropriate string of the chromosome ID""",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    args = parser.parse_args()
    args.dummy_chromosome_hdf5_data = os.path.abspath(args.dummy_chromosome_hdf5_data)
    args.cleaned_gene_annotation = os.path.abspath(args.cleaned_gene_annotation)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Begin reading files:
    # Get the genes:
    gene_pandaframe = import_filtered_genes(args.cleaned_gene_annotation, logger)
    # Subset the pandas file to get the appropriate chromosome only
    # NOTE this can 'succeed' even if your chromosome ID is incorrect, which
    # would give you an unexpected pandaframe, make sure you check that this
    # object is what it should be after subsetting by chromosome
    gene_pandaframe = gene_pandaframe.loc[
        gene_pandaframe["Chromosome"] == args.dummy_chromosome_id
    ]

    logger.info("Initializing GeneData from cleaned annotation file")
    # Wrap as GeneData
    gene_data = GeneData(
        gene_pandaframe, str("Sample_Genome_Name_" + str(args.dummy_chromosome_id))
    )  # The last argument of GeneData is just a genome ID

    logger.info("Initializing DensityData from HDF5 file")

    # NOTE will initialize the post processed file in the same directory that
    # the original density data file is in.
    processed_one_chromosome_density_data = DensityData.verify_h5_cache(
        args.dummy_chromosome_hdf5_data, gene_data, logger
    )
