#!/usr/bin/env/python

"""
Retrieve info of a list of genes
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
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(
        description="""output to std_out information of each gene in a
        user-provided list of genes using the DensityData.info_of_gene()
        function"""
    )

    parser.add_argument(
        "te_density_hdf5_result",
        type=str,
        help="""Path to the HDF5 file that contains the genes that the user wants
        to extract TE density data from""",
    )

    parser.add_argument(
        "cleaned_gene_data_file",
        type=str,
        help="""Path to the cleaned gene data
        file that was produced prior to running the pipeline, this is necessary
        to initialize the DensityData obj""",
    )

    parser.add_argument(
        "chromosome_to_subset_gene_data",
        type=str,
        help="""The cleaned gene data
        file may contain information of genes from multiple chromosomes,
        because the TE Density data corresponds to one chromosome, please
        specify the appropriate chromosome identifier for the density data you
        are trying to access, so that we may appropriately subset the gene
        data""",
    )

    parser.add_argument(
        "genome_id",
        type=str,
        help="""Please specify the genome ID, use the same one you used as an
        argument when running TE Density""",
    )

    parser.add_argument(
        "list_of_genes",
        type=str,
        help="Path to list of genes files",
    )

    parser.add_argument(
        "--window_idx",
        type=int,
        default=1,
        help="""Index of the window that you want to access information from,
        windows go from lowest to highest, default TE density settings yield 20
        windows (500, 10000, 500) (start, stop, step).""",
    )

    parser.add_argument(
        "--n_te_types",
        type=int,
        default=5,
        help="""Number of TE types that you want to display when showing the
        top and bottom TE categories for density relative to your gene of
        interest""",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.te_density_hdf5_result = os.path.abspath(args.te_density_hdf5_result)
    args.cleaned_gene_data_file = os.path.abspath(args.cleaned_gene_data_file)
    args.list_of_genes = os.path.abspath(args.list_of_genes)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Read pandas dataframe from cleaned genes file
    full_genome_gene_data = import_filtered_genes(args.cleaned_gene_data_file, logger)

    # MAGIC, chromosome category is inherent to the pandas data frame
    full_genome_gene_data_subsetted = full_genome_gene_data.loc[
        full_genome_gene_data["Chromosome"] == args.chromosome_to_subset_gene_data
    ]

    # Initialize GeneData obj
    specific_chromosome_gene_data = GeneData(
        full_genome_gene_data_subsetted, args.genome_id
    )

    # Initialize DensityData obj
    processed_density_data = DensityData.verify_h5_cache(
        args.te_density_hdf5_result, specific_chromosome_gene_data, logger
    )

    # Read list of genes file from user and begin printing output to console
    with open(args.list_of_genes, "r", encoding="utf-8") as in_file:
        all_genes = [gene.strip() for gene in in_file]
    for gene in all_genes:
        print(
            processed_density_data.info_of_gene(gene, args.window_idx, args.n_te_types)
        )
