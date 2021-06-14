#!/usr/bin/env/python

"""
Compare TE Density values between orthologs
"""

__author__ = "Scott Teresi"

import argparse
import os
import logging
import coloredlogs
import pandas as pd
import numpy as np
from collections import namedtuple

import time

from transposon.gene_data import GeneData
from transposon.density_data import DensityData
from transposon.import_filtered_genes import import_filtered_genes
from examples.Rice_Synteny.src.bargraphs import graph_barplot_density_differences


def read_ortholog_table(rice_ortholog_table):
    """
    Read a pandaframe from disk

    Args:
        rice_ortholog_table (str): Path to pandaframe saved on disk which is in
        .tsv format with the columns specified below

    Returns:
        dataframe (pandas.DataFrame): Pandas dataframe of the file that was
        provided
    """
    dataframe = pd.read_csv(
        rice_ortholog_table,
        sep="\t",
        header="infer",
        dtype={
            "OrgA_Chromosome": str,
            "Glaberrima": str,
            "OrgB_Chromosome": str,
            "Sativa": str,
            "E_Value": "float64",
        },
        na_values=["NA"],
    )
    return dataframe


def save_ortholog_table(ortholog_table, filename, output_dir, logger):
    """
    Saves an ortholog table (pandas.DataFrame) to disk

    Args:
        ortholog_table (pandas.DataFrame): A pandas dataframe of ortholog
        relationships between genes

        filename (str): String representing the name of the panda frame the
        user intends to save

        output_dir (str): The path to the output directory

        logger (logging.Logger): Object to pass logging commands through

    Returns:
        None
    """
    logger.info("Saving data to: %s" % (os.path.join(output_dir, filename)))
    ortholog_table.to_csv(
        os.path.join(output_dir, filename), sep="\t", header="True", index=False
    )


def identify_indices_of_syntelogs(
    orthologs, sativa_dd_obj, glaberrima_dd_obj, output_dir, logger
):
    """
    Add the indices for each gene in the ortholog table as a column to the
    ortholog table. Also saves this processed version to disk.

    Args:
        orthologs (pandas.DataFrame): A pandas dataframe of the ortholog
        relatioships between genes

        sativa_dd_obj (DensityData): An instance of DensityData in this case
        representing the information for O. sativa.

        glaberrima_dd_obj (DensityData): An instance of DensityData in this case
        representing the information for O. glaberrima.

        output_dir (str): The path to the output directory

        logger (logging.Logger): Object to pass logging commands through


    Returns:
        syntelog_matches (pandas.DataFrame): A purified version of the input
        argument 'orthologs' that only contains entries where there is a valid
        value locating a syntelog pair in each of their respective TE density
        files. No null or missing values.
    """
    # NB MAGIC the genome names / column names are hard-coded
    glaberrima_orth_genes = orthologs["Glaberrima"].tolist()
    sativa_orth_genes = orthologs["Sativa"].tolist()

    Index_and_Chromosome = namedtuple("Index_and_Chromosome", "Index Chromosome")

    super_dict = {}
    sativa_indices_and_genes = {
        dd_gene: Index_and_Chromosome(dd_index, sativa_dd_obj.unique_chromosome_id)
        for dd_index, dd_gene in enumerate(sativa_dd_obj.gene_list)
    }
    super_dict.update(sativa_indices_and_genes)
    orthologs["Sativa_Indices"] = [
        super_dict.get(gene_name, None) for gene_name in sativa_orth_genes
    ]

    super_dict = {}
    glaberrima_indices_and_genes = {
        dd_gene: Index_and_Chromosome(dd_index, glaberrima_dd_obj.unique_chromosome_id)
        for dd_index, dd_gene in enumerate(glaberrima_dd_obj.gene_list)
    }
    super_dict.update(glaberrima_indices_and_genes)
    orthologs["Glaberrima_Indices"] = [
        super_dict.get(gene_name, None) for gene_name in glaberrima_orth_genes
    ]

    syntelog_matches = orthologs[
        orthologs[["Sativa_Indices", "Glaberrima_Indices"]].notnull().all(1)
    ]  # the all command wants them both to be true

    # MAGIC name for file
    save_ortholog_table(syntelog_matches, "Syntelog_Matches", output_dir, logger)
    return syntelog_matches


def verify_te_type_sequence(density_data_one, density_data_two):
    """
    Simple function to verify that the actual order of the TE types in the
    genomes are the same, as this would cause an error during the comparison
    stage if they had different indices in the HDF5. This could be due to
    different TE annotation sources but should not be an issue for this
    example. I am including this note for future work as it could be an easy
    thing to overlook.

    Args:
        density_data_one (DensityData): A density data instance representing a
        single genome

        density_data_two (DensityData): A density data instance representing a
        single genome

    Returns: None, raises errors if detects errors as described above
    """
    if density_data_one.order_list != density_data_two.order_list:
        raise ValueError(
            """The TE Orders in %s do not match the TE Orders in
                         %s"""
            % (density_data_one, density_data_two)
        )
    if density_data_one.super_list != density_data_two.super_list:
        raise ValueError(
            """The TE Superfamilies in %s do not match the TE Superfamilies in
                         %s"""
            % (density_data_one, density_data_two)
        )


def verify_window_sequence(density_data_one, density_data_two):
    """
    Simple function to verify that the actual order of the window values in the
    genomes are the same, as this would cause an error during the comparison
    stage if they had different indices in the HDF5. This could be due to
    running the program with different window settings. I am including this
    note for future work as it could be an easy thing to overlook.

    Args:
        density_data_one (DensityData): A density data instance representing a
        single genome

        density_data_two (DensityData): A density data instance representing a
        single genome

    Returns: None, raises errors if detects errors as described above
    """
    if density_data_one.window_list != density_data_two.window_list:
        raise ValueError(
            """The windows in %s do not match the windows in
                         %s"""
            % (density_data_one, density_data_two)
        )


def verify_te_index_dictionary_sequence(density_data_one, density_data_two):
    """


    This could be due to
    running the program with different window settings. I am including this
    note for future work as it could be an easy thing to overlook.

    Args:
        density_data_one (DensityData): A density data instance representing a
        single genome

        density_data_two (DensityData): A density data instance representing a
        single genome

    Returns: None, raises errors if detects errors as described above
    """
    if density_data_one.order_index_dict != density_data_two.order_index_dict:
        raise ValueError(
            """The order index dictionary in %s do not match the order index
            dictionary in %s"""
            % (density_data_one, density_data_two)
        )

    if density_data_one.super_index_dict != density_data_two.super_index_dict:
        raise ValueError(
            """The super index dictionary in %s do not match the super index
            dictionary in %s"""
            % (density_data_one, density_data_two)
        )


def graph_histogram_syntelog_differences(
    syntelog_table_w_indices,
    sativa_dd_obj,
    glaberrima_dd_obj,
    te_type_iterable,
    te_hdf5_grouping,
    te_index_dict,
    direction,
    output_dir,
    logger,
):
    """

    Args:
        syntelog_table_w_indices (pandas.DataFrame): A purified version of the
        orthology table that only contains entries where there is a valid
        value locating a syntelog pair in each of their respective TE density
        files. No null or missing values. Intended to be received from the
        function 'identify_indices_of_syntelogs'

        sativa_dd_obj (DensityData): A density data instance representing the
        information for O. sativa

        glaberrima_dd_obj (DensityData): A density data instance representing the
        information for O. glaberrima

        te_type_iterable (list): Iterable of TE types as string

        te_hdf5_grouping (str): String representing the section of the HDF5 we
        need to access, acceptable arguments are 'RHO_ORDERS_LEFT',
        'RHO_ORDERS_RIGHT', 'RHO_SUPERFAMILIES_LEFT', and
        'RHO_SUPERFAMILIES_RIGHT'

        te_index_dict (dictionary): Dictionary of the te types and their
        indices in the HDF5

        direction (str): String representing whether or not the dat represents
        the upstream or downstrem values

        output_dir (str): Path to string for output directory

        logger (logging.Logger): Object to pass logging commands through

    Returns: None, calls the graphing code with the appropriate input data
    """
    sativa_indices = syntelog_table_w_indices["Sativa_Indices"].tolist()
    glaberrima_indices = syntelog_table_w_indices["Glaberrima_Indices"].tolist()
    verify_te_type_sequence(sativa_dd_obj, glaberrima_dd_obj)
    verify_window_sequence(sativa_dd_obj, glaberrima_dd_obj)
    verify_te_index_dictionary_sequence(sativa_dd_obj, glaberrima_dd_obj)

    # Iterate over the TEs for the given major grouping
    for given_te_type in te_type_iterable:

        for given_window_idx, given_window_val in enumerate(
            glaberrima_dd_obj.window_list  # pick one window because they both
            # match
        ):
            list_of_differences = []

            for glaberrima_val, sativa_val in zip(glaberrima_indices, sativa_indices):
                if glaberrima_val.Chromosome == sativa_val.Chromosome:

                    diff_val = (
                        glaberrima_dd_obj.data_frame[te_hdf5_grouping][
                            te_index_dict[given_te_type],
                            given_window_idx,
                            glaberrima_val.Index,
                        ]
                        - sativa_dd_obj.data_frame[te_hdf5_grouping][
                            te_index_dict[given_te_type],
                            given_window_idx,
                            sativa_val.Index,
                        ]
                    )  # NB remember shape is type, window, gene
                    list_of_differences.append(diff_val)

            # NB, MAGIC
            # Generate count of values that are 0.
            # Then remove density differences of zero from the list that we
            # use to generate the histogram, generate a count of 0 values so
            # that we can add to the legend of the plot.
            count_of_zeroes = list_of_differences.count(0)
            non_zero_differences = [i for i in list_of_differences if i != 0]

            graph_barplot_density_differences(
                non_zero_differences,
                given_te_type,
                given_window_val,
                direction,
                count_of_zeroes,
                os.path.join(output_dir, "graphs"),
                logger,
            )


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    data_dir = os.path.abspath(os.path.join(dir_main, "../../", "Domestication_Data"))
    parser = argparse.ArgumentParser(
        description="compare density values between orthologs"
    )

    parser.add_argument(
        "ortholog_input_file",
        type=str,
        help="parent path to ortholog file",
    )

    parser.add_argument(
        "sativa_chromosome_1",
        type=str,
        help="Path to Sativas first chromosome of HDF5 data",
    )

    parser.add_argument(
        "glaberrima_chromosome_1",
        type=str,
        help="Path to Glaberrimas first chromosome of HDF5 data",
    )

    parser.add_argument(
        "sativa_gene_data", type=str, help="Path to Sativas filtered gene data file"
    )

    parser.add_argument(
        "glaberrima_gene_data",
        type=str,
        help="Path to Glaberrimas filtered gene data file",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=data_dir,
        help="parent directory to output results",
    )
    args = parser.parse_args()
    args.ortholog_input_file = os.path.abspath(args.ortholog_input_file)

    args.sativa_chromosome_1 = os.path.abspath(args.sativa_chromosome_1)
    args.glaberrima_chromosome_1 = os.path.abspath(args.glaberrima_chromosome_1)

    args.sativa_gene_data = os.path.abspath(args.sativa_gene_data)
    args.glaberrima_gene_data = os.path.abspath(args.glaberrima_gene_data)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Begin reading files:
    logger.info("Reading ortholog input file: %s" % (args.ortholog_input_file))
    orthologs = read_ortholog_table(args.ortholog_input_file)

    orthologs = orthologs.loc[
        (orthologs["OrgA_Chromosome"] == "1") & (orthologs["OrgB_Chromosome"] == "1")
    ]

    logger.info(
        "Reading gene datas %s & %s"
        % (args.sativa_gene_data, args.glaberrima_gene_data)
    )
    all_genes_sativa = import_filtered_genes(args.sativa_gene_data, logger)
    all_genes_glaberrima = import_filtered_genes(args.glaberrima_gene_data, logger)

    # MAGIC filter the gene pandaframes to only have one chromosome because we
    # are working with only one chromosome of HDF5
    chrom1_sativa = all_genes_sativa.loc[all_genes_sativa["Chromosome"] == "1"]
    chrom1_glaberrima = all_genes_glaberrima.loc[
        all_genes_glaberrima["Chromosome"] == "1"
    ]

    # Wrap the pandafames as GeneData
    sativa_gene_data = GeneData(chrom1_sativa, "Sativa")
    glaberrima_gene_data = GeneData(chrom1_glaberrima, "Glaberrima")

    # NB rename for clarity
    sativa_hdf5 = args.sativa_chromosome_1
    glaberrima_hdf5 = args.glaberrima_chromosome_1

    processed_sativa_density_data = DensityData.verify_h5_cache(
        sativa_hdf5, sativa_gene_data, logger
    )
    processed_glaberrima_density_data = DensityData.verify_h5_cache(
        glaberrima_hdf5, glaberrima_gene_data, logger
    )

    syntelog_table_w_indices = identify_indices_of_syntelogs(
        orthologs,
        processed_sativa_density_data,
        processed_glaberrima_density_data,
        args.output_dir,
        logger,
    )

    graph_histogram_syntelog_differences(
        syntelog_table_w_indices,
        processed_sativa_density_data,
        processed_glaberrima_density_data,
        processed_sativa_density_data.order_list,
        "RHO_ORDERS_LEFT",
        processed_sativa_density_data.order_index_dict,
        "Upstream",
        args.output_dir,
        logger,
    )

    graph_histogram_syntelog_differences(
        syntelog_table_w_indices,
        processed_sativa_density_data,
        processed_glaberrima_density_data,
        processed_sativa_density_data.order_list,
        "RHO_ORDERS_RIGHT",
        processed_sativa_density_data.order_index_dict,
        "Downstream",
        args.output_dir,
        logger,
    )

    graph_histogram_syntelog_differences(
        syntelog_table_w_indices,
        processed_sativa_density_data,
        processed_glaberrima_density_data,
        processed_sativa_density_data.super_list,
        "RHO_SUPERFAMILIES_LEFT",
        processed_sativa_density_data.super_index_dict,
        "Upstream",
        args.output_dir,
        logger,
    )

    graph_histogram_syntelog_differences(
        syntelog_table_w_indices,
        processed_sativa_density_data,
        processed_glaberrima_density_data,
        processed_sativa_density_data.super_list,
        "RHO_SUPERFAMILIES_RIGHT",
        processed_sativa_density_data.super_index_dict,
        "Downstream",
        args.output_dir,
        logger,
    )
