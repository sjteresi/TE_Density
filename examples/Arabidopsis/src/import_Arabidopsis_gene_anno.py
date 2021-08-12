"""
Filter a gene annotation file for the TE Density algorithm
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs

from transposon import check_nulls


def write_cleaned_genes(gene_pandaframe, output_dir, old_filename, logger):
    file_name = os.path.join(
        output_dir,
        ("Cleaned_" + os.path.splitext(os.path.basename(old_filename))[0]) + ".tsv",
    )  # MAGIC to get proper extension

    logger.info("Writing cleaned gene file to: %s" % file_name)
    gene_pandaframe.to_csv(file_name, sep="\t", header=True, index=True)


def import_genes(genes_input_path):
    """Import genes file.

    Args:
        input_dir (command line argument) Specify the input directory of the gene
        annotation data, this is the same as the TE annotation directory
    """

    col_names = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Score",
        "Strand",
        "Frame",
        "FullName",
    ]

    col_to_use = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Strand",
        "FullName",
    ]

    gene_data = pd.read_csv(
        genes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        dtype={
            "Stop": "float64",
            "Start": "float64",
            "Chromosome": str,
            "Strand": str,
            "Fullname": str,
            "Feature": str,
            "Software": str,
        },
        comment="#",
    )

    # rows in annotation
    gene_data = gene_data[gene_data.Feature == "gene"]  # drop non-gene rows
    # gene_data.reset_index(inplace=True)  # reset index so we can have proper

    gene_data["Gene_Name"] = gene_data["FullName"].str.extract(r"ID=(.*?);")
    gene_data = gene_data.drop(columns=["FullName", "Software"])
    gene_data["Length"] = gene_data.Stop - gene_data.Start + 1

    gene_data.sort_values(by=["Chromosome", "Start"], inplace=True)
    check_nulls(gene_data, logger)
    gene_data = drop_nulls(gene_data, logger)

    # Set the gene name as the index
    gene_data.set_index("Gene_Name", inplace=True)

    print(gene_data)
    return gene_data


def get_nulls(my_df, logger):
    """
    Print out the row IDs where the null values exist

    Args:
        my_df (Pandaframes): Pandaframe to check null values in
        logger
    """
    nas = my_df[my_df.isna().any(axis=1)]
    logger.warning("Rows where null exist: %s" % nas)


def drop_nulls(my_df, logger):
    """
    Drop null values inside a Pandaframe

    Args:
        my_df (Pandaframes): Pandaframe to drop null values
    """
    nas = my_df[my_df.isna().any(axis=1)]
    if not nas.empty:
        logger.warning("Dropping rows with at least one Null value!")
        my_df = my_df.dropna(axis=0, how="any")
    return my_df


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Reformat gene annotation file")
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.join(
        dir_main, "../../../../", "TE_Data/filtered_input_data"
    )
    parser.add_argument(
        "gene_input_file", type=str, help="Parent path of gene annotation file"
    )

    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=output_default,
        help="Parent directory to output results",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.gene_input_file = os.path.abspath(args.gene_input_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Execute
    cleaned_genes = import_genes(args.gene_input_file, logger)
    write_cleaned_genes(cleaned_genes, args.output_dir, args.gene_input_file, logger)
