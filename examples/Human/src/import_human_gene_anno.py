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


def write_cleaned_genes(gene_pandaframe, output_dir, genome_name, logger):
    file_name = os.path.join(
        output_dir, ("Cleaned_Chr7_13_" + genome_name + "_Genes.tsv")
    )

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
    gene_data["Gene_Name"] = gene_data["FullName"].str.extract(r";gene_name=(.*?);")

    # NOTE
    gene_data.drop_duplicates(subset=["Gene_Name"], keep=False, inplace=True)
    # Drop duplicate gene  names for the human data set, will add explicit
    # function to fix duplicate gene names in future so that the code doesn't
    # crash (Crashes TE_Density due to shape error). TODO add error handling
    # and elegant renaming function to fix duplicate gene names, not relevant
    # to this example though

    gene_data = gene_data.drop(columns=["FullName", "Software"])

    gene_data.Strand = gene_data.Strand.astype(str)

    gene_data["Length"] = gene_data.Stop - gene_data.Start + 1

    # MAGIC I only want the 7th and 13th chromosome
    chromosomes_i_want = ["chr7", "chr13"]
    gene_data = gene_data.loc[gene_data["Chromosome"].isin(chromosomes_i_want)]

    gene_data.sort_values(by=["Chromosome", "Start"], inplace=True)
    check_nulls(gene_data, logger)
    gene_data = drop_nulls(gene_data, logger)
    gene_data.set_index("Gene_Name", inplace=True)

    return gene_data


def drop_nulls(my_df, logger):
    """
    Drop null values inside a PandaFrame

    Args:
        my_df (Pandas.Data.Frame):

    Returns:
        my_df (Pandas.Data.Frame): Without the offending rows containing null
        values, also has a call to the logger to show the user which rows were
        gotten rid of

    """
    nas = my_df[my_df.isna().any(axis=1)]
    if not nas.empty:
        logger.warning("Dropping rows with at least one Null value!")
        my_df = my_df.dropna(axis=0, how="any")
    return my_df


def return_duplicate_indices(dataframe):
    return dataframe[dataframe.index.duplicated(keep=False)]


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
    write_cleaned_genes(cleaned_genes, args.output_dir, "Human", logger)
