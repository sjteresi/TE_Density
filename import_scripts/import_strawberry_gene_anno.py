"""
Filter a gene annotation file for the TE Density algorithm
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs


def write_cleaned_genes(gene_pandaframe, output_dir, genome_name, logger):
    file_name = os.path.join(
        output_dir, ("Cleaned_" + genome_name + "_gene_annotation.tsv")
    )

    logger.info("Writing cleaned gene file to: %s" % file_name)
    gene_pandaframe.to_csv(file_name, sep="\t", header=True, index=True)


def import_genes(genes_input_path, contig_del=False):
    """Import genes file.

    Args:
        input_dir (command line argument) Specify the input directory of the gene
        annotation data, this is the same as the TE annotation directory

        contig_drop (bool): logical whether to drop rows with a contig as the
        chromosome id
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

    Gene_Data = pd.read_csv(
        genes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        dtype={"Start": "float64", "Stop": "float64"},
        comment="#",
    )

    # rows in annotation
    Gene_Data = Gene_Data[Gene_Data.Feature == "gene"]  # drop non-gene rows

    # clean the names and set as the index (get row wrt name c.f. idx)

    Gene_Data["Gene_Name"] = Gene_Data["FullName"].str.extract(r"ID=(.*?);")

    Gene_Data.set_index("Gene_Name", inplace=True)
    Gene_Data = Gene_Data.drop(columns=["FullName", "Software"])

    Gene_Data.Strand = Gene_Data.Strand.astype(str)
    Gene_Data["Length"] = Gene_Data.Stop - Gene_Data.Start + 1

    if contig_del:
        Gene_Data = Gene_Data[~Gene_Data.Chromosome.str.contains("contig", case=False)]

    Gene_Data.sort_values(by=["Chromosome", "Start"], inplace=True)
    return Gene_Data


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

    write_cleaned_genes(cleaned_genes, args.output_dir, "Strawberry", logger)
