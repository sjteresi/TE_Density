#!/usr/bin/env python

"""
Reformat the human TE file into a format conducive to the Transposon_Data class
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import logging
import os
import coloredlogs

from examples.Human.src.replace_human_TE_names import te_annot_renamer


def write_cleaned_TEs(te_pandaframe, output_dir, genome_name, logger):
    file_name = os.path.join(
        output_dir, ("Cleaned_Chr7_13_" + genome_name + "_TEs.tsv")
    )

    logger.info("Writing cleaned TE file to: %s" % file_name)
    te_pandaframe.to_csv(file_name, sep="\t", header=True, index=False)


def import_human_tes(human_TE_file, logger):
    """
    We want a pandas object with Chromosome, Start, Stop, Strand, Order, SuperFamily and Length as columns
    """
    data = pd.read_csv(
        human_TE_file,
        header="infer",
        sep="\t",
        dtype={
            "genoStart": "float64",
            "genoEnd": "float64",
            "strand": str,
            "genoName": str,
            "repClass": str,
            "repFamily": str,
        },
    )
    data.drop(
        columns=[
            "#bin",
            "swScore",
            "milliDiv",
            "milliIns",
            "milliDel",
            "id",
            "genoLeft",
            "repStart",
            "repEnd",
            "repLeft",
            "repName",
        ],
        inplace=True,
    )
    data.rename(
        columns={
            "genoName": "Chromosome",
            "repClass": "Order",
            "repFamily": "SuperFamily",
            "strand": "Strand",
            "genoStart": "Start",
            "genoEnd": "Stop",
        },
        inplace=True,
    )

    data["Length"] = data.Stop - data.Start + 1
    # NOTE only grabbing specific chromosomes
    chromosomes_i_want = ["chr7", "chr13"]  # MAGIC
    data = data.loc[data["Chromosome"].isin(chromosomes_i_want)]

    data = te_annot_renamer(data)

    return data


if __name__ == "__main__":
    """Command line interface to calculate density."""

    parser = argparse.ArgumentParser(description="Reformat TE annotation file")
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.join(
        dir_main, "../../../../", "TE_Data/filtered_input_data"
    )
    parser.add_argument(
        "TE_input_file", type=str, help="Parent path of TE annotation file"
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
    args.TE_input_file = os.path.abspath(args.TE_input_file)
    args.output_dir = os.path.abspath(args.output_dir)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    cleaned_tes = import_human_tes(args.TE_input_file, logger)
    write_cleaned_TEs(cleaned_tes, args.output_dir, "Human", logger)
