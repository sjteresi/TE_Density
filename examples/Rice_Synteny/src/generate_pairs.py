#!/usr/bin/env python3

"""
Master code file. Control filtration of syntelog data and generate summary
table
"""

__author__ = "Scott Teresi"

import argparse
import os

import logging
import coloredlogs


from import_syntelogs import import_syntelogs
from import_syntelogs import Syntelog_Data

from transposon import check_nulls


def process(
    syntelog_input_file, data_output_path,
):
    # Import the synteny data from raw file
    logger.info("Importing syntelogs: %s" % syntelog_input_file)
    syntelogs = import_syntelogs(syntelog_input_file)
    check_nulls(syntelogs, logger)

    # Wrap the data
    logger.debug("Wrapping Syntelog_Data...")
    instance_Syntelog_Data = Syntelog_Data(syntelogs)
    file_to_save = os.path.join(data_output_path, "set_syntelogs.tsv")
    logger.info("Writing syntelog data to disk: %s" % file_to_save)
    instance_Syntelog_Data.save_to_disk(file_to_save)


if __name__ == "__main__":
    """Command line interface to link syntelogs together."""

    parser = argparse.ArgumentParser(description="Filter syntelogs")
    path_main = os.path.abspath(__file__)
    parser.add_argument(
        "syntelog_input_file", type=str, help="parent path of syntelog file"
    )
    parser.add_argument(
        "--output_directory",
        "-o",
        type=str,
        help="parent path of output directory",
        default=os.path.join(path_main, "../../../../examples/Rice_Synteny/results"),
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.syntelog_input_file = os.path.abspath(args.syntelog_input_file)
    args.output_directory = os.path.abspath(args.output_directory)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    # Process
    logger.info("Starting filtration...")
    process(
        args.syntelog_input_file, args.output_directory,
    )
