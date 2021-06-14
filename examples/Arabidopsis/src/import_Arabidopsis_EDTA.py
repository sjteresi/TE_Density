"""
Filter a EDTA-created TE annotation to the appropriate format for TE
Density algorithm
"""

__author__ = "Scott Teresi"

import pandas as pd
import argparse
import os
import logging
import coloredlogs

from examples.Arabidopsis.src.replace_names_Arabidopsis import te_annot_renamer


def check_nulls(my_df, logger):
    """Check the TE dataframe for ANY null values in ANY rows

    Args:
        my_df (pandas.core.DataFrame): Pandas dataframe of TE values from TE
            annotation
    """
    Bool = my_df.isnull().values.any()
    if Bool:
        logger.critical("You have null values in your dataframe!")
        logger.critical("Here are the null values in the output:")
        null_columns = my_df.columns[my_df.isnull().any()]
        print((my_df[my_df.isnull().any(axis=1)][null_columns].head()))


def write_cleaned_transposons(te_pandaframe, output_dir, old_filename, logger):
    file_name = os.path.join(
        output_dir,
        ("Cleaned_" + os.path.splitext(os.path.basename(old_filename))[0]) + ".tsv",
    )  # MAGIC to get proper extension

    logger.info("Writing cleaned TE file to: %s" % file_name)
    te_pandaframe.to_csv(file_name, sep="\t", header=True, index=False)


def import_transposons(tes_input_path, te_annot_renamer, logger):
    """Import TE file and read as a dataframe in Pandas

    Args:
        tes_input_path (str): string of the file path to the TE annotation

        logger (logging obj): The object to call logs and info
    """
    col_names = [
        "Chromosome",
        "Software",
        "Feature",
        "Start",
        "Stop",
        "Score",
        "Strand",
        "Phase",
        "Attribute",
    ]

    te_data = pd.read_csv(
        tes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        comment="#",
        dtype={"Start": "float64", "Stop": "float64"},
    )

    # Drop extraneous columns
    te_data.drop(columns=["Score", "Software", "Phase", "Feature"], inplace=True)

    # Create Order and SuperFamily column from Attribute column
    # Because that column contains the detailed TE information
    # Then remove old Attribute column
    te_data["Attribute"] = te_data["Attribute"].str.extract(r"Classification=(.*?);")
    te_data[["Order", "SuperFamily"]] = te_data.Attribute.str.split("/", expand=True)
    te_data.drop(columns=["Attribute"], inplace=True)
    te_data.Order = te_data.Order.astype(str)
    te_data.SuperFamily = te_data.SuperFamily.astype(str)
    te_data.Strand = te_data.Strand.astype(str)

    # Call renamer
    te_data = te_annot_renamer(te_data)

    # Declare data types
    te_data["Length"] = te_data.Stop - te_data.Start + 1
    check_nulls(te_data, logger)

    te_data.sort_values(by=["Chromosome", "Start"], inplace=True)

    return te_data


if __name__ == "__main__":

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

    # Execute
    cleaned_transposons = import_transposons(
        args.TE_input_file, te_annot_renamer, logger
    )
    write_cleaned_transposons(
        cleaned_transposons, args.output_dir, args.TE_input_file, logger
    )
