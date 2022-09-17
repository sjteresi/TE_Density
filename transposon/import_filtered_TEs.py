import pandas as pd
from transposon import check_nulls


def import_filtered_TEs(tes_input_path, logger):
    """
    Import a pre-filtered TE file to the pipeline

    Args:
        tes_input_path (str): String to path of pre-filtered TE file

        logger (logging.Logger): object to log messages

    Returns:
        transposon_data (pandas.core.frame.DataFrame): A pandas dataframe
        representing the preprocessed transposon annotation file
    """
    try:
        transposon_data = pd.read_csv(
            tes_input_path,
            header="infer",
            sep="\t",
            dtype={
                "Start": "float64",
                "Stop": "float64",
                "Length": "float64",
                "Chromosome": str,
                "Strand": str,
                "Order": str,
                "SuperFamily": str,
            },
        )
    except Exception as err:
        msg = ("Error occurred while trying to read preprocessed TE "
               "annotation file into a Pandas dataframe, please refer "
               "to the README as to what information is expected")
        logger.critical(msg)
        raise err

    check_nulls(transposon_data, logger)

    # Sort for legibility
    transposon_data.sort_values(by=["Chromosome", "Start"], inplace=True)

    logger.info("import of pre-filtered transposon annotation... success!")

    return transposon_data
