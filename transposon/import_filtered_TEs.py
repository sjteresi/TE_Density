import pandas as pd
from transposon import check_nulls

def import_filtered_TEs(tes_input_path, logger):
    """
    Import a pre-filtered TE file to the pipeline

    Args:
        tes_input_path (str): String to path of pre-filtered TE file

        logger (logging.Logger): object to log messages
    """
    transposon_data = pd.read_csv(
        tes_input_path,
        header="infer",
        engine='python',
        sep="\t",
        dtype={"Start": "float64", "Stop": "float64", "Length": "float64",
               'Chromosome': str, 'Strand': str, 'Order': str, 'SuperFamily':
               str },
    )
    check_nulls(transposon_data, logger)
    # Sort for legibility
    transposon_data.sort_values(by=["Chromosome", "Start"], inplace=True)

    return transposon_data
