import pandas as pd


def check_nulls(my_df, logger):
    """Check the TE dataframe for ANY null values in ANY rows

    Args:
        my_df (pandas.core.DataFrame): Pandas dataframe of TE values from TE
            annotation
    """
    Bool = my_df.isnull().values.any()
    if Bool:
        logger.critical("You have null values in your dataframe!")


def import_transposons(tes_input_path, te_annot_renamer, contig_del, logger):
    """Import TE file and read as a dataframe in Pandas
        Args:
            tes_input_path (str): string of the file path to the TE annotation

            te_annot_renamer (function containing a dictionary and other methods):
                imported from separate file within the repository. This file
                performs the more specific filtering steps on the TEs such as
                changing the annotation details for specific TE types.

            contig_drop (bool): logical whether to drop rows with a contig as the
                chromosome id

            logger (logging obj)
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
        "Attribute",
    ]

    col_to_use = ["Chromosome", "Software", "Feature", "Start", "Stop", "Strand"]

    TE_Data = pd.read_csv(
        tes_input_path,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        dtypes={'Start': 'float64', 'Stop': 'float64', 'Chromosome': str,
                'Strand': str},
        usecols=col_to_use,
    )

    TE_Data = TE_Data[~TE_Data.Chromosome.str.contains("#")]  # remove comment
    # rows in annotation
    TE_Data[["Order", "SuperFamily"]] = TE_Data.Feature.str.split("/", expand=True)

    TE_Data = TE_Data.drop(["Feature", "Software"], axis=1)
    TE_Data = te_annot_renamer(TE_Data)  # NOTE call to the cleaner
    TE_Data["Length"] = TE_Data.Stop - TE_Data.Start + 1
    if contig_del:
        TE_Data = TE_Data[~TE_Data.Chromosome.str.contains("contig", case=False)]
    check_nulls(TE_Data, logger)
    return TE_Data
