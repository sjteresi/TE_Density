import pandas as pd
from transposon import check_nulls, check_strand


def import_filtered_genes(genes_input_path, logger):
    """
    Import the preprocessed gene annotation file. Read it as pandas dataframe

    genes_input_path (str): Path to cleaned input annotation file of genes

    logger (logging.Logger): Logging object

    Returns:
        gene_data (pandas.core.frame.DataFrame): A pandas dataframe
        representing the preprocessed gene annotation file
    """
    try:
        gene_data = pd.read_csv(
            genes_input_path,
            header="infer",
            sep="\t",
            dtype={
                "Start": "float64",
                "Stop": "float64",
                "Length": "float64",
                "Chromosome": str,
                "Strand": str,
                "Feature": str,
                "Gene_Name": str,
            },
        )
    except Exception as err:
        msg = """
            Error occurred while trying to read preprocessed gene
            annotation file into a Pandas dataframe, please refer
            to the README as to what information is expected
            input file: %s
            """
        logger.critical(msg, genes_input_path)
        raise err

    gene_data.set_index("Gene_Name", verify_integrity=True, inplace=True)
    check_nulls(gene_data, logger)
    check_strand(gene_data, logger)

    # NOTE edit Strand '.' values to be sense orientation as
    # described in check_strand()
    gene_data["Strand"].replace(to_replace={".": "+"}, inplace=True)

    # Sort for legibility
    gene_data.sort_values(by=["Chromosome", "Start"], inplace=True)

    logger.info("import of preprocessed gene annotation... success!")
    return gene_data
