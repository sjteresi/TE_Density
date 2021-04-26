import pandas as pd
from transposon import check_nulls


def import_filtered_genes(genes_input_path, logger):
    """
    TODO fill
    """
    gene_data = pd.read_csv(
        genes_input_path,
        header="infer",
        engine="python",
        sep="\t",
        index_col="Gene_Name",
        dtype={
            "Start": "float64",
            "Stop": "float64",
            "Length": "float64",
            "Chromosome": str,
            "Strand": str,
            "Feature": str,
        },
    )

    check_nulls(gene_data, logger)
    # Sort for legibility
    gene_data.sort_values(by=["Chromosome", "Start"], inplace=True)

    return gene_data
