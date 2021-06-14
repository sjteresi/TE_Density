#!/usr/bin/env python3

__author__ = "Scott Teresi"

import logging
import pandas as pd


def import_syntelogs(syntelog_input_file):
    """
    Import the syntelogs from the raw file and manage data filtration
    """

    col_names = [
        "OrgA_Chromosome",
        "OrgA_Gene_Region",
        "OrgA_Start",
        "OrgA_Stop",
        "OrgB_Chromosome",
        "OrgB_Gene_Region",
        "OrgB_Start",
        "OrgB_Stop",
        "E_Value",
        "Diagonal_Score",
        "Web_Link",
    ]

    col_to_use = [
        "OrgA_Chromosome",
        "OrgA_Gene_Region",
        "OrgB_Chromosome",
        "OrgB_Gene_Region",
        "E_Value",
        "Diagonal_Score",
    ]

    gene_data = pd.read_csv(
        syntelog_input_file,
        sep="\t+",
        header=None,
        engine="python",
        names=col_names,
        usecols=col_to_use,
        comment="#",
        dtype={
            "OrgA_Chromosome": str,
            "OrgA_Gene_Region": str,
            "OrgB_Chromosome": str,
            "OrgB_Gene_Region": str,
            "E_Value": "float64",
            "Diagonal_Score": "int32",
        },
    )

    # Get the correct name for the genes
    # MAGIC to split the name correctly
    gene_data["OrgA_Gene_Region"] = (
        gene_data["OrgA_Gene_Region"].str.split("\|\|").str[3]
    )
    gene_data["OrgB_Gene_Region"] = (
        gene_data["OrgB_Gene_Region"].str.split("\|\|").str[3]
    )

    # Remove rows that have transcript in the name because not worth dealing
    # with for example
    gene_data = gene_data[~gene_data["OrgA_Gene_Region"].str.contains("transcript")]
    gene_data = gene_data[~gene_data["OrgB_Gene_Region"].str.contains("transcript")]

    # Get the correct name for the gene names
    # MAGIC to split the name correctly
    gene_data["OrgA_Gene_Region"] = (
        gene_data["OrgA_Gene_Region"].str.split("CDS:").str[1]
    )
    gene_data["OrgB_Gene_Region"] = (
        gene_data["OrgB_Gene_Region"].str.split("CDS:").str[1]
    )

    gene_data["OrgA_Gene_Region"] = gene_data["OrgA_Gene_Region"].str.split(".").str[0]
    gene_data["OrgB_Gene_Region"] = gene_data["OrgB_Gene_Region"].str.split("-").str[0]

    # SynMap returns the transcript name for Sativa which can have slight
    # differences with the gene name, namely the letter g is replaced with
    # the letter t.
    # NB fix to get the gene name
    gene_data["OrgB_Gene_Region"] = gene_data["OrgB_Gene_Region"].str.replace("t", "g")

    # Get the correct name for the chromosome
    # MAGIC
    gene_data["OrgA_Chromosome"] = gene_data["OrgA_Chromosome"].str.split("_").str[1]
    gene_data["OrgB_Chromosome"] = gene_data["OrgB_Chromosome"].str.split("_").str[1]

    # This step is important, it could differ if your data input is different.
    gene_data.rename(
        columns={"OrgA_Gene_Region": "Glaberrima", "OrgB_Gene_Region": "Sativa"},
        inplace=True,
    )
    # Trim E-values less than 0.05
    # MAGIC
    gene_data = gene_data.loc[gene_data["E_Value"] < 0.05]

    gene_data.drop(
        columns=["Diagonal_Score"],
        inplace=True,
    )

    # I only want pairs where the chromosomes are equal
    gene_data = gene_data.loc[
        gene_data["OrgA_Chromosome"] == gene_data["OrgB_Chromosome"]
    ]

    chromosome_list = [str(i) for i in range(1, 12 + 1)]
    gene_data = gene_data.loc[gene_data["OrgA_Chromosome"].isin(chromosome_list)]

    return gene_data


class Syntelog_Data(object):
    """
    Wrappers for input data, multiple syntelog pairs.

    Used to provide a common interface and fast calculations with numpy.
    """

    def __init__(self, syntelog_dataframe, logger=None):
        """Initialize.

        Args:
            syntelog_dataframe (DataFrame): syntelog data frame.
        """
        self._logger = logger or logging.getLogger(__name__)
        self.dataframe = syntelog_dataframe

    def save_to_disk(self, filename):
        """
        Save the syntelogs to disk in a 2-column format.
        Arabidopsis in left-hand column, blueberry in right-hand column.

        Args:
            filename (str): path for the file
        """
        self.dataframe.to_csv(filename, sep="\t", header=True, index=False)

    def __repr__(self):
        """
        String representation for developer.
        """
        return "Syntelog_Data{}".format(self.dataframe)
