#!/usr/bin/env python3

__author__ = "Scott Teresi"

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

    syntelog_pandaframe = pd.read_csv(
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
    syntelog_pandaframe["OrgA_Gene_Region"] = (
        syntelog_pandaframe["OrgA_Gene_Region"].str.split("\|\|").str[3]
    )
    syntelog_pandaframe["OrgB_Gene_Region"] = (
        syntelog_pandaframe["OrgB_Gene_Region"].str.split("\|\|").str[3]
    )

    # Remove rows that have transcript in the name because not worth dealing
    # with for example
    syntelog_pandaframe = syntelog_pandaframe[
        ~syntelog_pandaframe["OrgA_Gene_Region"].str.contains("transcript")
    ]
    syntelog_pandaframe = syntelog_pandaframe[
        ~syntelog_pandaframe["OrgB_Gene_Region"].str.contains("transcript")
    ]

    # Get the correct name for the gene names
    # MAGIC to split the name correctly
    syntelog_pandaframe["OrgA_Gene_Region"] = (
        syntelog_pandaframe["OrgA_Gene_Region"].str.split("CDS:").str[1]
    )
    syntelog_pandaframe["OrgB_Gene_Region"] = (
        syntelog_pandaframe["OrgB_Gene_Region"].str.split("CDS:").str[1]
    )

    syntelog_pandaframe["OrgA_Gene_Region"] = (
        syntelog_pandaframe["OrgA_Gene_Region"].str.split(".").str[0]
    )
    syntelog_pandaframe["OrgB_Gene_Region"] = (
        syntelog_pandaframe["OrgB_Gene_Region"].str.split("-").str[0]
    )

    # SynMap returns the transcript name for Sativa which can have slight
    # differences with the gene name, namely the letter g is replaced with
    # the letter t.
    # NB fix to get the gene name
    syntelog_pandaframe["OrgB_Gene_Region"] = syntelog_pandaframe[
        "OrgB_Gene_Region"
    ].str.replace("t", "g")

    # Get the correct name for the chromosome
    # MAGIC
    syntelog_pandaframe["OrgA_Chromosome"] = (
        syntelog_pandaframe["OrgA_Chromosome"].str.split("_").str[1]
    )
    syntelog_pandaframe["OrgB_Chromosome"] = (
        syntelog_pandaframe["OrgB_Chromosome"].str.split("_").str[1]
    )

    # This step is important, it could differ if your data input is different.
    syntelog_pandaframe.rename(
        columns={"OrgA_Gene_Region": "Glaberrima", "OrgB_Gene_Region": "Sativa"},
        inplace=True,
    )
    # Trim E-values less than 0.05
    # MAGIC
    syntelog_pandaframe = syntelog_pandaframe.loc[syntelog_pandaframe["E_Value"] < 0.05]

    syntelog_pandaframe.drop(
        columns=["Diagonal_Score"],
        inplace=True,
    )

    # I only want pairs where the chromosomes are equal
    syntelog_pandaframe = syntelog_pandaframe.loc[
        syntelog_pandaframe["OrgA_Chromosome"] == syntelog_pandaframe["OrgB_Chromosome"]
    ]

    chromosome_list = [str(i) for i in range(1, 12 + 1)]
    syntelog_pandaframe = syntelog_pandaframe.loc[
        syntelog_pandaframe["OrgA_Chromosome"].isin(chromosome_list)
    ]

    return syntelog_pandaframe
