#!/usr/bin/env python3

"""
Contains information relevant to the specific strawberry genomes.
"""

__author__ = "Scott Teresi"

import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class GenomeData(object):
    """Wraps two dataframes, the transposable element dataframe and the gene
    dataframe, both of which are derived from annotation files.

    Used for generating descriptive statistics and graphs

    """

    def __init__(
        self, genome_id, gene_dataframe, transposon_dataframe, genome_size=None
    ):
        """Initialie the genome information, which is an amalgamation of gene
        and transposon data.
        Args:
            genome_id (str): String of the genome id.
            gene_dataframe (pandas.DataFrame): gene dataframe
            transposon_dataframe (pandas.DataFrame): transposon dataframe
            genome_size (float): give the approximate genome size in MB
        """

        self.genome_id = genome_id
        self.gene_dataframe = gene_dataframe
        self.transposon_dataframe = transposon_dataframe
        self.genome_size = genome_size  # genome size in Mb

        self.g_starts = self.gene_dataframe.Start
        self.t_starts = self.transposon_dataframe.Start
        self.g_stops = self.gene_dataframe
        self.t_stops = self.transposon_dataframe

        self.total_G_lengths = self.gene_dataframe.Length.sum()
        self.total_T_lengths = self.transposon_dataframe.Length.sum()

        self.total_G_lengths_MB = self.gene_dataframe.Length.sum() / 1000000000
        self.total_T_lengths_MB = self.transposon_dataframe.Length.sum() / 1000000000

        if self.genome_size != None:
            self.total_other_lengths_MB = (
                self.genome_size - self.total_G_lengths_MB - self.total_T_lengths_MB
            )

        self.num_of_TEs = len(self.transposon_dataframe.index)
        self.num_of_genes = len(self.gene_dataframe.index)

        self.subgenomes = None

    @property
    def Chromosomes(self):
        """
        Returns a list of chromosomes, also checks to see if the chromosomes
        between the gene_dataframe and the transposon_dataframe are equal.

        """
        comparison = (
            self.gene_dataframe.Chromosome.unique()
            == self.transposon_dataframe.Chromosome.unique()
        )

        if comparison.all():
            return self.gene_dataframe.Chromosome.unique().tolist()
        else:
            raise ValueError(
                """The chromosomes in the transposon dataframe do
                             not match the chromosomes in the gene
                             dataframe."""
            )

    def subgenome_chromosome(self):
        pass

    def split(self, df, group):
        gb = df.groupby(group)
        return [gb.get_group(x) for x in gb.groups]

    def order_transposon_subset(self, my_order):
        # x =  self.transposon_dataframe[self.transposon_dataframe['Order']==my_order]
        # return x
        # return self.transposon_dataframe.loc[self.transposon_dataframe['Order']==my_order]
        return self.transposon_dataframe.loc[
            self.transposon_dataframe["Order"] == my_order
        ].copy(deep=True)

    @staticmethod
    def _concat_subgenomes(grouped_genes, grouped_TEs, genome_chrom_list, selection):
        if selection == "TEs":
            return pd.concat(
                [
                    sub_TE
                    for sub_gene, sub_TE in zip(grouped_genes, grouped_TEs)
                    if sub_gene.Chromosome.unique()
                    and sub_TE.Chromosome.unique() in genome_chrom_list
                ]
            )

        if selection == "Genes":
            return pd.concat(
                [
                    sub_gene
                    for sub_gene, sub_TE in zip(grouped_genes, grouped_TEs)
                    if sub_gene.Chromosome.unique()
                    and sub_TE.Chromosome.unique() in genome_chrom_list
                ]
            )

    def Cam_subgenomes(self):
        """
        Returns a list of GenomeData objects that are the Camarosa subgenomes.
        """
        F_vesca = ["Fvb1-4", "Fvb2-2", "Fvb3-4", "Fvb4-3", "Fvb5-1", "Fvb6-1", "Fvb7-2"]
        F_nipponica = [
            "Fvb1-3",
            "Fvb2-1",
            "Fvb3-3",
            "Fvb4-2",
            "Fvb5-4",
            "Fvb6-2",
            "Fvb7-1",
        ]
        F_iinumae = [
            "Fvb1-2",
            "Fvb2-4",
            "Fvb3-2",
            "Fvb4-4",
            "Fvb5-4",
            "Fvb6-3",
            "Fvb7-3",
        ]
        F_viridis = [
            "Fvb1-1",
            "Fvb2-3",
            "Fvb3-1",
            "Fvb4-1",
            "Fvb5-2",
            "Fvb6-4",
            "Fvb7-4",
        ]
        grouped_genes = self._split(self.gene_dataframe, "Chromosome")
        grouped_TEs = self._split(self.transposon_dataframe, "Chromosome")

        subgenome_dict = {}

        F_viridis_Genes = self._concat_subgenomes(
            grouped_genes, grouped_TEs, F_viridis, "Genes"
        )
        F_viridis_TEs = self._concat_subgenomes(
            grouped_genes, grouped_TEs, F_viridis, "TEs"
        )

        F_nipponica_Genes = self._concat_subgenomes(
            grouped_genes, grouped_TEs, F_nipponica, "Genes"
        )
        F_nipponica_TEs = self._concat_subgenomes(
            grouped_genes, grouped_TEs, F_nipponica, "TEs"
        )

        F_iinumae_Genes = self._concat_subgenomes(
            grouped_genes, grouped_TEs, F_iinumae, "Genes"
        )
        F_iinumae_TEs = self._concat_subgenomes(
            grouped_genes, grouped_TEs, F_iinumae, "TEs"
        )

        F_vesca_Genes = self._concat_subgenomes(
            grouped_genes, grouped_TEs, F_vesca, "Genes"
        )
        F_vesca_TEs = self._concat_subgenomes(
            grouped_genes, grouped_TEs, F_vesca, "TEs"
        )

        subgenome_dict["F_Viridis"] = GenomeData(
            "F_viridis", F_viridis_Genes, F_viridis_TEs
        )

        subgenome_dict["F_Nipponica"] = GenomeData(
            "F_nipponica", F_nipponica_Genes, F_nipponica_TEs
        )

        subgenome_dict["F_Iinumae"] = GenomeData(
            "F_iinumae", F_iinumae_Genes, F_iinumae_TEs
        )

        subgenome_dict["F_Vesca"] = GenomeData("F_Vesca", F_vesca_Genes, F_vesca_TEs)

        return subgenome_dict

    @property
    def whole_genome_percent_gene_lengths(self):
        """

        """
        return self.total_G_lengths / 1000000000 / self.genome_size

    @property
    def whole_genome_percent_transposon_lengths(self):
        """

        """
        return self.total_T_lengths / 1000000000 / self.genome_size

    @property
    def whole_genome_percent_other_lengths(self):
        """

        """
        return (
            1
            - whole_genome_percent_gene_lengths
            - whole_genome_percent_transposon_lengths
        )

    def number_of_elements_per_grouping(self, grouping):
        """
        Grouping can be Order or SuperFamily

        Maybe make a check to see if people pass that in right


        Returns a dictionary of counts per grouping
        """
        my_grouping_counts = getattr(self.transposon_dataframe, grouping).value_counts()
        return my_grouping_counts.to_dict()

    @property
    def transposon_Order_number_dictionary(self):
        return self.number_of_elements_per_grouping("Order")

    @property
    def transposon_SuperFam_number_dictionary(self):
        return self.number_of_elements_per_grouping("SuperFamily")

    @property
    def order_sum_sequence_len_dictionary(self):
        return self.transposon_dataframe.groupby(["Order"]).Length.sum().to_dict()

    @property
    def order_sum_seq_len_dict_MB(self):
        return (
            self.transposon_dataframe.groupby(["Order"])
            .Length.sum()
            .divide(1000000000)
            .to_dict()
        )

    @property
    def superfam_sum_sequence_len_dictionary(self):
        return self.transposon_dataframe.groupby(["SuperFamily"]).Length.sum().to_dict()

    @property
    def superfam_sum_seq_len_dict_MB(self):
        return (
            self.transposon_dataframe.groupby(["SuperFamily"])
            .Length.sum()
            .divide(1000000000)
            .to_dict()
        )

    @property
    def average_order_length(self):
        return self.transposon_dataframe.groupby(["Order"]).Length.mean()

    @property
    def average_superfam_length(self):
        return self.transposon_dataframe.groupby(["SuperFamily"]).Length.mean()

    @property
    def median_order_length(self):
        return self.transposon_dataframe.groupby(["Order"]).Length.median()

    @property
    def median_superfam_length(self):
        return self.transposon_dataframe.groupby(["SuperFamily"]).Length.median()

    @property
    def orders_as_percent_sequences(self):
        """
        """
        my_dict = {}
        for key, val in self.order_sum_sequence_len_dictionary.items():
            my_dict[key] = val / 1000000000 / self.genome_size
        return my_dict

    def __repr__(self):
        """Printable representation."""

        info = """
               Genome ID: {self.genome_id}
               Genome Size: {self.genome_size}
               Genome Chromosomes: {self.Chromosomes}
               """
        return info.format(self=self)


class SubgenomeData(GenomeData):
    """
    Supply a regular genedataframe with a regular transposon dataframe, with a
    list of the chromosomes belonging to a given subgenome, and return a
    modified GenomeData object
    """

    def __init__(
        self,
        genome_id,
        gene_dataframe,
        transposon_dataframe,
        chromosome_grouping,
        subgenome_identity,
        genome_size=None,
    ):
        """
        Args:


            chromosome_grouping (list of strings): A user supplied list of strings of
                the chromosomes for that given subgenome_identity.
            subgenome_identity (string): The name of the subgenome
        """
        self.gene_dataframe = gene_dataframe.loc[
            gene_dataframe["Chromosome"].isin(chromosome_grouping)
        ]
        self.transposon_dataframe = transposon_dataframe.loc[
            transposon_dataframe["Chromosome"].isin(chromosome_grouping)
        ]
        self.genome_id = genome_id
        self.chromosome_grouping = chromosome_grouping
        self.subgenome_identity = subgenome_identity
        self.genome_size = self.get_subgenome_genome_size(
            gene_dataframe, transposon_dataframe, chromosome_grouping
        )
        super().__init__(
            genome_id, self.gene_dataframe, self.transposon_dataframe, self.genome_size
        )

    @staticmethod
    def _split(df, group):
        gb = df.groupby(group)
        return [gb.get_group(x) for x in gb.groups]

    def get_subgenome_genome_size(
        self, gene_dataframe, transposon_dataframe, chromosome_grouping
    ):

        chrom_g_groups = self._split(gene_dataframe, "Chromosome")
        chrom_T_groups = self._split(transposon_dataframe, "Chromosome")

        counter = 0
        for g_chrom, T_chrom in zip(chrom_g_groups, chrom_T_groups):
            if (
                g_chrom.Chromosome.unique()
                and T_chrom.Chromosome.unique() in chromosome_grouping
            ):
                max_stop = max(g_chrom.Stop.max(), T_chrom.Stop.max())
                counter += max_stop
        self.genome_size = counter / 1000000000
        return self.genome_size

    def __repr__(self):
        """
        Printable object representation
        """
        info = """
               Subgenome Identity: {self.subgenome_identity}
               Chromosome Grouping: {self.chromosome_grouping}
               Gene Dataframe Head: {self.gene_dataframe}
               TE Dataframe Head: {self.transposon_dataframe}
               """
        return info.format(self=self)
