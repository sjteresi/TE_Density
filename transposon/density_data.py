#!/usr/bin/env/python

"""
Read TE Density H5 files
"""

__author__ = "Scott Teresi"

import h5py
import numpy as np
import os
import shutil


class DensityData:
    def __init__(self, input_h5, gene_data, logger, sense_swap=True):
        """
        input_h5 (str): Path to h5 file of TE Density output.
        gene_data (GeneData):
        """
        new_filename = input_h5.replace(".h5", "_SenseSwapped.HDF5")
        if not os.path.isfile(new_filename):
            shutil.copyfile(input_h5, new_filename)  # copy because we
            # need to swap values
        self.data_frame = h5py.File(new_filename, "r+")

        self.key_list = list(self.data_frame.keys())
        self.gene_list = self.data_frame["GENE_NAMES"][:]
        self.num_genes = len(self.gene_list)
        self.chromosomes = self.data_frame["CHROMOSOME_ID"][:]
        self.unique_chromosomes = list(set(self.chromosomes))
        if len(self.unique_chromosomes) != 1:
            raise ValueError(
                "There are multiple unique chromosomes in this density data."
            )
        self.unique_chromosome_id = self.unique_chromosomes[0]  # MAGIC
        # self.name = self.chromosomes.unique()
        self.windows = self.data_frame["WINDOWS"]  # not int
        self.window_list = [int(i) for i in self.windows[:]]  # list of ints

        self.order_list = list(self.data_frame["ORDER_NAMES"][:])
        self.super_list = list(self.data_frame["SUPERFAMILY_NAMES"][:])

        # NB. Shape for these is (type of TE, window, gene)
        self.left_orders = self.data_frame["RHO_ORDERS_LEFT"]
        self.intra_orders = self.data_frame["RHO_ORDERS_INTRA"]
        self.right_orders = self.data_frame["RHO_ORDERS_RIGHT"]

        self.left_supers = self.data_frame["RHO_SUPERFAMILIES_LEFT"]
        self.intra_supers = self.data_frame["RHO_SUPERFAMILIES_INTRA"]
        self.right_supers = self.data_frame["RHO_SUPERFAMILIES_RIGHT"]

        self.genome_id = gene_data.genome_id

        gene_data.data_frame.Strand.replace({"+": 1, "-": 0}, inplace=True)
        strands_as_numpy = gene_data.data_frame.Strand.to_numpy(copy=False)
        zero_gene_list = np.where(strands_as_numpy == 0)[0]  # gets indices of 0s
        # NB the 0s are the antisense
        genes_to_swap = gene_data.data_frame.iloc[zero_gene_list, :].index.tolist()

        if sense_swap:
            self._swap_strand_vals(genes_to_swap)

    def _index_of_gene(self, gene_string):
        """Return the index of a gene given the name of the gene
        Args:
            gene_string (str): string representing the name of the gene
        Returns:
            Returns an index of the gene in the H5 dataset
        """
        if gene_string not in self.gene_list:
            raise IndexError(
                """The gene '%s' is not in the density data,
                please verify that the list of genes to swap sense and
                antisense values does not have any genes that are not present
                in the density data (h5 file). The gene list is: %s"""
                % (gene_string, self.gene_list)
            )
        return np.where(self.gene_list == gene_string)[0][0]  # MAGIC

    def nonzero_indices(self):
        pass
        # print(np.nonzero(data.left_orders[:]))

    @property
    def order_index_dict(self):
        """Returns a dictionary of TE order names as keys and indices as values"""
        order_dict = {}
        for i in range(len(self.order_list)):
            order_dict[self.order_list[i]] = i
        # TODO remove the revision groupings only for the dotplots
        # Removed temporarily for other coding
        # order_dict.pop("S_Revision")  # MAGIC pop the revision set
        return order_dict

    @property
    def super_index_dict(self):
        """Returns a dictionary of TE superfamily names as keys and indices as
        values"""
        super_dict = {}
        for i in range(len(self.super_list)):
            super_dict[self.super_list[i]] = i
        # TODO remove the revision groupings only for the dotplots
        # Removed temporarily for other coding
        # super_dict.pop("O_Revision")  # MAGIC pop the revision set
        return super_dict

    def _swap_strand_vals(self, gene_names):
        """Switch density values for the genes in which it is antisense due to
        the fact that antisense genes point in the opposite direction to sense
        genes

        Args:
            gene_names(list of str):
        """
        for name in gene_names:
            index_to_switch = self._index_of_gene(name)

            # SUPERFAMILIES
            left_val_super = self.data_frame["RHO_SUPERFAMILIES_LEFT"][
                :, :, index_to_switch
            ]
            right_val_super = self.data_frame["RHO_SUPERFAMILIES_RIGHT"][
                :, :, index_to_switch
            ]
            # ORDERS
            left_val_order = self.data_frame["RHO_ORDERS_LEFT"][:, :, index_to_switch]

            right_val_order = self.data_frame["RHO_ORDERS_RIGHT"][:, :, index_to_switch]

            # Reassign
            self.data_frame["RHO_SUPERFAMILIES_RIGHT"][
                :, :, index_to_switch
            ] = left_val_super

            self.data_frame["RHO_SUPERFAMILIES_LEFT"][
                :, :, index_to_switch
            ] = right_val_super
            self.data_frame["RHO_ORDERS_RIGHT"][:, :, index_to_switch] = left_val_order

            self.data_frame["RHO_ORDERS_LEFT"][:, :, index_to_switch] = right_val_order

    def __repr__(self):
        """String representation for developer."""
        info = """
                DensityData Genome ID: {self.genome_id}
                DensityData shape key: (type of TE, windows, genes)
                DensityData left_order shape: {self.left_orders.shape}
                DensityData intra_order shape: {self.intra_orders.shape}
                DensityData right_order shape: {self.right_orders.shape}
                DensityData left_supers shape: {self.left_supers.shape}
                DensityData intra_supers shape: {self.intra_supers.shape}
                DensityData right_supers shape: {self.right_supers.shape}
            """
        return info.format(self=self)
