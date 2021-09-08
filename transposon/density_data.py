#!/usr/bin/env/python

"""
Read TE Density H5 files
"""

__author__ = "Scott Teresi"

import h5py
import numpy as np
import os
import shutil
import re
from collections import namedtuple

DensitySlice = namedtuple(
    "DensitySlice", ["slice", "direction", "window_val", "te_type"]
)


class DensityData:
    def __init__(self, input_h5, gene_data, logger, sense_swap=True):
        """
        input_h5 (str): Path to h5 file of TE Density output.
        gene_data (GeneData):
        """
        new_filename = input_h5.replace(".h5", "_SenseSwapped.HDF5")
        # if not os.path.isfile(new_filename):
        shutil.copyfile(input_h5, new_filename)  # copy because we
        # need to swap values later
        self.data_frame = h5py.File(new_filename, "r+")

        # Remove the revision groupings from the dataset as they are an
        # artifact from accurately merging TEs
        self.data_frame = self.remove_revision_sets()

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

        # TODO
        # Do i need to call self.data_frame.close() here?

    def remove_revision_sets(self):
        for te_order_idx, te_order_name in enumerate(self.data_frame["ORDER_NAMES"][:]):
            if te_order_name == "S_Revision":
                order_idx_to_del = te_order_idx
                order_names_less_revision = np.delete(
                    self.data_frame["ORDER_NAMES"], [order_idx_to_del]
                )
                del self.data_frame["ORDER_NAMES"]
                self.data_frame.create_dataset(
                    "ORDER_NAMES", data=order_names_less_revision
                )

        for te_super_idx, te_super_name in enumerate(
            self.data_frame["SUPERFAMILY_NAMES"][:]
        ):
            if te_super_name == "O_Revision":
                super_idx_to_del = te_super_idx
                super_names_less_revision = np.delete(
                    self.data_frame["SUPERFAMILY_NAMES"], [super_idx_to_del]
                )
                del self.data_frame["SUPERFAMILY_NAMES"]
                self.data_frame.create_dataset(
                    "SUPERFAMILY_NAMES", data=super_names_less_revision
                )

        #
        order_datasets = ["RHO_ORDERS_LEFT", "RHO_ORDERS_INTRA", "RHO_ORDERS_RIGHT"]
        super_datasets = [
            "RHO_SUPERFAMILIES_LEFT",
            "RHO_SUPERFAMILIES_INTRA",
            "RHO_SUPERFAMILIES_RIGHT",
        ]

        for order_dataset in order_datasets:
            data_less_revision = np.delete(
                self.data_frame[order_dataset], order_idx_to_del, 0
            )
            del self.data_frame[order_dataset]
            self.data_frame.create_dataset(order_dataset, data=data_less_revision)

        for super_dataset in super_datasets:
            data_less_revision = np.delete(
                self.data_frame[super_dataset], super_idx_to_del, 0
            )
            del self.data_frame[super_dataset]
            self.data_frame.create_dataset(super_dataset, data=data_less_revision)

        return self.data_frame

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
        return order_dict

    @property
    def super_index_dict(self):
        """Returns a dictionary of TE superfamily names as keys and indices as
        values"""
        super_dict = {}
        for i in range(len(self.super_list)):
            super_dict[self.super_list[i]] = i
        return super_dict

    def yield_all_slices(self):
        """
        Yields a DensitySlice (named tuple) object for each TE
        type/direction/window combination for all genes in the DensityData obj.
        """
        directions = ["Upstream", "Downstream"]
        for direction in directions:
            for window_idx, window_val in enumerate(self.window_list):
                for te_type, te_order_idx in self.order_index_dict.items():
                    if direction == "Upstream":
                        yield DensitySlice(
                            self.left_orders[te_order_idx, window_idx, :],
                            direction,
                            window_val,
                            te_type,
                        )
                    if direction == "Downstream":
                        yield DensitySlice(
                            self.right_orders[te_order_idx, window_idx, :],
                            direction,
                            window_val,
                            te_type,
                        )

                # Iterate over SuperFamilies now
                for te_type, te_super_idx in self.super_index_dict.items():
                    if direction == "Upstream":
                        yield DensitySlice(
                            self.left_supers[te_super_idx, window_idx, :],
                            direction,
                            window_val,
                            te_type,
                        )
                    if direction == "Downstream":
                        yield DensitySlice(
                            self.right_supers[te_super_idx, window_idx, :],
                            direction,
                            window_val,
                            te_type,
                        )

    @classmethod
    def verify_h5_cache(cls, h5_file, gene_data_instance, logger):
        """
        Constructor
        Verify that previously value swapped H5 values are saved to disk and import
        those instead of swapping the values once more

        Args:
            h5_file (str): Path to string of TE density HDF5 file saved to disk
            gene_data_instance (GeneData): An instance of GeneData
            logger (logging.Logger):

        Returns:
            processed_density_data (DensityData): An instance of DensityData that
            combines information from the TE density HDF5 data and GeneData
        """
        if os.path.exists(h5_file.replace(".h5", "_SenseSwapped.HDF5")):
            logger.info("Previous sense swapped data exists, reading...")
            processed_density_data = cls(
                h5_file, gene_data_instance, logger, sense_swap=False
            )
        else:
            logger.info("Writing new sense swapped DensityData...")
            processed_density_data = cls(h5_file, gene_data_instance, logger)
        return processed_density_data

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

    @staticmethod
    def supply_density_data_files(path_to_folder):
        """
        Iterate over a folder containing the H5 files of TE Density output and
        return a list of absolute file paths.
        These may be used later to instantiate DensityData

        Args:
            path_to_folder (str): path to the folder containing multiple h5 files
            of density data, which are the outputs following the successful
            completion of the pipeline.

        Returns:
            raw_file_list (list of str): A list containing the absolute paths to
            each relevant H5 file of density data
        """
        raw_file_list = []  # init empty list to store filenames
        for root, dirs, files in os.walk(path_to_folder):
            for a_file_object in files:
                # N.B very particular usage of abspath and join.
                a_file_object = os.path.abspath(os.path.join(root, a_file_object))
                if a_file_object.endswith(".h5"):  # MAGIC file ext
                    raw_file_list.append(a_file_object)
        return raw_file_list

    @classmethod
    def from_list_gene_data_and_hdf5_dir(
        cls, list_of_gene_data, HDF5_folder, file_substring, logger
    ):
        """
        Returns a list of DensityData instances when given a list of GeneData
        and a directory of H5 files

        Args:
            list_of_gene_data (list): List of GeneData instances, each GeneData
                represents the information of a single chromosome, the
                chromosome ID must match with the string ID of the chromosome
                which is given by the .h5 file.

            HDF5_folder (str): Path to directory containing results (.h5) files
            following the successful computation of TE Density

            file_substring (str): MAGIC substring with which to identify the
                genome and chromosome IDs. Generally, substring should
                be "GenomeName_(.*?).h5" so that it can correctly grab the
                chromosome ID from the filename

            logger (logging.logger): Obj to log information to

        Returns: processed_density_data (list of DensityData instances).
        """
        processed_density_data = []
        for raw_hdf5_data_file in DensityData.supply_density_data_files(HDF5_folder):
            current_hdf5_file_chromosome = re.search(file_substring, raw_hdf5_data_file)
            if current_hdf5_file_chromosome is None:
                logger.warning(
                    """Using your regex pattern: %s, unable to identify
                    chromosome IDs from directory: %s. Please refer to the
                    documentation of the classmethod
                    'from_list_gene_data_and_hdf5_dir' for more information"""
                    % (file_substring, HDF5_folder)
                )
            current_hdf5_file_chromosome = current_hdf5_file_chromosome.group(
                1
            )  # MAGIC to extract info from
            # regex search obj
            for gene_data_obj in list_of_gene_data:
                if gene_data_obj.chromosome_unique_id == current_hdf5_file_chromosome:
                    dd_data_obj = cls.verify_h5_cache(
                        raw_hdf5_data_file, gene_data_obj, logger
                    )
                    processed_density_data.append(dd_data_obj)
        return processed_density_data

    def info_of_gene(self, gene_id, window_idx, n_te_types=5):
        """
        gene_id (str): String representing the name of the gene to report
        information

        window_idx (int): Integer representing the index of the window that you
        want to display information for, the smallest window starts at 0.

        n_te_types (int): Defaults to 5, integer representing how many TE types
        to show for the greatest and least values.
        """
        gene_index = self._index_of_gene(gene_id)
        window_val = self.window_list[window_idx]

        # TODO candidate to clean up below, utility function for those not
        # well-versed in indexing HDF5

        # NB sort in reverse order, the last items in this array contains
        # the greatest density value
        # ORDERS
        sorted_order_left_indices = np.argsort(
            self.left_orders[:, window_idx, gene_index]
        )
        sorted_order_intragenic_indices = np.argsort(
            self.intra_orders[:, 0, gene_index]
        )
        sorted_order_right_indices = np.argsort(
            self.right_orders[:, window_idx, gene_index]
        )

        # SUPERFAMILIES
        sorted_super_left_indices = np.argsort(
            self.left_supers[:, window_idx, gene_index]
        )
        sorted_super_intragenic_indices = np.argsort(
            self.intra_supers[:, 0, gene_index]
        )
        sorted_super_right_indices = np.argsort(
            self.right_supers[:, window_idx, gene_index]
        )

        # NB sorted array, containing actual values now
        # ORDERS
        sorted_order_left = self.left_orders[:, window_idx, gene_index][
            sorted_order_left_indices
        ]
        sorted_order_intra = self.intra_orders[:, 0, gene_index][
            sorted_order_intragenic_indices
        ]
        sorted_order_right = self.right_orders[:, window_idx, gene_index][
            sorted_order_right_indices
        ]
        # SUPERFAMILIES
        sorted_super_left = self.left_supers[:, window_idx, gene_index][
            sorted_super_left_indices
        ]
        sorted_super_intra = self.intra_supers[:, 0, gene_index][
            sorted_super_intragenic_indices
        ]
        sorted_super_right = self.right_supers[:, window_idx, gene_index][
            sorted_super_right_indices
        ]

        info = f"""
        -------------------------------------------------
        Gene Index and ID: {gene_index, gene_id}
        Window Index and Value: {window_idx, window_val}
        Chosen N TE Types to Display: {n_te_types}
        TE Orders in Annotation: {self.order_list}
        TE SuperFamilies in Annotation: {self.super_list}
        NOTE, there can be multiple entries with 0 as a density value, so if
        zero are returned, there very well might be more.
        -------------------------------------------------

        TOP {n_te_types} TE ORDERS:
            Upstream:
                {np.array(self.order_list)[sorted_order_left_indices[-n_te_types:]]}
                {sorted_order_left[-n_te_types:]}
            Intragenic:
                {np.array(self.order_list)[sorted_order_intragenic_indices[-n_te_types:]]}
                {sorted_order_intra[-n_te_types:]}
            Downstream:
                {np.array(self.order_list)[sorted_order_right_indices[-n_te_types:]]}
                {sorted_order_right[-n_te_types:]}

        BOTTOM {n_te_types} TE ORDERS:
            Upstream:
                {np.array(self.order_list)[sorted_order_left_indices[:n_te_types]]}
                {sorted_order_left[:n_te_types]}
            Intragenic:
                {np.array(self.order_list)[sorted_order_intragenic_indices[:n_te_types]]}
                {sorted_order_intra[:n_te_types]}
            Downstream:
                {np.array(self.order_list)[sorted_order_right_indices[:n_te_types]]}
                {sorted_order_right[:n_te_types]}

        -------------------------------------------------

        TOP {n_te_types} SUPERFAMILIES:
            Upstream:
                {np.array(self.super_list)[sorted_super_left_indices[-n_te_types:]]}
                {sorted_super_left[-n_te_types:]}
            Intragenic:
                {np.array(self.super_list)[sorted_super_intragenic_indices[-n_te_types:]]}
                {sorted_super_intra[-n_te_types:]}
            Downstream:
                {np.array(self.super_list)[sorted_super_right_indices[-n_te_types:]]}
                {sorted_super_right[-n_te_types:]}

        BOTTOM {n_te_types} SUPERFAMILIES:
            Upstream:
                {np.array(self.super_list)[sorted_super_left_indices[:n_te_types]]}
                {sorted_super_left[:n_te_types]}
            Intragenic:
                {np.array(self.super_list)[sorted_super_intragenic_indices[:n_te_types]]}
                {sorted_super_intra[:n_te_types]}
            Downstream:
                {np.array(self.super_list)[sorted_super_right_indices[:n_te_types]]}
                {sorted_super_right[:n_te_types]}

        -------------------------------------------------
        """
        return info

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
