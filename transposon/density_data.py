#!/usr/bin/env/python

"""
Read TE Density H5 files
"""

__author__ = "Scott Teresi"

import h5py
import numpy as np
import pandas as pd
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
        # MAGIC must have '.h5' as extension
        new_filename = input_h5.replace(".h5", "_SenseSwapped.HDF5")
        shutil.copyfile(input_h5, new_filename)  # copy because we
        # need to swap values later
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
        # NB the 0s are the antisense

        if sense_swap:
            zero_gene_list = np.where(strands_as_numpy == 0)[0]  # gets indices of 0s
            genes_to_swap = gene_data.data_frame.iloc[zero_gene_list, :].index.tolist()
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

    @property
    def window_index_dict(self):
        """Returns a dictionary of window ints as keys and indices as
        values"""
        window_dict = {}
        for i in range(len(self.window_list)):
            window_dict[self.window_list[i]] = i
        return window_dict

    def _verify_direction_string(self, direction):
        """
        Verify whether

        """
        acceptable_directions = ["Upstream", "Intra", "Downstream"]
        if direction not in acceptable_directions:
            raise ValueError(
                """The supplied direction string: %s, must be found
                within the following list: %s, in order to subset the
                             arrays appropriately."""
                % (direction, acceptable_directions)
            )

    def _verify_te_category_string(self, te_category):
        # TODO get docstring
        acceptable_te_categories = ["Order", "Superfamily"]
        if te_category not in acceptable_te_categories:
            raise ValueError(
                """The supplied te_category string: %s, must be found
                within the following list: %s, in order to subset the
                arrays appropriately."""
                % (te_category, acceptable_te_categories)
            )

    def _verify_window_val(self, direction, window_val):
        # TODO get docstring
        if direction == "Intra" and window_val is not None:
            raise ValueError(
                """The user requested the corresponding dataset
                             for 'Intra' values, and supplied '%s' for the
                             'window_val' argument. The 'window_val' argument
                             ought to be left to its default arg of 'None', as
                             there is no window for intra values. Raising an
                             error to notify user of unintended usage."""
                % window_val
            )
        elif direction == "Intra" and window_val is None:
            return

        acceptable_window_values = self.window_list
        if window_val not in acceptable_window_values:
            raise ValueError(
                """The supplied window value: %s, must be found
                within the following list: %s, in order to subset the
                upstream and downstream arrays appropriately."""
                % (window_val, acceptable_window_values)
            )

    def _verify_te_name(self, te_category, te_name):
        # TODO get docstring
        if te_category == "Order":
            acceptable_te_name_list = self.order_list
        if te_category == "Superfamily":
            acceptable_te_name_list = self.super_list
        if te_name not in acceptable_te_name_list:
            raise ValueError(
                """The supplied TE name string : %s, must be found
                within the following list: %s, in order to subset the
                arrays appropriately."""
                % (direction, acceptable_te_name_list)
            )

    @staticmethod
    def _verify_unique_chromosome_pandaframe(gene_info_pandas, chrom_col):
        # MAGIC to check chromosome identity
        if len(gene_info_pandas[chrom_col].unique()) > 1:
            raise ValueError(
                """Your gene pandaframe has too many unique
                             pseudomolecule in it, it should only have genes
                             belonging to one pseudomolecule."""
            )

    def _verify_self_chromosome_match_w_pandaframe(self, gene_info_pandas, chrom_col):
        if self.unique_chromosome_id != gene_info_pandas[chrom_col].unique()[0]:
            raise ValueError(
                """The pseudomolecule of DensityData instance does
                    not match pseudomolecule of the supplied pandas
                    object."""
            )

    def add_hdf5_indices_to_gene_data(
        self,
        gene_info_pandas,
        gene_name_col="Gene_Name",
        chrom_col="Chromosome",
        index_col="Index_Val",
    ):
        """
        Take a pandas dataframe of genes (rows) and pseudomolecule identities
        and add a column that contains each gene's index value in the
        DensityData (HDF5) dataset.

        Args:
            gene_info_pandas (pandas.core.frame.DataFrame): A pandas dataframe
                with AT LEAST the three columns defined as the default args:
                gene_name_col, chrom_col, index_col.

            gene_name_col (str): A string that represents the column that
                contains the genes in the users pandas dataframe.

            chrom_col (str): A string that represents the column that
                contains the identities of the pseudomolecules in the users
                pandas dataframe.

            index_col (str): A string that represents the column that contains
                the integer values of the index values for each gene in the
                HDF5. The column name defaults to 'Index_Val'.

        Returns:
            gene_info_pandas (pandas.core.frame.DataFrame): Returns the
                original dataframe, but with a new column. The column's
                name is supplied by the argument 'index_col'. The column
                contains integer values of the index values for each gene in
                the HDF5.
        """
        DensityData._verify_unique_chromosome_pandaframe(gene_info_pandas, chrom_col)
        self._verify_self_chromosome_match_w_pandaframe(gene_info_pandas, chrom_col)

        # NB get around setting with copy warning by creating a deep copy,
        # not an issue for performance.
        gene_info_pandas = gene_info_pandas.copy(deep=True)
        gene_info_pandas[index_col] = gene_info_pandas.apply(
            lambda x: self._index_of_gene(x[gene_name_col]), axis=1
        )
        return gene_info_pandas

    def add_te_vals_to_gene_info_pandas(
        self,
        gene_info_pandas,
        te_category,
        te_name,
        direction,
        window_val,
        gene_name_col="Gene_Name",
        chrom_col="Chromosome",
        index_col="Index_Val",
    ):
        """
        Take a pandas dataframe that already has HDF5 index values for each
        gene (row), and add column to the dataframe that has TE Density values
        for each gene (row) according to a user-supplied TE, direction, and
        window value. This step comes after the method 'add_hdf5_indices_to_gene_data'.

        Args:
            gene_info_pandas (pandas.core.frame.DataFrame): A pandas dataframe
                with AT LEAST the three columns defined as the default args:
                gene_name_col, chrom_col, index_col.

            te_category (str): A string of a TE category, must be either
                'Order' or 'Superfamily'.

            te_name (str): A string representing the valid name of a TE group
                that is in the HDF5.

            direction (str): A string representing whether the user wants TE
                Density data for 'Upstream' or 'Downstream'. Must be either
                'Upstream' or 'Downstream'.

            window_val (int): An integer representing a valid window value that
                the user wants the TE Density data for

            gene_name_col (str): A string that represents the column that
                contains the genes in the users pandas dataframe.

            chrom_col (str): A string that represents the column that
                contains the identities of the pseudomolecules in the users
                pandas dataframe.

            index_col (str): A string that represents the column that contains
                the integer values of the index values for each gene in the
                HDF5. These values can be acquired using the
                'add_hdf5_indices_to_gene_data' function.

        Returns:
            gene_info_pandas (pandas.core.frame.DataFrame): Returns the
                original dataframe, but with a new column. The column's
                identifer is "te_name + '_' + window_val + '_' + direction".
                The column contains floats of TE Density values for that gene
                (row).
        """
        DensityData._verify_unique_chromosome_pandaframe(gene_info_pandas, chrom_col)
        self._verify_self_chromosome_match_w_pandaframe(gene_info_pandas, chrom_col)

        te_string_iterable = [te_name, str(window_val), direction]
        te_column_string = "_".join(te_string_iterable)

        # NB get around setting with copy warning by creating a deep copy,
        # not an issue for performance.
        gene_info_pandas = gene_info_pandas.copy(deep=True)
        gene_info_pandas[te_column_string] = gene_info_pandas.apply(
            lambda x: self.get_specific_slice(
                te_category, te_name, direction, window_val, x[index_col]
            ).slice,
            axis=1,
        )
        return gene_info_pandas

    def get_specific_slice(
        self, te_category, te_name, direction, window_val=None, gene_indices=slice(None)
    ):
        """
        Return a DensitySlice obj for a combination of TE category, TE name,
            window, direction, and optionally a set of gene indices. This
            method used in the methods 'add_hdf5_indices_to_gene_data'
            and 'add_te_vals_to_gene_info_pandas' will likely be the most
            useful to users in accessing the TE Density data.

        Args:

            te_category (str): A string of a TE category, must be either
                'Order' or 'Superfamily'.

            te_name (str): A string representing the valid name of a TE group
                that is in the HDF5.

            window_val (int): An integer representing a valid window value that
                the user wants the TE Density data for

            direction (str): A string representing whether the user wants TE
                Density data for 'Upstream' or 'Downstream'. Must be either
                'Upstream' or 'Downstream'.

            gene_indices (np.array, list of int, or slice): The indices (of
                genes) that the user wants TE Density values for. Defaults to
                slice(None) which gives the values for ALL indices (genes).

        Returns:
            DensitySlice (DensitySlice): A DensitySlice obj that contains the
                TE Density values for the combination of args that the user
                supplied. Users can access their desired output data by
                accessing the '.slice' attribute of the DensitySlice obj.
                Will be a 1D array.
        """
        self._verify_te_category_string(te_category)
        self._verify_direction_string(direction)
        self._verify_window_val(direction, window_val)
        self._verify_te_name(te_category, te_name)

        if direction == "Intra" and te_category == "Order":
            slice_to_return = self.intra_orders[
                self.order_index_dict[te_name],
                0,  # MAGIC get first index for 'window' for intra, gives back
                # 1D array
                gene_indices,
            ]

        elif direction == "Intra" and te_category == "Superfamily":
            slice_to_return = self.intra_supers[
                self.super_index_dict[te_name],
                0,  # MAGIC get first index for 'window' for intra, gives back
                # 1D array
                gene_indices,
            ]

        elif direction == "Upstream" and te_category == "Order":
            slice_to_return = self.left_orders[
                self.order_index_dict[te_name],
                self.window_index_dict[window_val],
                gene_indices,
            ]

        elif direction == "Downstream" and te_category == "Order":
            slice_to_return = self.right_orders[
                self.order_index_dict[te_name],
                self.window_index_dict[window_val],
                gene_indices,
            ]

        elif direction == "Upstream" and te_category == "Superfamily":
            slice_to_return = self.left_supers[
                self.super_index_dict[te_name],
                self.window_index_dict[window_val],
                gene_indices,
            ]

        elif direction == "Downstream" and te_category == "Superfamily":
            slice_to_return = self.right_supers[
                self.super_index_dict[te_name],
                self.window_index_dict[window_val],
                gene_indices,
            ]
        else:
            raise ValueError()

        return DensitySlice(slice_to_return, direction, window_val, te_name)

    def yield_all_slices(self):
        # TODO think of refactoring, refactor other code to use the other
        # methods.
        """
        Yields a DensitySlice (named tuple) object for each TE
        type/direction/window combination for all genes in the DensityData obj.
        """
        directions = ["Upstream", "Downstream"]
        for direction in directions:
            for window_idx, window_val in enumerate(self.window_list):
                for te_type, te_order_idx in self.order_index_dict.items():
                    if "Revision" in te_type:  # MAGIC dont use the revision
                        # set, which is an artifact of calculating TE density
                        continue
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
                    if "Revision" in te_type:  # MAGIC dont use the revision
                        # set, which is an artifact of calculating TE density
                        continue
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

            # SWAP left and right superfamily values for antisense genes
            (
                self.data_frame["RHO_SUPERFAMILIES_LEFT"][:, :, index_to_switch],
                self.data_frame["RHO_SUPERFAMILIES_RIGHT"][:, :, index_to_switch],
            ) = (
                self.data_frame["RHO_SUPERFAMILIES_RIGHT"][:, :, index_to_switch],
                self.data_frame["RHO_SUPERFAMILIES_LEFT"][:, :, index_to_switch],
            )

            # SWAP left and right order values for antisense genes
            (
                self.data_frame["RHO_ORDERS_LEFT"][:, :, index_to_switch],
                self.data_frame["RHO_ORDERS_RIGHT"][:, :, index_to_switch],
            ) = (
                self.data_frame["RHO_ORDERS_RIGHT"][:, :, index_to_switch],
                self.data_frame["RHO_ORDERS_LEFT"][:, :, index_to_switch],
            )

    @staticmethod
    def _supply_density_data_files(path_to_folder):
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
        return sorted(raw_file_list)

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
                chromosome ID from the filename. The part in the parentheses ()
                should correspond to a regex group and the
                chromosome/pseudomolecule identifier for that dataset. For more
                information please visit: https://docs.python.org/3/library/re.html#

            logger (logging.logger): Obj to log information to

        Returns: processed_dd_data (list of DensityData instances).
        """

        # NB get the list of GeneData objs and sort by chromosome ID
        list_of_gene_data = sorted(
            list_of_gene_data,
            key=lambda gene_data: gene_data.chromosome_unique_id,
            reverse=False,
        )

        # NB get the list of TE Density output files in a dir and sort
        # alphabetically
        all_unprocessed_h5_files = sorted(
            DensityData._supply_density_data_files(HDF5_folder)
        )

        # NB get the hits of files that match your regex substring provided
        chromosome_ids_unprocessed_h5_files = [
            re.search(file_substring, x) for x in all_unprocessed_h5_files
        ]

        logger.info(
            """
            Using the user's provided regex string '%s' to match file
            objects and identify the proper pseudomolecule group for
            each file. Regex group 1 of this string must correspond to
            a pseudomolecule. This is needed to initialize DensityData.
            The user should verify that the pseudomolecule IDs derived from
            the GeneData correspond to the groups derived from the filename
            of the output .h5 data.
            """
            % file_substring
        )
        for re_search_obj, gene_data in zip(
            chromosome_ids_unprocessed_h5_files, list_of_gene_data
        ):
            try:
                logger.info(
                    "Pseudomolecule from GeneData is %s, Regex group 1 of %s is %s"
                    % (
                        gene_data.chromosome_unique_id,
                        re_search_obj,
                        re_search_obj.group(1),
                    )
                )
            except (IndexError, AttributeError) as err:
                logger.critical(
                    """
                    Unable to identify chromosome IDs from files
                    matching your provided regex
                    pattern: '%s' in the directory: %s.
                    Please refer to the documentation of the classmethod
                    'from_list_gene_data_and_hdf5_dir' in transposon/density_data.py
                    for more information. Regex group 1 must correspond to each
                    .h5 file's pseudomolecule/chromosome ID. The most common
                    pattern is "GENOMENAME_(*.?).h5"
                    """
                    % (file_substring, HDF5_folder)
                )
                raise err

        # NB check if we actually were able to identify any files matching the
        # user's supplied regex pattern, raise error and message if no hits
        if not any(chromosome_ids_unprocessed_h5_files):
            logger.critical(
                """Unable to identify files matching your provided regex
                pattern: %s in the directory: %s. \n
                Please refer to the documentation of the classmethod
                'from_list_gene_data_and_hdf5_dir' in transposon/density_data.py
                for more information"""
                % (file_substring, HDF5_folder)
            )
            raise ValueError

        # MAGIC to get substring (chromosome) from regex hit object
        chromosome_ids_unprocessed_h5_files = sorted(
            [x.group(1) for x in chromosome_ids_unprocessed_h5_files]
        )

        # NB get chromosome ID from list of GeneData
        chromosome_ids_gene_data = [
            gene_data.chromosome_unique_id for gene_data in list_of_gene_data
        ]

        # NB check if chromosome ID of h5 files matches with that of the gene
        # data, fail if they don't. This is needed to initialize DensityData
        if not chromosome_ids_unprocessed_h5_files == chromosome_ids_gene_data:
            logger.critical(
                """The strings of chromosomes in your unprocessed
                hdf5 files: %s, identified using your supplied
                regex pattern: '%s', do not match the
                chromosomes in the GeneData: %s."""
                % (
                    chromosome_ids_unprocessed_h5_files,
                    file_substring,
                    chromosome_ids_gene_data,
                )
            )
            raise ValueError

        # Initialize DensityData for each pseudomolecule
        processed_dd_data = [
            cls.verify_h5_cache(raw_hdf5_data_file, gene_data_obj, logger)
            for raw_hdf5_data_file, gene_data_obj, in zip(
                all_unprocessed_h5_files, list_of_gene_data
            )
        ]
        return processed_dd_data

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

        # TODO this does not avoid the O_Revision and S_Revision datasets,
        # future release will handle these artifacts more elegantly. Unable to
        # avoid at the moment.

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
