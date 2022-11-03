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

from transposon.gene_data import GeneData

DensitySlice = namedtuple(
    "DensitySlice", ["slice", "direction", "window_val", "te_type"]
)


# TODO refactor to use the MergeData
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

        self.gene_list = [
            gene.decode("utf-8") for gene in self.data_frame["GENE_NAMES"][:]
        ]
        self.num_genes = len(self.gene_list)
        self.chromosomes = [
            chromosome.decode("utf-8")
            for chromosome in self.data_frame["CHROMOSOME_ID"][:]
        ]
        self.unique_chromosomes = list(set(self.chromosomes))
        if len(self.unique_chromosomes) != 1:
            raise ValueError(
                "There are multiple unique chromosomes in this density data."
            )
        self.unique_chromosome_id = self.unique_chromosomes[0]  # MAGIC

        self.windows = self.data_frame["WINDOWS"]  # not int
        self.window_list = [int(i) for i in self.windows[:]]  # list of ints
        self.order_list = [
            order.decode("utf-8") for order in self.data_frame["ORDER_NAMES"][:]
        ]
        self.super_list = [
            superfam.decode("utf-8")
            for superfam in self.data_frame["SUPERFAMILY_NAMES"][:]
        ]

        # NB. Shape for these is (type of TE, window, gene)
        # NB these access the subarrays of the HDF5
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
        return self.gene_list.index(gene_string)

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
                """The supplied TE name string : '%s', must be found
                within the following list: %s, in order to subset the
                arrays appropriately."""
                % (te_name, acceptable_te_name_list)
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

    # NOTE, this is also a helper function, but it is used in a scheme to
    # instantiate multiple DensityData objects. Ask Mike if it should be
    # refactored to the utils script.
    @staticmethod
    def _supply_density_data_files(path_to_folder, pattern=".h5"):
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
        # TODO modify this function name and description now that I am also
        # using it for tsvs
        raw_file_list = []  # init empty list to store filenames
        for root, dirs, files in os.walk(path_to_folder):
            for a_file_object in files:
                # N.B very particular usage of abspath and join.
                a_file_object = os.path.abspath(os.path.join(root, a_file_object))
                if a_file_object.endswith(pattern):  # MAGIC file ext
                    raw_file_list.append(a_file_object)
        return sorted(raw_file_list)

    # TODO, ask Michael if this should be moved. I don't think so?
    @classmethod
    def from_list_genedata_dir_and_hdf5_dir(
        cls, genedata_cache_dir, HDF5_folder, logger
    ):
        """
        Returns a list of DensityData instances when given a list of GeneData
        and a directory of H5 files


        Args:
            genedata_cache_dir (list): Path to directory containing the cached
            genedata files that were made during the pipeline. These are
            basically .tsv files corresponding to each pseudomolecule of the
            user's gene annotation

            HDF5_folder (str): Path to directory containing results (.h5) files
            following the successful computation of TE Density

            logger (logging.logger): Obj to log information to

        Returns: processed_dd_data (list of DensityData instances).
        """
        # Sort the filepaths for the GeneDatas, this will make it so that
        # the pseudmolecules are sorted
        genedata_abs_paths = sorted(
            DensityData._supply_density_data_files(genedata_cache_dir, "GeneData.tsv")
        )
        # Initialize gene data from the sorted filepaths, now I have 'sorted'
        # GeneData objects
        list_of_gene_data = [GeneData.read(filepath) for filepath in genedata_abs_paths]

        # Sort the filepaths for the h5 data.
        all_unprocessed_h5_files = sorted(
            DensityData._supply_density_data_files(HDF5_folder)
        )

        if len(list_of_gene_data) != len(all_unprocessed_h5_files):
            logger.critical(
                """
                Your gene data cache files: %s do not correspond 1:1 to your
                output from TE Density: %s

                This function only looks for '*GeneData.tsv' files in a
                directory, as well as the '.h5' (NOT overlap) output files.
                It could help to make sure there are no extraneous files in
                each directory.
            """
                % (genedata_abs_paths, all_unprocessed_h5_files)
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

    # NOTE, this is legacy code. It has been replaced with
    # from_list_genedata_dir_and_hdf5_dir() which I think is more clean and
    # easier to interact with. Users kept getting tripped up with the regex
    # pattern and I forgot that users will have genomes with both 'chr_xyz' and
    # 'contig24_abc_xyz' as their pseudmolecule IDs, so a single regex pattern would
    # be hard to implement at times. The replacement class method,
    # from_list_genedata_dir_and_hdf5_dir(), foregos the usage of a regex
    # pattern and simply alphabetically sorts the filenames, this makes it easy
    # to match GeneData to HDF5 files.

    # More to come, Michael and I are still working on refactoring DensityData
    # and MergeData.
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

    def __repr__(self):
        """String representation for developer."""
        info = """
                DensityData Genome ID: {self.genome_id}
                DensityData Chromosome ID: {self.unique_chromosome_id}
                DensityData shape key: (type of TE, windows, genes)
                DensityData left_order shape: {self.left_orders.shape}
                DensityData intra_order shape: {self.intra_orders.shape}
                DensityData right_order shape: {self.right_orders.shape}
                DensityData left_supers shape: {self.left_supers.shape}
                DensityData intra_supers shape: {self.intra_supers.shape}
                DensityData right_supers shape: {self.right_supers.shape}
            """
        return info.format(self=self)
