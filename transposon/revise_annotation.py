#!/usr/bin/env python3

"""
Create a new transposon annotation file, one without overlapping TEs.
Overlapping TEs are only aggregated into one distinct if they are of the same
type.

This code is used in the context of the TE Density pipeline to ameliorate
the issue of double-counting base pairs when calculating TE density in a given
window. Often, in regular TE annotations, TEs can be fragmented and quite often
overlap with one another. While this may be biologically valid, it creates an
issue of accurately reporting the amount of base pairs in a given window, and
when calculating density (the number of base-pairs of TE / the base-pair window
size) a value can be reached that is over 1.00; similarly values under 1.00
could still be inflated. To this end, we only merge TEs that have the same identity
(Order / SuperFamily designations).

Finally, while this file is not an entry point for the user, they may elect to
not use a revised annotation. They can do so by supplying an optional
command-line argument to density.py. However by not using this file, the user
will be returned density values that are inflated, and in some cases over 1.00,
which can obfuscate interpretation of TE density results.
"""

__author__ = "Scott Teresi"

import pandas as pd
from tqdm import tqdm
import os


class Revise_Anno(object):
    """
    Contains methods to create an altered transposon annotation file (GFF
    format) so that TEs do not overlap with one another. A new TE entry is
    created that is the aggregation of the lowest Start value and the highest
    Stop value of the aggregated TEs.
    """

    def __init__(self, transposon_anno, h5_cache_loc,
                 logger, genome_name):
        """
        Initialize

        Args:
            transposon_anno (pandas.core.DataFrame):

            self.h5_cache_loc ():

            logger ():
        """
        self.transposon_data = transposon_anno
        self.logger = logger
        self.superfam_cache_loc = os.path.join(h5_cache_loc,
                                               (genome_name +
                                               '_superfam_revision_cache.h5'))
        self.order_cache_loc = os.path.join(h5_cache_loc,
                                               '_',
                                               genome_name,
                                               '_order_revision_cache.h5')
        self.nameless_cache_loc = os.path.join(h5_cache_loc,
                                               '_',
                                               genome_name,
                                               '_nameless_revision_cache.h5')
        self.updated_te_annotation = None

    def create_superfam(self):
        self._reset_starting_vals()
        self.transposon_data['Order'] = 'Revision'
        self.iterate_call_merge('SuperFamily')
        self._write(self.updated_te_annotation, self.superfam_cache_loc)

    def create_order(self):
        self._reset_starting_vals()
        self.transposon_data['SuperFamily'] = 'Revision'
        self.iterate_call_merge('Order')
        self._write(self.updated_te_annotation, self.order_cache_loc)

    def create_nameless(self):
        self._reset_starting_vals()
        # Change the Order and SuperFamily column vals
        self.transposon_data['Order'] = 'Total_TE_Density'
        self.transposon_data['SuperFamily'] = 'Total_TE_Density'
        # Call the constructors for nameless groups and save result
        self.iterate_call_merge('Nameless')
        self._write(self.updated_te_annotation, self.nameless_cache_loc)

    def _merge_all(self):
        self._verify_files()
        supers = self._read(self.superfam_cache_loc)
        orders = self._read(self.order_cache_loc)
        nameless = self._read(self.nameless_cache_loc)
        self.whole_te_annotation = pd.concat([supers,
                                              orders,
                                              nameless])
        self.whole_te_annotation.sort_values(by=['Start', 'Chromosome'],
                                             inplace=True)

    def save_updated_te_annotation(self, filename, header=True, index=False):
        """
        Save the annotation.

        Args:
            filename (str)
        """
        self._merge_all()
        self.whole_te_annotation.to_csv(filename, sep='\t', header=header,
                                        index=index)
        self.logger.info('Annotation has been saved...')

    def _verify_files(self):
        """
        This exists because it is possible someone may try to access
        self.whole_te_annotation without having run each of the 3 creation
        commands (in reality those commands are run outside of human use,
        inside verify_cache.py), but just to be sure that all 3 creation
        commands were run so that self.whole_te_annotation is correct, I wrote
        this verify command. Long story short, this check is here to make sure
        we can appropriately merge the revised elements to create a complete
        revised annotation.
        """
        # Superfamilies
        if os.path.exists(self.superfam_cache_loc):
            return
        else:
            msg = ('Revised SuperFamily H5 file DNE: {self.superfam_cache_loc}')
            self.logger.critical(msg.format(self=self))
            raise FileNotFoundError(msg.format(self=self))

        # Orders
        if os.path.exists(self.order_cache_loc):
            return
        else:
            msg = ('Revised Order H5 file DNE: {self.order_cache_loc}')
            self.logger.critical(msg.format(self=self))
            raise FileNotFoundError(msg.format(self=self))

        # Nameless
        if os.path.exists(self.nameless_cache_loc):
            return
        else:
            msg = ('Revised Nameless H5 file DNE: {self.nameless_cache_loc}')
            self.logger.critical(msg.format(self=self))
            raise FileNotFoundError(msg.format(self=self))

    def _reset_starting_vals(self):
        self.chrom_specific_frame_dict = None
        self.updated_te_annotation = None
        self.complete_chromosomes_dict = {}
        self.seed_frame = None
        self.search_frame = None

    def iterate_call_merge(self, te_type_to_split):
        """
        This is an entry point, define self variables and initiate the search
        space functions.

        Create the updated TE annotation for each TE grouping. Concatenate all
        of them when done to produce the finalized version of the transposon
        annotation. Iterates over the separate TE frames
        """
        with tqdm(total=len(self.chromosome_groups), position=0, ncols=80,
                  desc='Revising a chromosome for ' + te_type_to_split) as pbar:
            print()

            if te_type_to_split == 'Nameless':
                te_type_to_split = 'Order'  # Can be Order or Superfamily,
                # doesn't matter, just wanted to have nameless for the tqdm
                # object above.
            for chromosome_of_data in self.chromosome_groups:
                self.chrom_specific_frame_dict = {}
                for te_frame in self.te_groupings(chromosome_of_data,
                                                  te_type_to_split):
                    chromosome = te_frame.Chromosome.unique()[0]  # Magic number
                    self.current_te_identity = te_frame[te_type_to_split].unique()[0]  # Magic number

                    self.seed_frame = te_frame.copy(deep=True).sort_values(by=['Start'])
                    self.search_frame = te_frame.copy(deep=True).sort_values(by=['Start'])
                    self.seed_max_index = int(self.seed_frame.index.values.max())
                    self.chrom_specific_frame_dict[self.current_te_identity] = pd.DataFrame()

                    self.call_merge()  # NOTE function call, recursion starts

                pbar.update(1)
                print()  # want to make the progress manager clearer
                self.logger.info('Done with chromosome ' + chromosome + '...')
                self.concat_single_chrom(chromosome)
            self.concat_all_complete_chrom()

    @staticmethod
    def split(dataframe, group):
        """Return list of dataframes with each element being a subset of the df.

        I use this function to split by chromosome so that we may later do
        chromosome element-wise operations.

        This function is also used in revise_annotation.py to split on transposon
        identities.
        """
        grouped_df = dataframe.groupby(group)
        return [grouped_df.get_group(x) for x in grouped_df.groups]

    @staticmethod
    def te_groupings(dataframe, te_type):
        return Revise_Anno.split(dataframe, te_type)

    @property
    def chromosome_groups(self):
        """
        Split by chromosome
        """
        return self.split(self.transposon_data, 'Chromosome')

    def call_merge(self):
        """
        Grab one row at a time (1 TE) from the seed frame and call the merging
        function.
        """

        if self.seed_frame.empty:
            # When the dataframe is out of elements, exit the loop of calls
            return
        else:
            seed_element = self.seed_frame.iloc[0].to_frame().T  # Magic number
            seed_idx = seed_element.index.values[0]  # Magic number
            self.merge_by_like(seed_idx, seed_element)  # NOTE function call


    @staticmethod
    def adjust_length(panda_dataframe):
        """
        Redefine the length value for TEs, the value is no longer correct for a
        given TE after aggregating TEs.

        Args:
            panda_dataframe (pandas.core.DataFrame): A pandas dataframe that
            represents the TE data. TE lengths are incorrect.

        Returns:
            panda_dataframe (pandas.core.DataFrame): A pandas dataframe that
            represents the TE data. TE lengths are correct.
        """
        panda_dataframe['Length'] = panda_dataframe['Stop'] - panda_dataframe['Start'] + 1
        return panda_dataframe

    def concat_all_complete_chrom(self):
        """
        Concatenates all of the individual chromosomes for a given genome into
        one pandas.core.DataFrame. Each individual chromosome has previously
        been concatenated together elsewhere.
        """
        self.logger.info('Concatenating all chromosomes...')
        to_concat = [te_anno_dataframe.sort_values(by=['Start']) for
                     chromosome_id, te_anno_dataframe in
                     self.complete_chromosomes_dict.items()]
        self.updated_te_annotation = pd.concat(to_concat, ignore_index=True)

    @staticmethod
    def save_for_dev(pandaframe, filename):
        pandaframe.to_csv(filename, sep='\t', header=True, index=False)

    def concat_single_chrom(self, chromosome):
        """
        Loop over the dictionary of pandadataframes that are on a type-by-type
        basis and concatenate them all into the new transposon annotation.

        Concatenates all of the chromosome groupings for ONE chromosome

        """
        self.logger.info('Concatenating chromosome ' + chromosome + '...')
        to_concat = [te_anno_dataframe for te_id, te_anno_dataframe in
                     self.chrom_specific_frame_dict.items()]
        # Concatenate all of the data that was previously collected on a
        # TE-type specific basis
        full_chrom = pd.concat(to_concat, ignore_index=True)
        full_chrom = self.adjust_length(full_chrom)
        # Place the completed chromosome in a dictionary
        self.complete_chromosomes_dict[chromosome] = full_chrom

    def merge_by_like(self, seed_idx, seed_row, to_drop=True):
        """
        Given a single TE (seed row), identify and merge all other TEs that
        happen to overlap with the seed TE. Add the updated row to a new
        dataframe and move on to the next element for the search space.

        Args:
            seed_idx  (numpy.int64): The integer of the seed element's index.

            seed_row (pandas.core.frame.DataFrame): A Pandas DataFrame
            containing one row only. This one row is the current seed element.

            to_drop (bool): Whether or not to drop the indices of the seed_row
            from the seed & search frame. to_drop defaults to True on the first
            run through (so the seed element is always dropped from the seed
            and search frame), however during the recursive search, after the
            seed element has been updated to have the new stop, we set it to
            false, we do cannot drop that element again, and instead just do
            the search.
        """
        seed_start = int(seed_row.Start.values[0])
        seed_stop = int(seed_row.Stop.values[0])

        self.search_start = self.search_frame.Start.to_numpy(copy=False)
        self.search_stop = self.search_frame.Stop.to_numpy(copy=False)

        if to_drop:
            # Remove seed element from search and seed array
            # This always occurs on the first run
            self.search_frame.drop(index=[seed_idx], inplace=True)
            self.seed_frame.drop(index=[seed_idx], inplace=True)

        # Generate the array of indices of hits
        hit_scan_overlap_array = self.hit_scan_overlapping(seed_start, seed_stop)

        if len(hit_scan_overlap_array) > 0:  # there are hits, that are non-seed
            seed_stop = self.determine_seed_stop(seed_stop, hit_scan_overlap_array)
            # Remove the hits (the rows) from the search frame.
            self.search_frame = Revise_Anno.clear_array_by_index(self.search_frame,
                                                                 hit_scan_overlap_array)
            # Remove the hits (the rows) from the seed frame.
            self.seed_frame = Revise_Anno.clear_array_by_index(self.seed_frame,
                                                               hit_scan_overlap_array)

            # seed_row.Stop = seed_stop NOTE this was old, change?
            seed_row = Revise_Anno.set_seed_stop(seed_row, seed_stop)

            # Call search again with updated information

            # NOTE recursive function call here.
            # I do this because the updated TE (new stop position) may have new
            # hits itself, so we have to recursively keep searching until 0
            # hits.
            self.merge_by_like(seed_idx, seed_row, to_drop=False)

        elif len(hit_scan_overlap_array) == 0:
            # Once there are no more hits we go about adding the updated final
            # element into the dataframe
            # Add the row in its virgin or final updated form to the dataframe
            seed_row = Revise_Anno.set_seed_stop(seed_row, seed_stop)
            self.update_data_frame(seed_row)

            # Seed element has laready been dropped from search and seed array
            # Get ready to move on to the next seed element
            if seed_idx == self.seed_max_index:
                # Return because we are done with all elements
                return
            else:
                # Move to next element
                self.call_merge()

    def hit_scan_overlapping(self, seed_start, seed_stop):
        """
        Produce an array of indices from the search array where each index is
        the location of an element that satisfies both of the search
        parameters. I define the search parameters as a TE that has a start
        inside seed TE start and a stop value that is outside of the seed TE.

        Args:
            seed_start (int): The integer of the seed start value

            seed_stop (int): The integer of the seed stop value

        Returns:
            hit_scan_overlap_array (list):
        """
        hit_scan_overlap_array = self.search_frame[(self.search_frame.Start
                                                    >= seed_start) &
                                                   (self.search_frame.Start
                                                    <=
                                                    seed_stop)].index.to_list()
        return hit_scan_overlap_array

    @staticmethod
    def set_seed_stop(seed_row, new_seed_stop):
        """
        Set the seed stop value

        Args:
            seed_row (pandas.core.frame.DataFrame): A Pandas DataFrame
            containing one row only. This one row is the current seed element.

            new_seed_stop (int): The new value for the seed stop.
        """
        seed_row.Stop = new_seed_stop
        seed_row = seed_row.astype({'Start': 'int32', 'Stop': 'int32'})
        return seed_row

    def update_data_frame(self, seed_row):
        """
        Add the updated entry (seed_row) to the dataframe.

        Args:
            seed_row (pandas.core.frame.DataFrame): A Pandas DataFrame
            containing one row only. This one row is the current seed element.
            It should have been previously updated to possess the correct stop
            value.
        """
        to_concat = [self.chrom_specific_frame_dict[self.current_te_identity],
                     seed_row]
        updated_frame = pd.concat(to_concat)
        self.chrom_specific_frame_dict[self.current_te_identity] = updated_frame

    def determine_seed_stop(self, seed_stop, array_of_hits):
        """
        Identify the maximum seed stop value

        Args:
            seed_stop (int): The value for the seed stop.

            array_of_hits (list-like): An array of indices representing rows
            that overlapped with the seed element.

        Returns:
            seed_stop (int): The value for the seed stop.
        """
        for hit_index in array_of_hits:
            if self.search_frame.loc[hit_index, ].Stop > seed_stop:
                seed_stop = self.search_frame.loc[hit_index, ].Stop
        return seed_stop

    @staticmethod
    def clear_array_by_index(dataframe, array_of_hits):
        """
        Remove rows from a given array. Done to remove hits so that proceeding
        searches don't hit the same element multiple times.

        Args:
            dataframe (pandas.core.frame.DataFrame): A data to remove hits from.
            Can be both the search frame or the seed frame.

            array_of_hits (list-like): An array of indices representing rows
            that overlapped with the seed element.

        Returns:
            dataframe(pandas.core.frame.DataFrame): The updated dataframe,
            indices have been removed.
        """
        dataframe.drop(index=array_of_hits, axis=0, inplace=True)
        return dataframe

    @staticmethod
    def _write(dataframe, filename, key='default'):
        """Write a Pandaframe to disk.

        Args:
            dataframe ():
            filename (str): a string of the filename to write.
            key (str): identifier for the group (dataset) in the hdf5 obj.
        """
        dataframe.to_hdf(filename, key=key, mode='w')

    @staticmethod
    def _read(filename, key='default'):
        """Read from disk. Returns a Pandaframe from an hdf5 file

        Args:
            filename (str): a string of the filename to write.
            key (str): identifier for the group (dataset) in the hdf5 obj.
        """
        return pd.read_hdf(filename, key=key)
