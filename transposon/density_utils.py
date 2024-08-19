#!/usr/bin/env/python

"""
Helper functions to interrogate DensityData
"""

__author__ = "Scott Teresi"

from transposon.density_data import DensityData
from transposon.density_data import DensitySlice
import pandas as pd
import numpy as np
import logging


def add_hdf5_indices_to_gene_data_from_list_hdf5(
    cleaned_genes,
    list_processed_dd_instance,
    gene_name_col="Gene_Name",
    chrom_col="Chromosome",
    index_col="Index_Val",
):
    """
    Args:
        cleaned_genes (pandas.core.frame.DataFrame): The output of
            import_filtered_genes()

        list_processed_dd_instance (list of DensityData instances)

        gene_name_col (str): Pandas column ID for genes

        chrom_col (str): Pandas column ID for chromosomes/pseudomolecules

        index_col (str): Name for NEW Pandas column, columns is for the HDF5
            indices per gene

    Returns:
        gene_data_w_indices (pandas.core.frame.DataFrame): This is similar in
            format to the cleaned_genes argument, but has an additional column
            which is the indices of the gene inside its HDF5 file. That way the
            user can easily get the TE Density values relevant to that gene.
    """
    to_concat = []
    # MAGIC 'Chromosome' name for column
    for chrom, dataframe in cleaned_genes.groupby(chrom_col):
        # NB, reset index is done because I want Gene_Name to act as a column,
        # not as a pandas index object, it is initially read in as an index
        # though because of the import_filtered_genes function
        dataframe.reset_index(inplace=True, drop=False)
        for processed_dd_datum in list_processed_dd_instance:
            if processed_dd_datum.unique_chromosome_id == chrom:
                gene_data_w_indices = add_hdf5_indices_to_gene_data(
                    processed_dd_datum, dataframe, gene_name_col, chrom_col, index_col
                )
                to_concat.append(gene_data_w_indices)
    gene_data_w_indices = pd.concat(to_concat)
    return gene_data_w_indices


def add_te_vals_to_gene_info_pandas_from_list_hdf5(
    gene_frame_with_indices,
    list_processed_dd_instance,
    te_group,
    te_name,
    direction,
    window,
    gene_name_col="Gene_Name",
    chrom_col="Chromosome",
    index_col="Index_Val",
):
    """
    TODO
    """
    to_concat = []
    for chrom, dataframe in gene_frame_with_indices.groupby(chrom_col):
        for processed_dd_datum in list_processed_dd_instance:
            if processed_dd_datum.unique_chromosome_id == chrom:
                x = add_te_vals_to_gene_info_pandas(
                    processed_dd_datum,
                    dataframe,
                    te_group,
                    te_name,
                    direction,
                    window,
                    gene_name_col,
                    chrom_col,
                    index_col,
                )
                to_concat.append(x)
    gene_frame_w_ind_te_vals = pd.concat(to_concat)
    return gene_frame_w_ind_te_vals


# TODO do I need the gene_name_col arg?
def add_te_vals_to_gene_info_pandas(
    dd_instance,
    gene_info_pandas,
    te_category,
    te_name,
    direction,
    window_val,
    gene_name_col="Gene_Name",
    chrom_col="Chromosome",
    index_col="Index_Val",
    logger=logging.getLogger(__name__),
):
    """
    Take a pandas dataframe that already has HDF5 index values for each
    gene (row), and add column to the dataframe that has TE Density values
    for each gene (row) according to a user-supplied TE, direction, and
    window value. This step comes after the method 'add_hdf5_indices_to_gene_data'.

    Args:
        dd_instance (DensityData)

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

    OPTIONAL Args:
        gene_name_col (str): Column ID for the genes in the user's pandas
            dataframe.

        chrom_col (str): Column ID for the chromosomes in the user's pandas
            dataframe

        index_col (str): Column ID fo the HDF5 indices of each gene in the
            user's pandas dataframe. Values can be acquired through the method
            add_hdf5_indices_to_gene_data_from_list_hdf5

    Returns:
        gene_info_pandas (pandas.core.frame.DataFrame): Returns the
            original dataframe, but with a new column. The column's
            identifier is "te_name + '_' + window_val + '_' + direction".
            The column contains floats of TE Density values for that gene
            (row).
    """
    # Set our new column name
    te_string_iterable = [te_name, str(window_val), direction]
    te_column_string = "_".join(te_string_iterable)

    # NB get around setting with copy warning by creating a deep copy,
    # not an issue for performance.
    gene_info_pandas = gene_info_pandas.copy(deep=True)

    # Do a bunch of verification steps to make sure the user didn't provide
    # inputs that aren't valid
    dd_instance._verify_te_category_string(te_category)
    dd_instance._verify_direction_string(direction)
    dd_instance._verify_window_val(direction, window_val)

    # Do more verification related to the pandas dataframe this time
    verify_uniq_chrom_pandaframe(gene_info_pandas, chrom_col)
    verify_chromosome_match_w_pandaframe(dd_instance, gene_info_pandas, chrom_col)

    # Check to see if the requested TE is indeed in the HDF5
    # This warning will occur a lot if people have fragmented TE annotations
    # where a certain TE type is not present on all psuedomolecules/scaffolds.
    # We don't want to raise an error here and force a crash but we do want to
    # communicate to the user that they are asking for a TE that is not
    # available. Since this function will often be used to generate one big
    # table for the entire genome, we will add a column of 0 to the table
    try:
        dd_instance._verify_te_name(te_category, te_name)
    except ValueError:
        logger.warning(
            f"""
            Logging warning to user:

            User has asked for a TE name: ({te_category} {te_name}) that is
            not in the HDF5 for {dd_instance.unique_chromosome_id}.
            Window: {window_val}, Direction: {direction}.

            Setting that column's TE Density to 0 for those genes on that scaffold.
            """
        )
        gene_info_pandas[te_column_string] = 0
        print(gene_info_pandas[te_column_string])
        return gene_info_pandas

    # Add the column of TE values for each gene to the pandas dataframe
    gene_info_pandas[te_column_string] = gene_info_pandas.apply(
        lambda x: get_specific_slice(
            dd_instance, te_category, te_name, direction, window_val, x[index_col]
        ).slice,
        axis=1,
    )
    return gene_info_pandas


def get_specific_slice(
    dd_instance,
    te_category,
    te_name,
    direction,
    window_val=None,
    gene_indices=slice(None),
):
    """
    Return a DensitySlice obj for a combination of TE category, TE name,
        window, direction, and optionally a set of gene indices. This
        method used in the methods 'add_hdf5_indices_to_gene_data'
        and 'add_te_vals_to_gene_info_pandas' will likely be the most
        useful to users in accessing the TE Density data.

    Args:
        dd_instance (DensityData)

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
    # Do a bunch of verification steps to make sure the user didn't provide
    # inputs that aren't valid, NOTE this doubles up on verifying with
    # add_te_vals_to_gene_info_pandas()
    dd_instance._verify_te_category_string(te_category)
    dd_instance._verify_direction_string(direction)
    dd_instance._verify_window_val(direction, window_val)

    # TODO consider adding a try except block here to catch the ValueError and
    # return a NaN array. But since this function is expected to be used
    # in a limited fasion on one specific chromosome, I will leave it as is for
    # now. Scott 05-08-2024
    dd_instance._verify_te_name(te_category, te_name)

    if direction == "Intra" and te_category == "Order":
        slice_to_return = dd_instance.intra_orders[
            dd_instance.order_index_dict[te_name],
            0,  # MAGIC get first index for 'window' for intra, gives back
            # 1D array
            gene_indices,
        ]

    elif direction == "Intra" and te_category == "Superfamily":
        slice_to_return = dd_instance.intra_supers[
            dd_instance.super_index_dict[te_name],
            0,  # MAGIC get first index for 'window' for intra, gives back
            # 1D array
            gene_indices,
        ]

    elif direction == "Upstream" and te_category == "Order":
        slice_to_return = dd_instance.left_orders[
            dd_instance.order_index_dict[te_name],
            dd_instance.window_index_dict[window_val],
            gene_indices,
        ]

    elif direction == "Downstream" and te_category == "Order":
        slice_to_return = dd_instance.right_orders[
            dd_instance.order_index_dict[te_name],
            dd_instance.window_index_dict[window_val],
            gene_indices,
        ]

    elif direction == "Upstream" and te_category == "Superfamily":
        slice_to_return = dd_instance.left_supers[
            dd_instance.super_index_dict[te_name],
            dd_instance.window_index_dict[window_val],
            gene_indices,
        ]

    elif direction == "Downstream" and te_category == "Superfamily":
        slice_to_return = dd_instance.right_supers[
            dd_instance.super_index_dict[te_name],
            dd_instance.window_index_dict[window_val],
            gene_indices,
        ]
    else:
        raise ValueError()

    return DensitySlice(slice_to_return, direction, window_val, te_name)


def add_hdf5_indices_to_gene_data(
    dd_instance,
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
        dd_instance (TODO)


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
    verify_uniq_chrom_pandaframe(gene_info_pandas, chrom_col)
    verify_chromosome_match_w_pandaframe(dd_instance, gene_info_pandas, chrom_col)

    # NB get around setting with copy warning by creating a deep copy,
    # not an issue for performance.
    gene_info_pandas = gene_info_pandas.copy(deep=True)
    gene_info_pandas[index_col] = gene_info_pandas.apply(
        lambda x: dd_instance._index_of_gene(x[gene_name_col]), axis=1
    )
    return gene_info_pandas


def yield_all_slices(dd_instance):
    # TODO think of refactoring, refactor other code to use the other
    # methods.
    # NOTE I, Scott, hate this method and it was the best I could figure out at
    # the time for interrogating the blueberry data. Only the bluberry example
    # code uses it. I think there are now better methods...
    """
    Yields a DensitySlice (named tuple) object for each TE
    type/direction/window combination for all genes in the DensityData obj.
    """
    directions = ["Upstream", "Downstream"]
    for direction in directions:
        for window_idx, window_val in enumerate(dd_instance.window_list):
            for te_type, te_order_idx in dd_instance.order_index_dict.items():
                if "Revision" in te_type:  # MAGIC dont use the revision
                    # set, which is an artifact of calculating TE density
                    continue
                if direction == "Upstream":
                    yield DensitySlice(
                        dd_instance.left_orders[te_order_idx, window_idx, :],
                        direction,
                        window_val,
                        te_type,
                    )
                if direction == "Downstream":
                    yield DensitySlice(
                        dd_instance.right_orders[te_order_idx, window_idx, :],
                        direction,
                        window_val,
                        te_type,
                    )

            # Iterate over SuperFamilies now
            for te_type, te_super_idx in dd_instance.super_index_dict.items():
                if "Revision" in te_type:  # MAGIC dont use the revision
                    # set, which is an artifact of calculating TE density
                    continue
                if direction == "Upstream":
                    yield DensitySlice(
                        dd_instance.left_supers[te_super_idx, window_idx, :],
                        direction,
                        window_val,
                        te_type,
                    )
                if direction == "Downstream":
                    yield DensitySlice(
                        dd_instance.right_supers[te_super_idx, window_idx, :],
                        direction,
                        window_val,
                        te_type,
                    )


def verify_gene_in_dd_instance(dd_instance, gene_name):
    """
    Check if gene name is actually in the DensityData instance. Raise an error
    if it is not.
    """
    # NOTE this check is semi redundant with index of gene, consider
    # refactoring
    if gene_name not in dd_instance.gene_list:
        raise ValueError(
            f"""Gene {gene_name} is not in your list of genes for
                this specific density data object:
            {dd_instance}"""
        )


# TODO refactor this to use the window value argument NOT the window index
# TODO introduce some sort of check to make sure their gene is actually in the
# dd instance
def info_of_gene(dd_instance, gene_id, window_idx, n_te_types=5):
    """
    gene_id (str): String representing the name of the gene to report
    information

    window_idx (int): Integer representing the index of the window that you
    want to display information for, the smallest window starts at 0.

    n_te_types (int): Defaults to 5, integer representing how many TE types
    to show for the greatest and least values.
    """
    verify_gene_in_dd_instance(dd_instance, gene_id)
    gene_index = dd_instance._index_of_gene(gene_id)
    window_val = dd_instance.window_list[window_idx]

    # TODO candidate to clean up below, utility function for those not
    # well-versed in indexing HDF5

    # TODO this does not avoid the O_Revision and S_Revision datasets,
    # future release will handle these artifacts more elegantly. Unable to
    # avoid at the moment.

    # NB sort in reverse order, the last items in this array contains
    # the greatest density value
    # ORDERS
    sorted_order_left_indices = np.argsort(
        dd_instance.left_orders[:, window_idx, gene_index]
    )
    sorted_order_intragenic_indices = np.argsort(
        dd_instance.intra_orders[:, 0, gene_index]
    )
    sorted_order_right_indices = np.argsort(
        dd_instance.right_orders[:, window_idx, gene_index]
    )

    # SUPERFAMILIES
    sorted_super_left_indices = np.argsort(
        dd_instance.left_supers[:, window_idx, gene_index]
    )
    sorted_super_intragenic_indices = np.argsort(
        dd_instance.intra_supers[:, 0, gene_index]
    )
    sorted_super_right_indices = np.argsort(
        dd_instance.right_supers[:, window_idx, gene_index]
    )

    # NB sorted array, containing actual values now
    # ORDERS
    sorted_order_left = dd_instance.left_orders[:, window_idx, gene_index][
        sorted_order_left_indices
    ]
    sorted_order_intra = dd_instance.intra_orders[:, 0, gene_index][
        sorted_order_intragenic_indices
    ]
    sorted_order_right = dd_instance.right_orders[:, window_idx, gene_index][
        sorted_order_right_indices
    ]
    # SUPERFAMILIES
    sorted_super_left = dd_instance.left_supers[:, window_idx, gene_index][
        sorted_super_left_indices
    ]
    sorted_super_intra = dd_instance.intra_supers[:, 0, gene_index][
        sorted_super_intragenic_indices
    ]
    sorted_super_right = dd_instance.right_supers[:, window_idx, gene_index][
        sorted_super_right_indices
    ]

    info = f"""
    -------------------------------------------------
    Gene Index and ID: {gene_index, gene_id}
    Window Index and Value: {window_idx, window_val}
    Chosen N TE Types to Display: {n_te_types}
    TE Orders in Annotation: {dd_instance.order_list}
    TE SuperFamilies in Annotation: {dd_instance.super_list}
    NOTE, there can be multiple entries with 0 as a density value, so if
    zero are returned, there very well might be more.
    -------------------------------------------------

    TOP {n_te_types} TE ORDERS:
        Upstream:
            {np.array(dd_instance.order_list)[sorted_order_left_indices[-n_te_types:]]}
            {sorted_order_left[-n_te_types:]}
        Intragenic:
            {np.array(dd_instance.order_list)[sorted_order_intragenic_indices[-n_te_types:]]}
            {sorted_order_intra[-n_te_types:]}
        Downstream:
            {np.array(dd_instance.order_list)[sorted_order_right_indices[-n_te_types:]]}
            {sorted_order_right[-n_te_types:]}

    BOTTOM {n_te_types} TE ORDERS:
        Upstream:
            {np.array(dd_instance.order_list)[sorted_order_left_indices[:n_te_types]]}
            {sorted_order_left[:n_te_types]}
        Intragenic:
            {np.array(dd_instance.order_list)[sorted_order_intragenic_indices[:n_te_types]]}
            {sorted_order_intra[:n_te_types]}
        Downstream:
            {np.array(dd_instance.order_list)[sorted_order_right_indices[:n_te_types]]}
            {sorted_order_right[:n_te_types]}

    -------------------------------------------------

    TOP {n_te_types} SUPERFAMILIES:
        Upstream:
            {np.array(dd_instance.super_list)[sorted_super_left_indices[-n_te_types:]]}
            {sorted_super_left[-n_te_types:]}
        Intragenic:
            {np.array(dd_instance.super_list)[sorted_super_intragenic_indices[-n_te_types:]]}
            {sorted_super_intra[-n_te_types:]}
        Downstream:
            {np.array(dd_instance.super_list)[sorted_super_right_indices[-n_te_types:]]}
            {sorted_super_right[-n_te_types:]}

    BOTTOM {n_te_types} SUPERFAMILIES:
        Upstream:
            {np.array(dd_instance.super_list)[sorted_super_left_indices[:n_te_types]]}
            {sorted_super_left[:n_te_types]}
        Intragenic:
            {np.array(dd_instance.super_list)[sorted_super_intragenic_indices[:n_te_types]]}
            {sorted_super_intra[:n_te_types]}
        Downstream:
            {np.array(dd_instance.super_list)[sorted_super_right_indices[:n_te_types]]}
            {sorted_super_right[:n_te_types]}

    -------------------------------------------------
    """
    return info


def verify_uniq_chrom_pandaframe(gene_info_pandas, chrom_col):
    # MAGIC to check chromosome identity
    if len(gene_info_pandas[chrom_col].unique()) > 1:
        raise ValueError(
            """Your gene pandaframe has too many unique
                pseudomolecules in it, it should only have genes
                belonging to one pseudomolecule."""
        )


def verify_chromosome_match_w_pandaframe(dd_instance, gene_info_pandas, chrom_col):
    if dd_instance.unique_chromosome_id != gene_info_pandas[chrom_col].unique()[0]:
        raise ValueError(
            """The pseudomolecule of DensityData instance does
                not match pseudomolecule of the supplied pandas
                object."""
        )
