"""
A collection of helpful dataframe functions.
"""
__author__ = "Scott Teresi"

import pandas as pd

# NOTE this file is a candidate for deletion later, currently these functions
# are not in use but they may be useful at a later date.


def get_nulls(my_df):
    """
    Print out a count of null values per column.
    Print out the row IDs where the null values exist

    Args:
        my_df (Pandaframes): Pandaframe to check null values in
    """
    null_columns = my_df.columns[my_df.isnull().any()]
    count_of_null = my_df[null_columns].isnull().sum()
    print("Counts of null values per column: " "\n", count_of_null, "\n")
    rows_where_null = my_df[my_df.isnull().any(axis=1)][null_columns].head()
    print("Rows where null exist: ", "\n", rows_where_null, "\n")


def drop_nulls(my_df, status=False):
    """
    Drop null values inside a Pandaframe

    Args:
        my_df (Pandaframes): Pandaframe to drop null values
    """
    if status:
        print("DROPPING ROWS WITH AT LEAST ONE NULL VALUE!!!")
    my_df = my_df.dropna(axis=0, how="any")
    return my_df


def swap_columns(dataframe, col_condition, col_1, col_2):
    """
    Swap the values of two designated columns for a row based on a column
    condition, return the dataframe

    Args:
        my_df (Pandaframes): Pandaframe to swap columns in.
    """
    dataframe.loc[col_condition, [col_1, col_2]] = dataframe.loc[
        col_condition, [col_2, col_1]
    ].values
    return dataframe
