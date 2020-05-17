#!/usr/bin/env python3

__author__ = "Scott Teresi"

"""
Output the relevant area (divisor) for density calculation.
"""


from transposon.density import validate_window
import pandas as pd


class NormMatrix(object):
    """
    Functions for calculating the relevant area (divisor) necessary for the
    normalization matrix.
    """

    @staticmethod
    def divisor_left(GeneData, windows):
        """
        Get the relevant area for the left side. Returns a numpy array of
        windows (rows) and genes (columns). When converting from pandas to
        numpy it removes the indices (which are the windows in my pandaframe)
        and removes the headers (which are the column name gene names in my
        pandaframe. The genes (columns) are in the order supplied from
        GeneData.names, which to my knowledge, is the exact same order as they
        appear in the annotation file.

        Left: Window length most of the time, in the edge cases where the window
        extends past 0, turning negative, we need to clip the value to 0 and make
        sure the relevant area is between 0 and the left window stop.

        Args:
            GeneData (GeneData): The wrapped gene annotation file.
            windows (iterable of ints): number of base pairs given by the windowing
            operation, used to calculate the relevant area (divisor).
        """
        # NOTE there might be a way we can have a function generate the values
        # on the creation of the dataframe below, but it would require
        # suppliying both the window and the genedatum to a separate function.
        # Look at intra, but it would be more complex because we need the
        # window

        # TODO update all of the args for each divisor function.
        # TODO examine performance
        genes_window_dataframe = pd.DataFrame(
                                            columns=[name for name in GeneData.names],
                                            index=windows)
        for window in windows:
            for name in GeneData.names:
                gene_datum = GeneData.get_gene(name)
                win_length = gene_datum.win_length(window)
                win_start = gene_datum.left_win_start(window)
                win_stop = gene_datum.left_win_stop
                win_length = validate_window(win_start, win_stop, win_length)
                genes_window_dataframe.at[window, name] = win_length
        return genes_window_dataframe.to_numpy(copy=False)

    @staticmethod
    def divisor_intra(GeneData, windows):
        """
        Get the relevant area for the intronic area. Returns a numpy array of
        windows (rows) and genes (columns). When converting from pandas to
        numpy it removes the indices (which are the windows in my pandaframe)
        and removes the headers (which are the column name gene names in my
        pandaframe. The genes (columns) are in the order supplied from
        GeneData.names, which to my knowledge, is the exact same order as they
        appear in the annotation file.

        Intra: Gene length, here we are considering intronic TEs, so the relevant
        area is the gene length.

        Args:
        """
        # NOTE
        # Technically the value repeats, as it makes a list for one row, and
        # then re-uses that list for the second row, BUT that is okay, because
        # it is always the same value (the gene length) and the divisor value
        # (the gene length) is independent of the window value for the intra
        # valculations
        genes_window_dataframe = pd.DataFrame(
                                             [[GeneData.get_gene(name).length for name in
                                             GeneData.names]],
                                             columns=[name for name in GeneData.names],
                                             index=windows)
        return genes_window_dataframe.to_numpy(copy=False)

    @staticmethod
    def divisor_right(GeneData, windows):
        """
        Get the relevant area for the right area. Returns a numpy array of
        windows (rows) and genes (columns). When converting from pandas to
        numpy it removes the indices (which are the windows in my pandaframe)
        and removes the headers (which are the column name gene names in my
        pandaframe. The genes (columns) are in the order supplied from
        GeneData.names, which to my knowledge, is the exact same order as they
        appear in the annotation file.

        Right: Window length, here there are no edge cases as the positions TEs or
        genes in the genome is bounded from [0, inf).

        Args:
        """
        genes_window_dataframe = pd.DataFrame(
                                            columns=[name for name in GeneData.names],
                                            index=windows)
        for window in windows:
            for name in GeneData.names:
                gene_datum = GeneData.get_gene(name)
                win_length = gene_datum.win_length(window)
                genes_window_dataframe.at[window, name] = win_length
        return genes_window_dataframe.to_numpy(copy=False)

    def __repr__(self):
        """
        String representation for developer.
        """
        info = """
               GeneData: {self._gene_data}
               Window: {self._window}
               """
        return info.format(self=self)
