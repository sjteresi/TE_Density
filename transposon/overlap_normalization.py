#!/usr/bin/env python3

"""
Create a matrix of the number of occurrences of overlaps for each base pair in
the genome. Tally the number of occurrences of the specific TE on the base pair,
record the identity of the TE.
"""

__author__ = "Scott Teresi"

import numpy as np
import pandas as pd

from transposon.transposon_data import TransposonData
from transposon.gene_data import GeneData

from statistics import mode

class O_NormMatrix(object):
    """
    Test
    """

    def __init__(self, wrapped_gene_data, wrapped_transposon_data):
        """
        Initialize.

        Args:
            wrapped_transposon_data (TransposonData): A TransposonData object
        """

        self.gs = wrapped_gene_data
        self.g_dframe = wrapped_gene_data.data_frame[['Chromosome',
                                                    'Start', 'Stop',
                                                    ]].copy(deep=True).sort_values(by=['Start'])
        self.ts = wrapped_transposon_data
        self.t_dframe = wrapped_transposon_data.data_frame[['Chromosome',
                                                            'Start', 'Stop',
                                                            'Order',
                                                            'SuperFamily']].copy(deep=True).sort_values(by=['Start'])
        self.calc_range_col()
        #self.t_dframe.Order.astype(str)
        #self.max_stop = int(self.t_dframe.Stop.max())
        #self.len_t_dframe = len(wrapped_transposon_data.data_frame)

    #def get_range_by_order(self):

        #orig = np.empty((self.max_stop,))
        #for i in range(0, self.len_t_dframe, 1):
            #to_cat = np.arange(int(self.ts.starts[i]), int(self.ts.stops[i]), 1)
            #orig = np.concatenate((orig, to_cat)).astype(int)
        #print(orig.dtype)
        #bincount = np.bincount(orig)
        #print(orig)
        #print(orig.shape)
        #print(self.max_stop)

    def calc_range_col(self):
        self.t_dframe['Range'] = [list(range(int(i), int(j+1))) for i, j in self.t_dframe[['Start','Stop']].values]

    def panda_method(self):
        # NOTE Unused
        for item in self.ts.order_name_set:
            self.t_dframe.loc[self.t_dframe['Order'] == item, item] = 1
            self.t_dframe.loc[self.t_dframe['Order'] != item, item] = 0

        for item in self.ts.superfamily_name_set:
            self.t_dframe.loc[self.t_dframe['SuperFamily'] == item, item] = 1
            self.t_dframe.loc[self.t_dframe['SuperFamily'] != item, item] = 0

    def matrix_method(self):
        unique_bp = set([a for b in self.t_dframe.Range.tolist() for a in b])
        unique_bp = sorted(unique_bp)

        superfam_dict = {}
        for superfam in self.ts.superfamily_name_set:
            superfam_dict[superfam] = [inner for outer
                                       in self.t_dframe.Range[self.t_dframe['SuperFamily']
                                                             ==
                                                             superfam].values.tolist()
                                       for inner in outer]

        colnames = self.ts.order_name_set.union(self.ts.superfamily_name_set)
        x = pd.DataFrame(index=unique_bp, columns=colnames)
        x.index.rename('Unique_BP', inplace=True)
        for bp_to_match in unique_bp:
            for superfam, bp_list in superfam_dict.items():
                x.loc[[bp_to_match], superfam] = bp_list.count(bp_to_match)


        x.to_csv('test.csv')


    def __repr__(self):
        """
        String representation for developer.
        """
        #info = "Transposon Data Pandaframe({self.transposons.data_frame})"
        info = "Count Matrix {self.count_matrix.shape}"
        return info.format(self=self)

