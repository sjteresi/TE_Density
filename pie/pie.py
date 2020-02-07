#!/usr/bin/env python3

"""
Generate some nice pie graphs for showing TE distributions

TODO
    - Get pie graph of TE distros

"""

__author__ = "Scott Teresi"

import os
import argparse
import pandas as pd
import time
import matplotlib.pyplot as plt
import numpy as np
import coloredlogs
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

from transposon.replace_names import TE_Renamer
from transposon.import_genes import import_genes
from transposon.import_transposons import import_transposons

#-------------------------------------------
# PandaFrame Helper Functions

def get_dtypes(my_df):
    print(my_df.dtypes)

def get_nulls(my_df):
    null_columns = my_df.columns[my_df.isnull().any()]
    count_of_null = my_df[null_columns].isnull().sum()
    print('Counts of null values per column: ' '\n', count_of_null, '\n')
    rows_where_null = my_df[my_df.isnull().any(axis=1)][null_columns].head()
    print('Rows where null exist: ', '\n', rows_where_null, '\n')

def drop_nulls(my_df, status=False):
    if status:
        print('DROPPING ROWS WITH AT LEAST ONE NULL VALUE!!!')
    my_df = my_df.dropna(axis = 0, how ='any')
    return my_df

def get_info(my_df):
    print(my_df.info())

def get_counts(my_df_col):
    print(my_df_col.value_counts())

def get_unique(my_df_col):
    return my_df_col.unique()

def swap_columns(df, col_condition, c1, c2):
    df.loc[col_condition, [c1, c2]] = \
    df.loc[col_condition, [c2, c1]].values
    return df

def split(df, group):
    """
    Returns a list of the dataframe with each element being a subset of the df.
    I use this function to split by chromosome so that we may later do
    chromosome element-wise operations.
    """
    gb = df.groupby(group)
    return [gb.get_group(x) for x in gb.groups]


#-------------------------------------------

def gene_names(sub_gene_data):
    """Return unique gene names for the input gene data (e.g. one chromosome)."""

    # MAGIC_NUMBER the gene name column is 'Gene_Name'
    gene_name_id = 'Gene_Name'
    names = sub_gene_data[gene_name_id].unique()
    return names


def sum_genes_and_TEs(Gene_Data, TE_Data):
    # The Genome Size is: 0.805489 gigabasepair
    #genome_size = 0.805489 # Cam
    genome_size = 0.24
    my_lengths = []
    gene_lengths = Gene_Data.Length.sum()
    TE_lengths = TE_Data.Length.sum()


    # Convert to Gbp
    gene_lengths = gene_lengths / 1000000000 / genome_size # normalize
    TE_lengths = TE_lengths / 1000000000 / genome_size # normalize
    other = 1 - gene_lengths - TE_lengths
    my_lengths.append(other)
    my_lengths.append(gene_lengths)
    my_lengths.append(TE_lengths)

    #print(f"Total Lengths of Genes: {gene_lengths} Gbp")
    #print(f"Total Lengths of TEs: {TE_lengths} Gbp")
    return my_lengths

def main_fig(TE_Data):
    family_labels = get_unique(TE_Data.Family) # list
    subfamily_labels = get_unique(TE_Data.SubFamily) # list

    print(family_labels)
    print(subfamily_labels)


    sizes = sum_genes_and_TEs(Gene_Data, TE_Data)
    total_num_TEs = len(TE_Data.index)
    #print(f"Total Rows of TEs: {total_num_TEs}")

    #--------------------------------------------#
    # Plotting
    def genome_content():
        # Genome Content Plot
        labels = ['Other', 'Genes', 'TEs']
        new_sizes = []
        for size in sizes:
            new_sizes.append(size * 100)
        colors = ['forestgreen', 'cornflowerblue', 'crimson']
        explode = (0, 0, 0.1) # explode 3rd slice
        patches, texts, autotext = plt.pie(new_sizes, colors=colors, shadow=True,
                                    startangle=120, explode=explode,
                                    autopct='%1.1f%%', labels=labels)
        plt.legend(patches, labels, loc='best')
        plt.title('Contribution to Genome Size')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig('H4_Genome_Content.png')
        plt.close()
        #plt.show()

    def te_fam_content():
        # TE Distro Plot
        my_family_counts = TE_Data.Family.value_counts()
        my_family_counts = my_family_counts.to_dict()
        family_sum = 0
        for key,val in my_family_counts.items():
            family_sum += val
        new_sizes = []
        for key,val in my_family_counts.items():
            my_family_counts[key] = val/family_sum*100 # normalize
            new_sizes.append(my_family_counts[key])

        labels = my_family_counts.keys()
        colors = ['seagreen', 'deepskyblue', 'lightcoral', 'darkviolet']
        explode = (0.1, 0, 0, 0, 0) # explode
        patches, texts, autotext = plt.pie(new_sizes, colors=colors, shadow=True,
                                    startangle=150, explode=explode,
                                    autopct='%1.1f%%', labels=labels)
        plt.legend(patches, labels, loc='best')
        plt.title('Percentage of Number of TEs Compared to Total Number of TEs (Family)')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig('H4_Family.png')
        plt.close()
        #plt.show()


    def te_subfam_content():
        # TE Distro Plot
        subfamily_labels = get_unique(TE_Data.SubFamily) # list
        my_subfamily_counts = TE_Data.SubFamily.value_counts()
        my_subfamily_counts = my_subfamily_counts.to_dict()
        subfamily_sum = 0


        my_subfamily_counts['Pure Unknown'] = my_subfamily_counts.pop('Unknown_SubFam')


        for key,val in my_subfamily_counts.items():
            subfamily_sum += val
        new_sizes = []
        for key,val in my_subfamily_counts.items():
            my_subfamily_counts[key] = val/subfamily_sum*100 # normalize
            new_sizes.append(my_subfamily_counts[key])

        labels = my_subfamily_counts.keys()
        colors = ['sienna', 'goldenrod', 'lightgreen', 'forestgreen',
                  'steelblue', 'darkslateblue', 'orchid', 'crimson', 'mistyrose' ]
        explode = (0, 0.1, 0.1, 0, 0, 0, 0, 0, 0) # explode
        patches, texts, autotext = plt.pie(new_sizes, colors=colors, shadow=True,
                                    startangle=210, explode=explode,
                                    autopct='%1.1f%%', labels=labels)
        #plt.legend(patches, labels, loc='best')
        plt.title('Percentage of Number of TEs Compared to Total Number of TEs (Subfamily)')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig('H4_Subfamily.png')
        plt.close()
        #plt.show()


    def graphs_chromosomes():
        # Loop
        #n = 2.9
        #fig.subplots_adjust(hspace = n, wspace = n)
        #n = 4
        colors = ['sienna', 'goldenrod', 'lightgreen', 'forestgreen',
                  'steelblue', 'darkslateblue', 'orchid', 'crimson', 'mistyrose' ]

        explode = (0, 0.1, 0.1, 0, 0, 0, 0, 0, 0) # explode

        all_c = get_unique(TE_Data.Chromosome) # list

        for xi in range(0, len(all_c), 4):
            fig, axes = plt.subplots(nrows = 2, ncols = 2)
            #print(f"XI is: {xi}")
            #raise ValueError

            for i, ax in enumerate(axes.ravel()):
                #print(f"i is: {i}")
                list_subset = all_c[xi: xi + i + 1]
                titles = get_unique(grouped_TEs[xi + i].Chromosome) # list
                #print(f"Title is: {titles} \n")
                size_vals =  []
                my_subfamily_counts = grouped_TEs[xi+i].SubFamily.value_counts()
                my_subfamily_counts = my_subfamily_counts.to_dict()
                subfamily_sum = 0
                my_subfamily_counts['Pure Unknown'] = my_subfamily_counts.pop('Unknown')

                for key,val in my_subfamily_counts.items():
                    subfamily_sum += val
                for key,val in my_subfamily_counts.items():
                    my_subfamily_counts[key] = val/subfamily_sum*100 # normalize
                    size_vals.append(my_subfamily_counts[key])

                labels = my_subfamily_counts.keys()
                ax.pie(size_vals,
                       explode = explode,
                       labels = None,
                       autopct = '%1.1f%%',
                       shadow = True,
                       startangle = 90,
                       radius = 1.7,
                       pctdistance=0.84,
                       textprops={'fontsize': 6.8})
                ax.set_title(titles[0], loc='left', pad=-10.0)
                plt.tight_layout()
            #plt.tight_layout(pad=n, w_pad = n, h_pad = n)
            fig.legend(labels=labels,
                       loc = 'lower right',
                       title='TEs',
                      borderaxespad=0.1)
            n = 0.9
            plt.subplots_adjust(right=0.75, hspace = n)
            plt.savefig(all_c[xi]+'_Cam_All.png')
            #plt.show()














    genome_content()
    te_fam_content()
    te_subfam_content()
    #graphs_chromosomes()

def check_groupings(grouped_genes,grouped_TEs):
    """
    Function to make sure the chromosome pairs of the genes and the TEs are
    in correct. Checks the first 10 lines of each pair. This is just to make
    sure that each pair of chromosomes are right. Correct subsetting would be
    managed by my custom split command, which has been tested.

    Args:
        grouped_genes (list of pandaframes): Gene dataframes separated by chromosome
        grouped_TEs (list of pandaframes): TE dataframes separated by chromosome
    """
    #print(grouped_TEs)
    #print(grouped_genes)
    zipped_set = zip(grouped_genes, grouped_TEs)
    try:
        for g_element, t_element in zipped_set:
                #print(type(g_element))
                assert g_element.Chromosome.iloc[0:10].values[0] == \
                t_element.Chromosome.iloc[0:10].values[0]
    except AssertionError as error:
        raise ValueError('Chromosomes do not match for the grouped_genes or grouped_TEs')

def validate_args(args, logger):
    """Raise if an input argument is invalid."""

    if not os.path.isdir(args.input_dir):
        logger.critical("argument 'input_dir' is not a directory")
        raise ValueError("%s is not a directory"%(abs_path))
    if not os.path.isdir(args.output_dir):
        logger.critical("argument 'output_dir' is not a directory")
        raise ValueError("%s is not a directory"%(abs_path))


if __name__ == '__main__':
    """Command line interface to calculate density."""

    parser = argparse.ArgumentParser(description="calculate TE density")
    path_main = os.path.abspath(__file__)
    parser.add_argument('input_dir', type=str,
                        help='parent directory of gene & transposon files')
    parser.add_argument('--output_dir', '-o', type=str,
                        default=os.path.join(path_main, '../..', 'results'),
                        help='parent directory to output results')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='set debugging level to DEBUG')
    args = parser.parse_args()
    args.output_dir = os.path.abspath(args.output_dir)
    args.input_dir = os.path.abspath(args.input_dir)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    logger.info("Start processing directory '%s'"%(args.input_dir))
    for argname, argval in vars(args).items():
        logger.debug("%-12s: %s"%(argname, argval))
    validate_args(args, logger)

    # FUTURE move this preprocessing to it's object
    logger.info("Importing transposons, this may take a moment...")
    TE_Data = import_transposons(args.input_dir, 'H4_TEs.gff')
    print(TE_Data.head())
