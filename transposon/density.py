#!/usr/bin/env python3

"""
Calculate transposable element density.

"""

__author__ = "Scott Teresi, Michael Teresi"

import os
import time
import argparse
import logging
from enum import IntEnum, unique
#from multiprocessing import Process
#import multiprocessing
#from threading import Thread

import numpy as np
import pandas as pd

# FUTURE enum.Flag is more appropriate but let us delineate first and revisit
@unique
class DensityConditions(IntEnum):
    """Enumerate the conditions where density is calculated with respect to.

    TODO scott, pls add a one liner explanation for each test

    Below, inclusive generally can be evaluated as <= or >=
    Exclusive would be < or >

    Conditions deal with the transposable element (TE) position wrt gene/windows.


        IN_GENE_ONLY: The TE is only inside the gene, that includes gene start
        and gene stop (inclusive). The density would be the total number of TE
        bases (TE length) divided by the Gene Length.


        IN_WINDOW_ONLY: The TE is between the window start or stop (inclusive)
        and the gene start or stop (exclusive). The density would be the total
        number of TE bases (TE length) divided by the size of the window.


        IN_WINDOW_AND_GENE: The TE is within the gene (inclusive on gene edges)
        and ends within the window (can end on window edges, inclusive). This
        test is sort of a fusion of IN_GENE_ONLY and IN_WINDOW_ONLY, we would
        need to calculate two values!
        Something that we would call Intra density and either upstream or downstream
        density. The intra value would just come from the portion that is inside,
        divided by the length of the gene. And the window portion would just be
        the portion that is inside the window. Reminder that a TE base that
        overlaps with the gene start or stop is considered inside the gene,
        where a TE base that ends on a window is considered inside the window

        IN_WINDOW_FROM_OUTSIDE: A TE extends into the window, but part of it is
        outside the window. Thus the only relevant part is the part that sits
        inside the window. So imagine a case where the TE extends past the
        WindowStop, the desired amount of TE Bases would be:
        (WindowStop - TE_Start) and then we would divide that by the WindowSize
        to get the density



    """

    IN_GENE_ONLY= 0
    IN_WINDOW_ONLY = 1
    IN_WINDOW_AND_GENE = 2
    IN_WINDOW_NOT_GENE = 3
    IN_WINDOW_AND_GENE_UP_DOWN_STREAM = 4  # TE in window or gene but start or stop is outside window?


#-------------------------------------------
# Directory Movement and Helpers

def cwd():
    """ Prints the current working directory """
    print(os.getcwd())

# def ch_main_path():  # TODO remove, not sure where this is needed anymore
#     """ Change the path to the code data folder """
#     os.chdir('/home/scott/Documents/Uni/Research/Projects/TE_Density/Code/')

def get_head(Data):
    """ Get the heading of the Pandaframe """
    try:
        print(Data.head())
    except AttributeError as e:
        raise AttributeError('Your input data was incorrect, check to make sure it is a Pandaframe')

def save_output(Data, Output_Name):
    """ Save the output of the Pandaframe """
    if Output_Name[-4:] != '.csv':
        raise NameError('Please make sure your filename has .csv in it!')
    global OUTPUT_DIR  # TODO remove global
    Data.to_csv(os.path.join(OUTPUT_DIR, Output_Name),
            header=True,
            sep=',')
    # TODO remove ch_main_path, not sure where this is needed anymore
    #ch_main_path() # change to Code directory

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

#-------------------------------------------
# Main Functions

def import_genes():
    """ Import Genes File """

    # TODO remove MAGIC NUMBER (perhaps just search by extension (gtf)?)
    gtf_filename = 'camarosa_gtf_data.gtf' # DECLARE YOUR DATA NAME
    #ch_input_data_path()

    col_names = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
                 'Score', 'Strand', 'Frame', 'FullName']

    col_to_use = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
                  'Strand', 'FullName' ]

    global INPUT_DIR  # TODO remove global
    Gene_Data = pd.read_csv(
            os.path.join(INPUT_DIR, gtf_filename),
            sep='\t+',
            header=None,
            engine='python',
            names = col_names,
            usecols = col_to_use)

    Gene_Data = Gene_Data[Gene_Data.Feature == 'gene']
        # drop non-gene rows

    Gene_Data[['Name1', 'Gene_Name']] = Gene_Data.FullName.str.split(';Name=', expand=True)
        # first step to fix names

    Gene_Data = Gene_Data.set_index('Gene_Name')
        # set new names as index

    Gene_Data = Gene_Data.drop(['FullName', 'Name1', 'Software'], axis = 1)
        # remove extraneous rows

    Gene_Data.Strand = Gene_Data.Strand.astype(str)
    Gene_Data.Start = Gene_Data.Start.astype(int) # Converting to int for space
    Gene_Data.Stop = Gene_Data.Stop.astype(int) # Converting to int for space
    # MAGIC NUMBER NOTE
    Gene_Data['Length'] = Gene_Data.Stop - Gene_Data.Start + 1 # check + 1
    # Scott thinks we should delete the + 1 for length

    # We will not swap Start and Stop for Antisense strands. We will do this
    # post-processing
    #col_condition = Gene_Data['Strand'] == '-'
    #Gene_Data = swap_columns(Gene_Data, col_condition, 'Start', 'Stop')
    return Gene_Data


def import_transposons():
    """ Import the TEs """
    # TODO remove MAGIC NUMBER (perhaps just search by extension (gtf)?)
    gff_filename = 'camarosa_gff_data.gff' # DECLARE YOUR DATA NAME
    #ch_input_data_path()

    col_names = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
        'Score', 'Strand', 'Frame', 'Attribute']

    col_to_use = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
                 'Strand']

    global INPUT_DIR  # TODO remove global
    TE_Data = pd.read_csv(
            os.path.join(INPUT_DIR, gff_filename),
            sep='\t+',
            header=None,
            engine='python',
            names = col_names,
            usecols = col_to_use)

    TE_Data[['Family', 'SubFamily']] = TE_Data.Feature.str.split('/', expand=True)
    TE_Data.SubFamily.fillna(value='Unknown', inplace=True) # replace None w U
        # step to fix TE names

    TE_Data = TE_Data.drop(['Feature', 'Software'], axis=1)

    TE_Data = drop_nulls(TE_Data) # Dropping nulls because I checked
    TE_Data.Strand = TE_Data.Strand.astype(str)
    TE_Data.Start = TE_Data.Start.astype(int) # Converting to int for space
    TE_Data.Stop = TE_Data.Stop.astype(int) # Converting to int for space
    TE_Data['Length'] = TE_Data.Stop - TE_Data.Start + 1 # check + 1
    TE_Data = TE_Data[TE_Data.Family != 'Simple_repeat'] # drop s repeat
    TE_Data = replace_names(TE_Data)
    return TE_Data

def replace_names(my_TEs):
    U = 'Unknown'
    master_family = {
        'RC?':'DNA',
        'RC':'DNA',
        'SINE?':U
    }

    U = 'Unknown_SubFam'
    master_subfamily = {
        'Uknown':U,
        'MuDr':'MULE',
        'MULE-MuDR':'MULE',
        'Pao':U,
        'Caulimovirus':U,
        'hAT-Tag1':'hAT',
        'hAT-Tip100':'hAT',
        'hAT-Charlie':'hAT',
        'Helitron':U,
        'unknown':U,
        'Maverick':U,
        'Harbinger':'PIF-Harbinger',
        'TcMar-Pogo':U,
        'CR1':'LINE',
        'hAT-Ac':'hAT',
        'L2':'LINE',
        'L1':'LINE',
        'Jockey':'LINE',
        'MuLE-MuDR':'MULE',
        'MuDR':'MULE',
        'Mutator':'MULE',
        'Micro_like':U
    }

    my_TEs.Family.replace(master_family, inplace=True)
    my_TEs.SubFamily.replace(master_subfamily, inplace=True)
    return my_TEs


class DensityInputs(object):
    """Contains data for calculating transposon density."""

    # NOTE for now, just store a simplified np.array
    # FUTURE store data frame and add properties for the different views
    def __init__(self, genes, transposons, window):
        """Initializer.

        Args:
            genes (numpy.array): Jx2, J # of genes, column1 start idx, column2 stop index
            transposons (numpy.array): Kx2, K # of transposons, column1 start idx, column2 stop idx
        """
        self.genes = gene
        self.transposons = transposons
        self.window = window



class DensityOutput(object):
    """Contains data on a density result."""
    # NOTE to scott from mike: use this to do a map/reduce strategy
    # basically, let the caller deal with formatting the result
    # don't collate them in the same function that calculates them

    def __init__(self, density, condition):
        self.density = density
        self.condition = condition

# NOTE we should probably just make the conditions a class b/c we are going to have
# tons of copy paste code otherwise
# also, although it is functional (data in, data out)
# one can consider the window part of the state that the functions share
# maybe contain the genes, transposons, write the conditions / rho as clasmethods
# then provide helpers to glue it together (e.g. loop over conditions w/ same interface)

def is_inside_only(genes, transposon):
    """Return where the TE is only inside the gene.

    Args:
        genes (numpy.array): Jx2, J # of genes, column1 start idx, column2 stop index
        transposons (numpy.array): array of 2, start idx, stop idx
    """

    # NOTE just use a loop for now above this for multiple TE
    # use np.newaxis in the future to broadcast
    # or can we load a matrix that size in ram?

    t_start = transposon[0]
    t_stop = transposon[1]

    #g_start = genes[0]
    #g_stop = genes[1]

    down = t_start >= genes[:,0]
    up = t_stop <= genes[:,1]

    return np.logical_and(up, down)


def in_left_window_only(genes, transposon, window):
    """ Return where the TE is only inside the window (no overlap with gene)
        On the left side

    Args:
        genes (numpy.array): Jx2, J # of genes, column1 start idx, column2 stop index
        transposons (numpy.array): array of 2, start idx, stop idx
        window (int): integer value of the current window
    """
    t_start = transposon[0]
    t_stop = transposon[1]

    #np.where(genes - window < 0, genes, 0)

    # NOTE no check for window < 0 being set to 0
    down = t_start >= genes[:,0] - window # test still passes
    up = t_stop < genes[:,0] # not inclusive

    return np.logical_and(up, down)


def rho(genes, transposons, passed_condition, window_start, window_stop):
    """Calculate the density where it passed the conditional test."""
    pass

    # TODO need to account for in window and intra density
    # not necessarily in this function
    # anything that overlaps the gene it has an intra density
    # anything that overlaps the window it has a window density
    # some can have both of these

    # no. bases / no. relevant bases

    # window density
    # no. bases of TE / window  (when TE is between gene / window)

    # intra (the relevant ares is the gene)
    # no. bases of TE / gene (when TE is inside gene)

def line_overlap(min_a, max_a, min_b, max_b):
    pass

def rho_intra(genes_start, genes_stop, genes_length, transposon_start, transposon_stop):
    """Intra density of a transposon wrt genes.

    Args:

    """
    assert genes_start.shape == genes_stop.shape
    assert genes_start.shape == genes_length.shape

    lower = np.minimum(genes_stop, transposon_stop)
    upper = np.maximum(genes_start, transposon_start)
    # MAGIC NUMBER stop is inclusive so add one, scott thinks no
    # NOTE we are using 0 indexing I believe. 
    te_overlaps =  np.maximum(0, lower - upper)
    densities = np.divide(
        te_overlaps,
        genes_length,
        out=np.zeros_like(te_overlaps, dtype='float'),
        where=genes_length!=0
    )

    assert densities.shape == genes_start.shape
    return densities

def rho_left_window(genes_start, genes_stop, window, transposon_start, transposon_stop, transposon_length):

    """Density to the left (downstream) of a gene.
    When TE is between gene and window

    Relevant things: No. TE bases / window
    GeneStart so that we may calculate the WindowStart
    WindowStart is unique to each gene. Calculated via Gstart - Window
    WindowStart is the leftmost window


    Args:
        window (int): integer value of the current window
    """
    # TODO make sure the edge cases work when the window should be negative and
    # reset to 0 for that instance
    # so far I can reset the window start for that instance, but I cannot fix
    # the window for the sake of division later when you try to check values.

    assert genes_start.shape == genes_stop.shape

    window_start = np.subtract(genes_start, window)
    window_start[window_start < 0] = 0 # this might be an effective way to force 0 

    lower_bound = np.maximum(window_start, transposon_start)
    upper_bound = np.minimum(genes_start, transposon_stop)
    te_overlaps =  np.maximum(0, upper_bound - lower_bound)
    densities = np.divide(
        te_overlaps,
        window,
        out=np.zeros_like(te_overlaps, dtype='float')
    )

    assert densities.shape == genes_start.shape
    return densities

def rho_right_window(genes_start, genes_stop, window, transposon_start, transposon_stop, transposon_length):

    """Density to the right (upstream) of a gene.
    When TE is between gene and window

    Relevant things: No. TE bases / window
    GeneStop so that we may calculate the WindowStop
    WindowStop is unique to each gene. Calculated via Gstop + Window
    WindowStop is the rightmost window

    Args:

    """
    assert genes_start.shape == genes_stop.shape

    window_stop = np.add(genes_stop, window)
    lower_bound = np.maximum(genes_stop, transposon_start)
    # lower bound gets TE starts to the right of gene stops
    upper_bound = np.minimum(window_stop, transposon_stop)
    # upper bound gets TE stops to the left of window stops
    te_overlaps =  np.maximum(0, upper_bound - lower_bound)
    densities = np.divide(
        te_overlaps,
        window,
        out=np.zeros_like(te_overlaps, dtype='float')
    )

    assert densities.shape == genes_start.shape
    return densities

def density_algorithm(genes, tes, window, increment, max_window):
    """
    te data frame has columns: SEE import_transposons
    te data frame has rows: each row is a temp


    """
    # NOTE create 2 structs to hold files / param & result

    # MICHAEL SHOULD LOOK AT THIS COMMAND
    try:
        get_unique(genes.Chromosome) == get_unique(tes.Chromosome)
    except:
        raise ValueError("You do not have the same chromosomes in your files")

    #Genes_W_Density = genes.copy(deep=True)

    # Create subsets
    # Creates a filtered list by chromosome to iterate over
    for ind_chromosome in get_unique(genes.Chromosome):
        g_filter = genes['Chromosome'] == ind_chromosome
        G = genes.where(g_filter)
        G = drop_nulls(G)

        t_filter = tes['Chromosome'] == ind_chromosome
        T = tes.where(t_filter)
        T = drop_nulls(T)

        # At this point I have a subset of genes to iterate over
        # This subset is on a chromosome by chromosome basis
        # I will perform my windowing operations on this subset
        # Perhaps in the future I can add multiprocessing for each chromosome

        #Genes_W_Density = subset_genes.copy(deep=True)
        #Genes_W_Density = init_empty_densities(Genes_W_Density, tes, window)
        #G['Inside'] = np.nan
        while window <= max_window:
            logging.debug(" Gene df shape:  {}".format(genes.values.shape))
            logging.debug(" TE df shape:  {}".format(tes.values.shape))
            # Perform the windowing operations
            # Multiple tests need to be done here for each window
            # All the tests must be run for that window and then rerun for the
            # next window
            # I foresee that for each gene file we will have an output wiht all
            # of the original gene data and the output will be 500_LTR_Upstream
            # The TE types (LTR) are given by the TE_Dataframe.Family and
            # TE_Dataframe.SubFamily columns

            # The init_empty_densities function was to add the appropriate
            # columns, we may not need to worry about that for now


            # All the commented code below are my attempts to do the work
            #-----------------------------


            get_head(G)
            save_output(G, 'Test_Inside.csv')
            raise ValueError


            window += increment


def init_empty_densities(my_genes, my_tes, window):
    Family_List = my_tes.Family.unique()
    SubFamily_List = my_tes.SubFamily.unique()
    Directions = ['_upstream', '_intra', '_downstream']

    for family in Family_List:
        for direction in Directions:
            col_name = (str(window) + '_' + family + direction)
            my_genes[col_name] = np.nan

    for subfamily in SubFamily_List:
        for direction in Directions:
            col_name = (str(window) + '_' + family + direction)
            my_genes[col_name] = np.nan
    my_genes['TEs_inside'] = np.nan
    return my_genes


if __name__ == '__main__':

    logging.basicConfig(level=os.environ.get("LOGLEVEL", "DEBUG"))
    parser = argparse.ArgumentParser(description="calculate TE density")
    parser.add_argument('input_dir', type=str,
                        help='parent directory of gene & transposon files')
    parser.add_argument('--output_dir', '-o', type=str,
                        default=os.path.abspath(__file__) + 'Output_Data',
                        help='parent directory to output results')
    args = parser.parse_args()

    global OUTPUT_DIR  # TODO remove global
    OUTPUT_DIR = args.output_dir
    global INPUT_DIR  # TODO remove global
    INPUT_DIR = args.input_dir

    Gene_Data = import_genes()
    get_head(Gene_Data)

    TE_Data = import_transposons()
    get_head(TE_Data)


    density_algorithm(
                    Gene_Data,
                    TE_Data,
                    window=1000,
                    increment=500,
                    max_window=10000
                    )
