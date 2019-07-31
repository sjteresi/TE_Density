# By Scott Teresi
#-------------------------------------------
# Imports
import numpy as np
import pandas as pd
import os
import time

from multiprocessing import Process
import multiprocessing
from threading import Thread
import logging
from collections import defaultdict
#-------------------------------------------
# Directory Movement and Helpers 

def cwd():
    """ Prints the current working directory """
    print(os.getcwd())

def ch_input_data_path():
    """ Change the path to the input data folder """
    os.chdir('/home/scott/Documents/Uni/Research/Raw_Strawberry/Camarosa/')

def ch_output_data_path():
    """ Change the path to the output data folder """
    os.chdir('/home/scott/Documents/Uni/Research/Projects/TE_Density/Output_Data/')


def ch_main_path():
    """ Change the path to the code data folder """
    os.chdir('/home/scott/Documents/Uni/Research/Projects/TE_Density/Code/')

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
    ch_output_data_path() # change to output directory
    Data.to_csv(Output_Name,
            header=True,
            sep=',')
    ch_main_path() # change to Code directory


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


#-------------------------------------------
# Main Functions

def import_genes():
    """ Import Genes File """
    my_data = 'camarosa_gtf_data.gtf' # DECLARE YOUR DATA NAME 
    ch_input_data_path()

    col_names = ["Chromosome", 'Software', 'Feature', 'Start', 'Stop', \
        'Score', 'Strand', 'Frame', 'FullName']

    col_to_use = ["Chromosome", 'Software', 'Feature', 'Start', 'Stop', \
        'FullName']

    Gene_Data = pd.read_csv(my_data,
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

    Gene_Data = Gene_Data.drop(['FullName', 'Name1'], axis = 1)
        # remove extraneous rows
 
    Gene_Data.Start = Gene_Data.Start.astype(int) # Converting to int for space
    Gene_Data.Stop = Gene_Data.Stop.astype(int) # Converting to int for space
    Gene_Data['Length'] = Gene_Data.Stop - Gene_Data.Start + 1 # check + 1
    return Gene_Data


def import_transposons():
    """ Import the TEs """
    my_data = 'camarosa_gff_data.gff' # DECLARE YOUR DATA NAME 
    ch_input_data_path()

    col_names = ["Chromosome", 'Software', 'Feature', 'Start', 'Stop', \
        'Score', 'Strand', 'Frame', 'Attribute']

    col_to_use = ["Chromosome", 'Software', 'Feature', 'Start', 'Stop']

    TE_Data = pd.read_csv(my_data,
            sep='\t+',
            header=None,
            engine='python',
            names = col_names,
            usecols = col_to_use)
            #dtype = {'Start':int})

    TE_Data[['Family', 'SubFamily']] = TE_Data.Feature.str.split('/', expand=True)
    TE_Data.SubFamily.fillna(value='Unknown', inplace=True) # replace None w U
        # step to fix TE names


    TE_Data = TE_Data.drop('Feature', axis=1)

    TE_Data = drop_nulls(TE_Data) # Dropping nulls because I checked
    TE_Data.Start = TE_Data.Start.astype(int) # Converting to int for space
    TE_Data.Stop = TE_Data.Stop.astype(int) # Converting to int for space
    TE_Data['Length'] = TE_Data.Stop - TE_Data.Start + 1 # check + 1 
    TE_Data = TE_Data[TE_Data.Family != 'Simple_repeat'] # drop s repeat
    TE_Data = replace_names(TE_Data)
    #get_head(TE_Data)
    #save_output(TE_Data, 'Cleaned_TEs.csv')
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

def density_algorithm(genes, tes, window, increment, max_window):
    try:
        get_unique(genes.Chromosome) == get_unique(tes.Chromosome)
    except:
        raise ValueError("You do not have the same chromosomes in your files")

    #Genes_W_Density = genes.copy(deep=True)

    # Create subsets
    # Creates a filtered list by chromosome to iterate over
    for ind_chromosome in get_unique(genes.Chromosome):
        g_filter = genes['Chromosome'] == ind_chromosome
        subset_genes = genes.where(g_filter)
        subset_genes = drop_nulls(subset_genes)

        t_filter = tes['Chromosome'] == ind_chromosome
        subset_tes = tes.where(t_filter)
        subset_tes = drop_nulls(subset_tes)

        while window <= max_window:
            Genes_W_Density = subset_genes.copy(deep=True)
            Genes_W_Density = init_empty_densities(Genes_W_Density, tes, window)
            #save_output(Genes_W_Density, 'Test_Nan.csv')
            #raise NameError
            for g_ind, g_row in subset_genes.iterrows():
                for t_ind, t_row in subset_tes.iterrows():
                    




    #print(TE_Data.Family.unique())
    #get_counts(TE_Data.Family)
    #save_output(TE_Data, 'Test.csv')
    #Genes_W_Density['LTR_left'] = np.nan                  

            #genes.loc[g_ind, "TE_Density LTR"] = 

            #g_row['Stop'] + window / 

    #print(genes.where(genes['Chromosome'] == 'Fvb1-1'))
            raise NameError
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
    return my_genes


if __name__ == '__main__':
    Gene_Data = import_genes()
    TE_Data = import_transposons()
    density_algorithm(
                    Gene_Data,
                    TE_Data,
                    window=500,
                    increment=500,
                    max_window=10000
                    )



