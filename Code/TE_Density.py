# By Scott Teresi
#-------------------------------------------
# Imports

import numpy as np
import pandas as pd
import os
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
    ch_output_data_path() # change to output directory
    Data.to_csv(Output_Name,
            header=True,
            sep=',')
    ch_main_path() # change to Code directory


#-------------------------------------------
# PandaFrame Helper Functions

def get_dtypes(my_df):
    print(TE_Data.dtypes)

def get_nulls(my_df):
    null_columns = my_df.columns[my_df.isnull().any()]
    count_of_null = my_df[null_columns].isnull().sum()
    print('Counts of null values per column: ' '\n', count_of_null, '\n')

    rows_where_null = my_df[my_df.isnull().any(axis=1)][null_columns].head()
    print('Rows where null exist: ', '\n', rows_where_null, '\n')

def drop_nulls(my_df):
    print('DROPPING ROWS WITH AT LEAST ONE NULL VALUE!!!')
    my_df = my_df.dropna(axis = 0, how ='any')
    return my_df


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

    TE_Data[['Family', 'Subfamily']] = TE_Data.Feature.str.split('/', expand=True)
        # step to fix TE names

    TE_Data = TE_Data.drop('Feature', axis=1)
    TE_Data.Subfamily = TE_Data.Subfamily.astype(str)
        # Some of the subfamilies are None, needs to be string

    TE_Data = drop_nulls(TE_Data) # Dropping nulls because I checked
    TE_Data.Start = TE_Data.Start.astype(int) # Converting to int for space
    TE_Data.Stop = TE_Data.Stop.astype(int) # Converting to int for space






if __name__ == '__main__':
    #Gene_Data = import_genes()
    #get_head(Gene_Data)
    import_transposons()
    



