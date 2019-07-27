# By Scott Teresi
#-------------------------------------------
# Imports
import numpy as np
import pandas as pd
import os
#-------------------------------------------

# Helper functions

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
    try:
        print(Data.head())
    except AttributeError as e:
        raise AttributeError('Your input data was incorrect, check to make sure it is a Pandaframe')

def save_output(Data, Output_Name):
    ch_output_data_path() # change to output directory
    Data.to_csv(Output_Name,
            header=True,
            sep=',')
    ch_main_path() # change to Code directory



def import_data():
    my_data = 'camarosa_gtf_data.gtf'
    ch_input_data_path()

    col_names = ["Chromosome",
            'Software',
            'Feature',
            'Start',
            'Stop',
            'Nonsense1',
            'Nonsense2',
            'Nonsense3',
            'FullName']

    col_to_use = ["Chromosome",
            'Software',
            'Feature',
            'Start',
            'Stop',
            'FullName']

    Gene_Data = pd.read_csv(my_data,
            sep='\t+',
            header=None,
            engine='python',
            names = col_names,
            usecols = col_to_use,)


    Gene_Data = Gene_Data[Gene_Data.Feature == 'gene'] # drop non-gene rows
    Gene_Data[['Name1', 'Gene_Name']] = Gene_Data.FullName.str.split(';Name=', expand=True) # first step to fix names
    Gene_Data = Gene_Data.set_index('Gene_Name') # set new names as index
    Gene_Data = Gene_Data.drop(['FullName', 'Name1'], axis = 1) # remove extraneous rows




    get_head(Gene_Data)










    save_output(Gene_Data, 'test.csv')















if __name__ == '__main__':
    import_data()



    #print(list(Gene_Data.columns))

    #get_head(Gene_Data)
