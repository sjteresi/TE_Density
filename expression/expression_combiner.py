#!/usr/bin/env python3

"""
This file takes the input gtf file of all the genes, and combines it with our
expression data to create a unified table for data statistics
"""

__author__ = "Scott Teresi"
import argparse
import os
import pandas as pd
import numpy as np
#-----------------------------------------------
def import_genes():
    os.chdir('/home/scott/Documents/Uni/Research/Raw_Strawberry/Camarosa/')
    gtf_filename = 'camarosa_gtf_data.gtf' # DECLARE YOUR DATA NAME

    col_names = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
                 'Score', 'Strand', 'Frame', 'FullName']

    col_to_use = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
                  'Strand', 'FullName' ]

    Gene_Data = pd.read_csv(
            gtf_filename,
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

    Gene_Data.Strand = Gene_Data.Strand.astype(str)
    Gene_Data.Start = Gene_Data.Start.astype(int) # Converting to int for space
    Gene_Data.Stop = Gene_Data.Stop.astype(int) # Converting to int for space
    Gene_Data['Length'] = Gene_Data.Stop - Gene_Data.Start + 1 # check + 1
    #col_condition = Gene_Data['Strand'] == '-' # do not swap, missing function
    #Gene_Data = swap_columns(Gene_Data, col_condition, 'Start', 'Stop')

    # NOTE Gene_Data now complete
    return Gene_Data

def import_exp_counts():
    os.chdir('/home/scott/Documents/Uni/Research/Raw_Strawberry/Camarosa/Expression')
    my_files = os.listdir()
    num = 0
    df_list = []
    for my_file in my_files:
        if my_file.endswith('.count'):
            num += 1
            count_name = my_file.split('.')[0]
            exp_data_subset = pd.read_csv(
                                    my_file,
                                    sep='\t+',
                                    header=None,
                                    engine='python',
                                    names=['Gene_Name', f"{count_name}_Fruit"],
                                    dtype={'Gene_Name':str,f"{count_name}_Fruit":int}
                                    )
            df_list.append(exp_data_subset)
    df_list = (df.set_index('Gene_Name') for df in df_list) # generator object
    df = next(df_list)
    for df_ in df_list:
        df = df.merge(df_, on='Gene_Name')
    return df

def merge_all(genes,expression):
    All_Data = pd.merge(left=genes,right=expression, how='left', \
                             left_on='Gene_Name', right_on='Gene_Name')
    return All_Data

if __name__ == '__main__':
    parser = False # NOTE change this line as desired for input data
    if parser:
        parser = argparse.ArgumentParser(description='CSV to PD_Dataframe')
        parser.add_argument('indir', type=str, help='Input dir for data')
        parser.add_argument('--outdir', type=str, help='Output dir for data')
        args = parser.parse_args()
        print(args.indir)
    else:
        pass

    Gene_Data = import_genes()
    expression_full = import_exp_counts()
    All_Data = merge_all(Gene_Data,expression_full)
    All_Data.to_csv('Camarosa_GenesWExpression.csv')

    #-----------------------------
