#!/usr/bin/env python3

"""
This file takes the input gtf file of all the genes, and combines it with our
expression data to create a unified table for data statistics
"""

__author__ = "Scott Teresi"
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

    # NOTE
    # The below criteria language is used to appropriately get around the:
    # SettingWithCopy Warning

    criteria = Gene_Data['Feature'] == 'exon'
    criteria_row_indices = Gene_Data[criteria].index
    Gene_Exons = Gene_Data.loc[criteria_row_indices, :]
    Gene_Exons['Gene_Name'] = Gene_Exons.FullName.str.split('-mRNA').str.get(0).str.split('ID=').str.get(1)
    Gene_Exons.drop('FullName', axis = 1, inplace = True)

    del Gene_Data

    Gene_Exons = Gene_Exons.set_index('Gene_Name')
    Gene_Exons.Start = Gene_Exons.Start.astype(int) # Converting to int for space
    Gene_Exons.Stop = Gene_Exons.Stop.astype(int) # Converting to int for space
    Gene_Exons['Exon_Length'] = Gene_Exons.Stop - Gene_Exons.Start + 1 # NOTE check
    Gene_Exons = Gene_Exons.groupby(Gene_Exons.index)['Exon_Length'].sum()
    Gene_Exons.to_csv('Gene_Exons.csv', header=['Exon_Length'])

    return Gene_Exons

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
    df_merged = pd.merge(left=genes,right=expression, how='left', \
                             left_on='Gene_Name', right_on='Gene_Name')
    return df_merged


def TPM(Genes_W_Exp):
    """
    Transcripts Per Million
    First make sure gene length is in KB

    Next normalize by gene length, divide each read count by gene length, that
    is RPK

    Following that, normalize for sequencing depth. Add up the read counts,
    already normalized for gene length, and get the total for each replicate
    (tissue).

    Then divide that total by 1000000. That gives us the Scaling Factor.

    Finally divide the read counts (that are already divided by the gene length) by the
    Scaling Factor
    """
    Genes_W_Exp['Exon_Length'] = Genes_W_Exp['Exon_Length'] / 1000
    # put exon length in kb

    f = lambda x : x / Genes_W_Exp['Exon_Length'] if x.name != 'Exon_Length' else x
    Genes_W_Exp = Genes_W_Exp.apply(f)
    # normalize by gene length 

    f = lambda x : x / (x.sum() / 1000000) if x.name != 'Exon_Length' else x
    Genes_W_Exp = Genes_W_Exp.apply(f)
    # divide the read counts by the scaling factor

    return Genes_W_Exp

if __name__ == '__main__':
    Gene_Data = import_genes()
    expression_full = import_exp_counts()
    Genes_Exp = merge_all(Gene_Data,expression_full)
    Genes_Exp_Full = TPM(Genes_Exp)
    Genes_Exp_Full.to_csv('All_Genes_TPM.csv')
    #-----------------------------
