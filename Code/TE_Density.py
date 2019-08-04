# By Scott Teresi
#-------------------------------------------
# Imports
import numpy as np import pandas as pd
import os
import time

#from multiprocessing import Process
#import multiprocessing
#from threading import Thread
#import logging
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

def swap_columns(df, col_condition, c1, c2):
    df.loc[col_condition, [c1, c2]] = \
    df.loc[col_condition, [c2, c1]].values
    return df

#-------------------------------------------
# Main Functions

def import_genes():
    """ Import Genes File """
    my_data = 'camarosa_gtf_data.gtf' # DECLARE YOUR DATA NAME 
    ch_input_data_path()

    col_names = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
                 'Score', 'Strand', 'Frame', 'FullName']

    col_to_use = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
                  'Strand', 'FullName' ]

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

    Gene_Data.Strand = Gene_Data.Strand.astype(str)
    Gene_Data.Start = Gene_Data.Start.astype(int) # Converting to int for space
    Gene_Data.Stop = Gene_Data.Stop.astype(int) # Converting to int for space
    Gene_Data['Length'] = Gene_Data.Stop - Gene_Data.Start + 1 # check + 1


    col_condition = Gene_Data['Strand'] == '-'
    Gene_Data = swap_columns(Gene_Data, col_condition, 'Start', 'Stop')
    return Gene_Data


def import_transposons():
    """ Import the TEs """
    my_data = 'camarosa_gff_data.gff' # DECLARE YOUR DATA NAME 
    ch_input_data_path()

    col_names = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
        'Score', 'Strand', 'Frame', 'Attribute']

    col_to_use = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
                 'Strand']

    TE_Data = pd.read_csv(my_data,
            sep='\t+',
            header=None,
            engine='python',
            names = col_names,
            usecols = col_to_use)

    TE_Data[['Family', 'SubFamily']] = TE_Data.Feature.str.split('/', expand=True)
    TE_Data.SubFamily.fillna(value='Unknown', inplace=True) # replace None w U
        # step to fix TE names


    TE_Data = TE_Data.drop('Feature', axis=1)

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

def density_algorithm(genes, tes, window, increment, max_window):
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






            #G['Inside'] = G.apply(TEs_localization, T = T, axis=1)
            #get_head(G)
            #get_head(T)
            #G['Inside'] = G.apply(lambda x: T[(T['Start'] >= x['Start']) & (T['Stop'] <= x['Stop'])]['Start'].count(), axis=1)
            
       
            #G.Inside = G.apply(lambda x: T[(T.Start >= x.Start) & (T.Stop <= x.Stop)]['Start'].count(), axis=1)
            #G.Inside = G.apply(TEs_localization, T=T, axis = 1)

            #print(T[(T.Start >= G.Start) & (T.Stop <=G.Stop)]['Start'].count())
            #print(G.equals(G2))


            A.apply(lambda x: B[(B['Start'] >= x['Start']) & (B['Stop'] <= x['Stop'])].apply(lambda y : (y['Start'] / y['Stop']) +1), axis=1)
            #G['Inside'] = G.apply(TEs_localization, TEs=T, axis= 1)
            #print(type(G) = (G.apply(lambda x: T[(T['Start'] >= x['Start']) & (T['Stop'] <= x['Stop'])], axis = 1)))
                                  #.apply(lambda y : (y['Start'] / y['Stop']) +1), axis=1)

            #G.Inside = G.apply(lambda x: T[(T.Start >= x['Start']) & (T['Stop'] <= x['Stop'])]['Start'].count(), axis=1)
            #G.Inside = G.apply(lambda x: T[(T.Start >= x['Start']) & (T['Stop'] <= x['Stop'])]['Start'].mul(1000), axis=1)



            #G.Inside = G.apply(lambda x: T[T(['Start'] >= x['Start']) & (T['Stop'] <= x['Stop'])]x['Start']+x['Stop'])
            #G.Inside = G.apply(lambda x: T[(T['Start'] >= x['Start']) & (T['Stop'] <= x['Stop'])]x['Start'] + x['Stop'], axis = 1)
            #print(T[T(['Start'] >= G['Start']) & (T['Stop'] <= G['Stop'])]['Start'].count())


            #for g_ind, g_row in G.iterrows():
                #g_row.apply(lambda g_row: TEs_localization(g_row, T=T), axis=1)
                #for t_ind, t_row in T.iterrows():
                    #g_row['TEs_inside'] += 1


                #print(g_row)
                #g_row.TEs_inside[T.Start < T.Stop] += 1
            #subset_genes.TEs_inside = subset_genes.mask(T['Start'] < 1500, subset_genes.TEs_inside += 1)


            get_head(G)
            save_output(G, 'Test_Inside.csv')
            raise ValueError


            #T.mask(subset_genes.Start < T.Start, 24)
            #get_head(T)

            #for g_ind, g_row in subset_genes.iterrows():
                #subset_genes.where([(g_row.Start <= T.Start) & (T.Stop <= g_row.Stop), 'TEs_inside'] += 1
            #print(subset_genes.apply(lambda row: row["Start"] + row["Stop"], axis = 1))
            #print(subset_genes.apply(add_them, df2=subset_tes, axis=1))
            #raise NameError
            #Genes_W_Density.loc[(subset_tes.Start.values >=
                                 #Genes_W_Density.Start.values) &
                                #(subset_tes.Stop.values <=
                                 #Genes_W_Density.Stop.values), 'TEs_inside'] += 1
            #Genes_W_Density = Genes_W_Density.apply(TEs_inside, T = subset_tes, axis = 1)

            #subset_genes.loc[(T.Start >= subset_genes.Start) & (T.Stop <= subset_genes.Stop), 'TEs_inside'] += 1
            #subset_genes = subset_genes.reset_index(drop=True)
            #T = T.reset_index(drop=True)
            #subset_genes.loc[(subset_genes.Start <= T.Start) & (T.Stop <= subset_genes.Stop), 'TEs_inside'] += 1
            #get_head(subset_genes)
            #save_output(subset_genes, 'Test_Inside.csv')
            #print(subset_genes.apply(TEs_inside, T=subset_tes, axis=1))





            window += increment

def TEs_localization(G,TEs):
    pass
    # input is x
    # TEs are y
    #adjusted_G = TEs.apply(lambda y: y['Inside'] if (G['Stop'] >= y['Start'])else None, axis =1)
    #adjusted_G = adjusted_G.dropna(axis=0, how='all')
    #if adjusted_G.empty:
        #adjusted_G = np.nan
    #return adjusted_G

   
    #return (T[(T.Start >= G.Start) & (T.Stop <= G.Stop)]['Stop'].div(1000))
    #return (T[(T.Start >= G.Start) & (T.Stop <= G.Stop)]['Stop'].div(1000))
    #return (T[(T.Stop <= G.Stop)](['Stop'].div(10000) + 1))
    #raise ValueError
    #return (T[(T.Start >= G.Start) & (T.Stop <= G.Stop)]['Start'].count())

        #count += 1

    #G.loc(axis=1)[(T.Start >= G.Start) & (T.Stop <= G.Stop), 'TEs_inside'] += 1
    #G.where((T.Start >= G.Start) & (T.Stop <= G.Stop), I
    #G.TEs_inside = I

    #(T.Start >= G.Start) & (T.Stop <= G.Stop) = G.TEs_inside += 1
        #G.TEs_inside += 1
    #return G
    #df.where(m, -df) == np.where(m, df, -df)
    #if T.Start >= G.Start and T.Stop <= G.Stop:
        #return T.Start + T.Start
    
    #return df.Start + df.Stop


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
    Gene_Data = import_genes()
    TE_Data = import_transposons()
    density_algorithm(
                    Gene_Data,
                    TE_Data,
                    window=10000,
                    increment=500,
                    max_window=10000
                    )

