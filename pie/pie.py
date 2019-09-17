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
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


#-------------------------------------------
# Directory Movement and Helpers

def cwd():
    """ Prints the current working directory """
    print(os.getcwd())


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

    Gene_Data = Gene_Data[Gene_Data.Feature == 'gene']  # drop non-gene rows

    # clean the names and set as the index (get row wrt name c.f. idx)
    Gene_Data[['Name1', 'Gene_Name']] = Gene_Data.FullName.str.split(';Name=', expand=True)
    Gene_Data.set_index('Gene_Name', inplace=True)
    Gene_Data = Gene_Data.drop(['FullName', 'Name1', 'Software'], axis = 1)

    Gene_Data.Strand = Gene_Data.Strand.astype(str)
    Gene_Data.Start = Gene_Data.Start.astype('uint32')
    Gene_Data.Stop = Gene_Data.Stop.astype('uint32')
    # NOTE is the gene index a closed | open interval?
    # Scott thinks we should delete the + 1 for length
    Gene_Data['Length'] = Gene_Data.Stop - Gene_Data.Start

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
    # NOTE see same comment on gene data types
    TE_Data.Start = TE_Data.Start.astype('uint32')
    TE_Data.Stop = TE_Data.Stop.astype('uint32')
    # NOTE see same comment on gene intervals / off-by-one
    TE_Data['Length'] = TE_Data.Stop - TE_Data.Start
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

def gene_names(sub_gene_data):
    """Return unique gene names for the input gene data (e.g. one chromosome)."""

    # MAGIC_NUMBER the gene name column is 'Gene_Name'
    gene_name_id = 'Gene_Name'
    names = sub_gene_data[gene_name_id].unique()
    return names


def sum_genes_and_TEs(Gene_Data, TE_Data):
    # The Genome Size is: 0.805489 gigabasepair
    genome_size = 0.805489
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
        plt.savefig('Genome_Content.png')
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
        explode = (0.1, 0, 0, 0) # explode
        patches, texts, autotext = plt.pie(new_sizes, colors=colors, shadow=True,
                                    startangle=150, explode=explode,
                                    autopct='%1.1f%%', labels=labels)
        plt.legend(patches, labels, loc='best')
        plt.title('Percentage of Number of TEs Compared to Total Number of TEs (Family)')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig('Family.png')
        plt.close()
        #plt.show()


    def te_subfam_content():
        # TE Distro Plot
        subfamily_labels = get_unique(TE_Data.SubFamily) # list
        my_subfamily_counts = TE_Data.SubFamily.value_counts()
        my_subfamily_counts = my_subfamily_counts.to_dict()
        subfamily_sum = 0


        my_subfamily_counts['Pure Unknown'] = my_subfamily_counts.pop('Unknown')


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
        plt.savefig('Subfamily.png')
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
            plt.savefig(all_c[xi]+'_All.png')
            #plt.show()














    genome_content()
    te_fam_content()
    te_subfam_content()
    graphs_chromosomes()

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
    #print("\ngene data...")
    #get_head(Gene_Data)
    # NOTE grouped_genes is a list of all the data frames
    grouped_genes = split(Gene_Data, 'Chromosome') # check docstring for my split func


    TE_Data = import_transposons()
    #print("\nTE data...")
    #get_head(TE_Data)
    grouped_TEs = split(TE_Data, 'Chromosome') # check docstring for my split func


    #check_groupings(grouped_genes, grouped_TEs)

    main_fig(TE_Data)




