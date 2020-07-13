#!/usr/bin/env python3

"""
Generate some descriptive graphs for the strawberry genome TE data.
"""

__author__ = "Scott Teresi"

import os
import argparse
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
import coloredlogs
import logging
from functools import reduce
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)



matplotlib.style.use('ggplot')


from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.import_genes import import_genes
from transposon.import_transposons import import_transposons
from transposon.replace_names import te_annot_renamer
from transposon.density import verify_TE_cache, verify_gene_cache
from transposon.genome_data import GenomeData, SubgenomeData


def get_unique(my_df_col):
    return my_df_col.unique()


def split(df, group):
    """
    Returns a list of the dataframe with each element being a subset of the df.
    I use this function to split by chromosome so that we may later do
    chromosome element-wise operations.
    """
    gb = df.groupby(group)
    return [gb.get_group(x) for x in gb.groups]


def graph_general_genome_stats_pie(GenomeData, output_dir, explode_bool=False):
    """
    Plot the Overall TE, Gene, and Other content in the genome (basepairs)

    Args:
        GenomeData (GenomeData instance): Contains the accessor methods for the
        transposon and gene data.

        output_dir (str): A command-line argument of where to output the
        resultant pie graph.

        explode_bool (boolean): A boolean of whether or not to 'explode' one of
        the pie graph segments, purely visual. Defaults to explode the TE
        section of the pie graph.
    """
    gen_size = GenomeData.genome_size
    G_lengths_MB_fraction = GenomeData.total_G_lengths_MB / gen_size
    T_lengths_MB_fraction = GenomeData.total_T_lengths_MB / gen_size
    other_lengths_MB_fraction = GenomeData.total_other_lengths_MB / gen_size
    sizes = [G_lengths_MB_fraction, T_lengths_MB_fraction,
             other_lengths_MB_fraction]
    sizes_as_percents = [size*100 for size in sizes]
    labels = ['Genes', 'TEs', 'Other']
    colors = ['forestgreen', 'crimson', 'cornflowerblue']
    explode = [0 for item in sizes_as_percents]
    if explode_bool:
        explode[1] = 0.1
    explode = tuple(explode)
    patches, texts, autotext = plt.pie(sizes, colors=colors, shadow=False,
                                       startangle=120, explode=explode,
                                       autopct='%1.1f%%', labels=labels)
    plt.legend(patches, labels, loc='best')
    plt.title('Contribution to Genome Size')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output_dir + '/' + GenomeData.genome_id + \
                '_General_Genome_Content.png')
    plt.close()


def graphs_diploid_chromosomes(TE_Data, output_dir, selection):
    # NOTE unused
    # Loop
    # n = 2.9
    # fig.subplots_adjust(hspace = n, wspace = n)
    # n = 4
    colors = ['sienna', 'goldenrod', 'lightgreen', 'forestgreen',
              'steelblue', 'darkslateblue', 'orchid', 'crimson', 'mistyrose',
              'darkmagenta']

    # explode = (0, 0.1, 0.1, 0, 0, 0, 0, 0, 0)  # explode

    grouped_TEs = split(TE_Data, 'Chromosome')
    for chromosome_chunk in grouped_TEs:
        # size_vals = []
        my_chromosome = str(chromosome_chunk.Chromosome.unique()[0])
        my_subfamily_counts = chromosome_chunk.SuperFamily.value_counts()
        my_subfamily_counts = my_subfamily_counts.to_dict()
        subfamily_sum = 0

        for key, val in my_subfamily_counts.items():
            subfamily_sum += val
        new_sizes = []
        for key, val in my_subfamily_counts.items():
            my_subfamily_counts[key] = val/subfamily_sum*100  # normalize
            new_sizes.append(my_subfamily_counts[key])

        labels = my_subfamily_counts.keys()
        colors = ['saddlebrown', 'goldenrod', 'lightgreen', 'forestgreen',
                  'steelblue', 'darkslateblue', 'orchid', 'crimson', 'mistyrose',
                  'gold', 'darkmagenta']
        explode = []
        for i in new_sizes:
            explode.append(0)
        # explode[0] = 0.1
        explode = tuple(explode)
        patches, texts, autotext = plt.pie(new_sizes, colors=colors, shadow=False,
                                           startangle=210, explode=explode,
                                           autopct='%1.1f%%', labels=labels)
        # plt.legend(patches, labels, loc='best')
        plt.title('Relative %s of TEs Compared to Total Number of TEs (SuperFamily)')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig(output_dir + '/' + selection +
                    '_TE_' + my_chromosome + '_SuperFam_Content.png')
        plt.close()


def graph_order_numbers_as_percents(GenomeData, output_dir, explode_bool=False):
    """
    Creates a pie chart displaying TE and gene numbers for
    the given GenomeData object. Displays the  number of a given TE order as a
    percentage of the total number of TEs. Does not consider the length of the
    TEs!

    Args:
        genome (GenomeData): GenomeData object, used to create a dictionary
        of the subgenomes. This can be constructed with genome.Cam_subgenomes()

        output_dir (str): A string of the output directory. Supplied as a
        commandline argument to the main code.

        explode_bool (boolean): A boolean of whether or not to 'explode' one of
        the pie graph segments, purely visual. Defaults to explode the LTR
        section of the pie graph.
    """
    orders = list(sorted(GenomeData.transposon_Order_number_dictionary.keys()))
    values = [GenomeData.transposon_Order_number_dictionary[order] for order in
              orders]
    normalized_numbers_as_percents = [(val /
                                      GenomeData.num_of_TEs)*100 for
                                      val in values]
    colors = ['seagreen', 'royalblue', 'khaki', 'crimson', 'peru']
    explode = [0 for item in normalized_numbers_as_percents]
    if explode_bool:
        explode[3] = 0.1
    explode = tuple(explode)
    patches, texts, autotext = plt.pie(normalized_numbers_as_percents,
                                       colors=colors, shadow=False,
                                       startangle=150, explode=explode,
                                       autopct='%1.1f%%', labels=orders)
    # plt.legend(patches, labels, loc='best')
    plt.title('# of TE Orders Relative to Total # of TEs (Order)')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output_dir + '/' + GenomeData.genome_id +
                '_WG_TE_Order_NumberPercent.png')
    plt.close()


def graph_order_sequences_as_percents(GenomeData, output_dir,
                                      explode_bool=False):
    """
    Creates a pie chart displaying TE sequences as percents of the entire
    TE_nome.

    Args:
        genome (GenomeData): GenomeData object, used to create a dictionary
        of the subgenomes. This can be constructed with genome.Cam_subgenomes()

        output_dir (str): A string of the output directory. Supplied as a
        commandline argument to the main code.

        explode_bool (boolean): A boolean of whether or not to 'explode' one of
        the pie graph segments, purely visual. Defaults to explode the LTR
        section of the pie graph.
    """

    orders = list(sorted(GenomeData.order_sum_seq_len_dict_MB.keys()))
    values = [GenomeData.order_sum_seq_len_dict_MB[order] for order in
              orders]
    normalized_sequences_as_percents = [(val /
                                        GenomeData.total_T_lengths_MB)*100 for
                                        val in values]
    colors = ['seagreen', 'royalblue', 'khaki', 'crimson', 'peru']
    explode = [0 for item in normalized_sequences_as_percents]
    if explode_bool:
        explode[3] = 0.1
    explode = tuple(explode)
    patches, texts, autotext = plt.pie(normalized_sequences_as_percents,
                                       colors=colors, shadow=False,
                                       startangle=150, explode=explode,
                                       autopct='%1.1f%%', labels=orders)
    # plt.legend(patches, labels, loc='best')
    plt.title('Basepairs of TE Orders Relative to Total Basepairs of TEs')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output_dir + '/' + GenomeData.genome_id +
                '_WG_TE_Order_SequencePercent.png')
    plt.close()


def graph_superfam_sequences_as_percents(GenomeData, output_dir,
                                         explode_bool=False):
    """
    Creates a pie chart displaying TE sequences as percents of the entire
    TE_nome.

    Args:
        genome (GenomeData): GenomeData object, used to create a dictionary
        of the subgenomes. This can be constructed with genome.Cam_subgenomes()

        output_dir (str): A string of the output directory. Supplied as a
        commandline argument to the main code.

        explode_bool (boolean): A boolean of whether or not to 'explode' one of
        the pie graph segments, purely visual. Defaults to explode the LTR
        section of the pie graph.
    """

    superfamilies = list(sorted(GenomeData.superfam_sum_seq_len_dict_MB.keys(),
                         key=str.lower))
    values = [GenomeData.superfam_sum_seq_len_dict_MB[superfam] for superfam in
              superfamilies]
    normalized_sequences_as_percents = [(val /
                                        GenomeData.total_T_lengths_MB)*100 for
                                        val in values]
    colors = ['mediumblue', 'dimgrey', 'darkorange', 'crimson', 'forestgreen',
              'darkmagenta', 'saddlebrown', 'peachpuff', 'cornflowerblue',
              'khaki', 'darkturquoise', 'maroon']
    explode = [0 for item in normalized_sequences_as_percents]
    if explode_bool:
        explode[3] = 0.1
    explode = tuple(explode)
    patches, texts, autotext = plt.pie(normalized_sequences_as_percents,
                                       colors=colors, shadow=False,
                                       startangle=150, explode=explode,
                                       autopct='%1.1f%%', labels=superfamilies)
    # plt.legend(patches, labels, loc='best')
    plt.title('Basepairs of TE Superfamilies Relative to Total BP of TEs')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output_dir + '/' + GenomeData.genome_id +
                '_WG_TE_Superfam_SequencePercent.png')
    plt.close()


def graph_superfam_numbers_as_percents(GenomeData, output_dir, explode_bool=False):
    """
    Creates a pie chart displaying TE and gene numbers for
    the given GenomeData object. Displays the  number of a given TE order as a
    percentage of the total number of TEs. Does not consider the length of the
    TEs!

    Args:
        genome (GenomeData): GenomeData object, used to create a dictionary
        of the subgenomes. This can be constructed with genome.Cam_subgenomes()

        output_dir (str): A string of the output directory. Supplied as a
        commandline argument to the main code.

        explode_bool (boolean): A boolean of whether or not to 'explode' one of
        the pie graph segments, purely visual. Defaults to explode the LTR
        section of the pie graph.
    """
    superfamilies = list(sorted(GenomeData.transposon_SuperFam_number_dictionary.keys(),
                         key=str.lower))
    values = [GenomeData.transposon_SuperFam_number_dictionary[superfam] for superfam in
              superfamilies]
    normalized_numbers_as_percents = [(val /
                                      GenomeData.num_of_TEs)*100 for
                                      val in values]
    colors = ['mediumblue', 'dimgrey', 'darkorange', 'crimson', 'forestgreen',
              'darkmagenta', 'saddlebrown', 'peachpuff', 'cornflowerblue',
              'khaki', 'darkturquoise', 'maroon']

    explode = [0 for item in normalized_numbers_as_percents]
    if explode_bool:
        explode[3] = 0.1
    explode = tuple(explode)
    patches, texts, autotext = plt.pie(normalized_numbers_as_percents,
                                       colors=colors, shadow=False,
                                       startangle=150, explode=explode,
                                       autopct='%1.1f%%', labels=superfamilies)
    # plt.legend(patches, labels, loc='best')
    plt.title('# of TE Superfamilies Relative to Total # of TEs')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output_dir + '/' + GenomeData.genome_id +
                '_WG_TE_Superfam_NumberPercent.png')
    plt.close()


def graph_bar_superfam_lengths_genome_wide(list_of_genome_data, file_name, output_dir):
    """
    Creates a stacked bar chart displaying TE and gene cumulative lengths for
    each subgenome.

    Args:
        genome (GenomeData): GenomeData object, used to create a dictionary
        of the subgenomes. This can be constructed with genome.Cam_subgenomes()

        file_name (str): A string of what you want to name the output graph.

        output_dir (str): A string of the output directory. Supplied as a
        commandline argument to the main code.
    """

    fig, ax = plt.subplots()
    colors = ['g', 'mediumblue', 'dimgrey', 'darkorange', 'crimson', 'forestgreen',
              'darkmagenta', 'saddlebrown', 'peachpuff', 'cornflowerblue',
              'khaki', 'darkturquoise', 'maroon']
    ax.set_prop_cycle(color=colors)
    for GenomeData in list_of_genome_data:
        i = 0
        main_label = [GenomeData.subgenome_identity]
        gene_lengths = [GenomeData.total_G_lengths_MB]
        width = 0.15
        ax.bar(main_label, gene_lengths, width, label='Genes', color=colors[i])
        margin_bottom = gene_lengths.pop()
        i += 1

        for key, val in sorted(GenomeData.superfam_sum_seq_len_dict_MB.items()):
            ax.bar(main_label, val, width, bottom=margin_bottom,
                   label=key, color=colors[i])
            i += 1
            margin_bottom += val

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='lower center')
    plt.ylabel('MB')
    plt.xlabel('Subgenome')
    plt.title('TE compositions of the four Camarosa subgenomes')
    plt.savefig(output_dir + '/' + file_name + '.png')
    plt.close()


def graph_cam_superfam_chromosomes(Camarosa, output_dir,
                                   F_vesca, F_iinumae,
                                   F_viridis, F_nipponica):
    # NOTE could use some heavy cleanup when I have the time

    gene_chrom_list = Camarosa.split(Camarosa.gene_dataframe, 'Chromosome')
    te_chrom_list = Camarosa.split(Camarosa.transposon_dataframe, 'Chromosome')

    all_chromosomes = []
    for genes, tes in zip(gene_chrom_list, te_chrom_list):
        all_chromosomes.append(GenomeData('Camarosa', genes, tes))

    # Loop
    # n = 2.9
    # fig.subplots_adjust(hspace = n, wspace = n)
    # n = 4
    colors = ['mediumblue', 'dimgrey', 'darkorange', 'crimson', 'forestgreen',
              'darkmagenta', 'saddlebrown', 'peachpuff', 'cornflowerblue',
              'khaki', 'darkturquoise', 'maroon']

    # explode = (0, 0.1, 0.1, 0, 0, 0, 0, 0, 0)  # explode

    for xi in range(0, len(all_chromosomes), 4):
        fig, axes = plt.subplots(nrows=2, ncols=2)
        # print(f"XI is: {xi}")

        for i, ax in enumerate(axes.ravel()):
            titles = all_chromosomes[xi+i].Chromosomes[0]
            # print(f"Title is: {titles} \n")

            if titles in F_vesca:
                all_chromosomes[xi+i].genome_id = 'F_vesca'
            if titles in F_iinumae:
                all_chromosomes[xi+i].genome_id = 'F_iinumae'
            if titles in F_viridis:
                all_chromosomes[xi+i].genome_id = 'F_viridis'
            if titles in F_nipponica:
                all_chromosomes[xi+i].genome_id = 'F_nipponica'

            superfamilies = list(sorted(
                all_chromosomes[xi+i].transposon_SuperFam_number_dictionary.keys(),
                key=str.lower))
            values = [all_chromosomes[xi+i].transposon_SuperFam_number_dictionary[superfam] for superfam in superfamilies]
            normalized_numbers_as_percents = [(val / all_chromosomes[xi+i].num_of_TEs)*100 for val in values]

            labels = superfamilies
            # NOTE labels is unused here
            ax.pie(values,
                   labels=None,
                   autopct='%1.1f%%',
                   shadow=False,
                   startangle=90,
                   radius=1.7,
                   pctdistance=0.84,
                   colors=colors,
                   textprops={'fontsize': 6.8})
            ax.set_title((all_chromosomes[xi+i].genome_id + ' (' + str(titles)
                         + ')'),
                         loc='left', pad=19)
            plt.tight_layout()
        # n = 5
        # plt.tight_layout(pad=n, w_pad = n, h_pad = n)
        fig.legend(labels=superfamilies,
                   loc='center right',
                   title='TEs',
                   borderaxespad=-0.2)
        n = 0.9
        plt.subplots_adjust(right=0.75, hspace=n)
        plt.savefig(output_dir + '/' + 'Camarosa' +
                    '_' + str(titles) + '_All_Chromosomes.png')
        plt.close()


def graph_average_LTR_length_H4_Cam(H4, Camarosa, output_dir):
    """
    Generate a bar graph of the average length of LTR retrotransposons. Y axis
    is the counts and x-axis is binned sizes in KB. Overlay the H4 and Camarosa
    genomes on top of one another in different translucent colors.
    """
    # print(Camarosa.transposon_dataframe.groupby(['Order']).Length.mean()['LTR'])
    # x = pd.cut(Camarosa.transposon_dataframe.groupby(['Order']).Length'LTR'], bins=5)
    # print(H4.transposon_dataframe.groupby(['Order']).Length)
    ltr_cam = Camarosa.order_transposon_subset('LTR')
    ltr_h4 = H4.order_transposon_subset('LTR')
    # print(ltr_cam.loc[ltr_cam['Length']==48384])
    ltr_cam['Length'] = ltr_cam['Length'].div(1000)
    ltr_h4['Length'] = ltr_h4['Length'].div(1000)
    bin_sizes = np.arange(0, 10, 0.25)

    # What is the ratio between the bin sizes values in Cam and H4
    ltr_cam_value_counts = pd.cut(ltr_cam.Length,
                                  bins=bin_sizes).value_counts()
    ltr_h4_value_counts = pd.cut(ltr_h4.Length,
                                 bins=bin_sizes).value_counts()

    # print(ltr_cam_value_counts[0] / ltr_h4_value_counts[0])

    plt.hist(ltr_cam.Length, bins=bin_sizes, alpha=0.5,
             histtype='bar', ec='black')
    plt.hist(ltr_h4.Length, bins=bin_sizes, alpha=0.5,
             histtype='bar', ec='black')
    plt.ylabel('Count')
    plt.legend(['Camarosa', 'H4'])
    plt.xlabel('Length of LTR retrotransposons in KB')
    plt.savefig(output_dir + '/' + Camarosa.genome_id +
                '_' + H4.genome_id + '_LTR_Bins.png')
    plt.close()


def graph_table_superfam_lengths_genome_wide(Camarosa, list_of_genome_data,
                                             output_dir):
    """

    Args:
        genome (GenomeData): GenomeData object, used to create a dictionary
        of the subgenomes. This can be constructed with genome.Cam_subgenomes()

        file_name (str): A string of what you want to name the output graph.

        output_dir (str): A string of the output directory. Supplied as a
        commandline argument to the main code.
    """

    fig, ax = plt.subplots()
    ax.axis('off')
    ax.axis('tight')

    frames = []
    for GenomeData in list_of_genome_data:
        superfam_and_genes = GenomeData.superfam_sum_seq_len_dict_MB
        superfam_and_genes['Genes'] = GenomeData.total_G_lengths_MB
        subset_df = \
        pd.DataFrame.from_dict(superfam_and_genes,
                               orient='index',
                               columns=[GenomeData.subgenome_identity])
        subset_df.index.name = 'SuperFamily'
        subset_df = subset_df / GenomeData.genome_size * 100
        frames.append(subset_df)

    cam_superfam_and_genes = Camarosa.superfam_sum_seq_len_dict_MB
    cam_superfam_and_genes['Genes'] = Camarosa.total_G_lengths_MB
    cam_superfam_and_genes['Other'] = Camarosa.total_G_lengths_MB
    Cam_data = pd.DataFrame.from_dict(cam_superfam_and_genes,
                                      orient='index',
                                      columns=[Camarosa.genome_id])
    Cam_data.index.name = 'SuperFamily'
    Cam_data = Cam_data / Camarosa.genome_size * 100
    frames.append(Cam_data)

    my_table = reduce(lambda x, y: pd.merge(x, y, on='SuperFamily'), frames)
    my_table = my_table.round(2)
    ax.table(cellText=my_table.values,
             colLabels=my_table.columns,
             loc='center', rowLabels=my_table.index)
    plt.title('Percent of Subgenome Size')
    fig.tight_layout()
    plt.savefig(output_dir + '/' + Camarosa.genome_id +
                '_' + 'Subgenome_Table.png')
    plt.close()


def graph_table_mean_TE_lengths(Camarosa, list_of_genome_data,
                                selection,
                                output_dir):
    # TODO Edit the arg list and try to reduce the number of ifs.
    """
    Args:
    """

    fig, ax = plt.subplots()
    ax.axis('off')
    ax.axis('tight')

    frames = []
    for GenomeData in list_of_genome_data:
        if selection == 'superfam':
            mean_length_dict = GenomeData.average_superfam_length.to_dict()
        if selection == 'order':
            mean_length_dict = GenomeData.average_order_length.to_dict()
        subset_df = \
        pd.DataFrame.from_dict(mean_length_dict,
                               orient='index',
                               columns=[GenomeData.subgenome_identity])
        if selection == 'superfam':
            subset_df.index.name = 'SuperFamily'
        if selection == 'order':
            subset_df.index.name = 'Order'
        frames.append(subset_df)

    if selection == 'superfam':
        mean_length_dict = Camarosa.average_superfam_length.to_dict()
    if selection == 'order':
        mean_length_dict = Camarosa.average_order_length.to_dict()
    Cam_data = pd.DataFrame.from_dict(mean_length_dict,
                                      orient='index',
                                      columns=[Camarosa.genome_id])
    if selection == 'superfam':
        Cam_data.index.name = 'SuperFamily'
    if selection == 'order':
        Cam_data.index.name = 'Order'
    frames.append(Cam_data)

    if selection == 'superfam':
        my_table = reduce(lambda x, y: pd.merge(x, y, on='SuperFamily'), frames)
    if selection == 'order':
        my_table = reduce(lambda x, y: pd.merge(x, y, on='Order'), frames)
    my_table = my_table.round(2)
    ax.table(cellText=my_table.values,
             colLabels=my_table.columns,
             loc='center', rowLabels=my_table.index)
    plt.title('Mean TE Length (bp)')
    fig.tight_layout()
    plt.savefig(output_dir + '/' + Camarosa.genome_id +
                '_' + selection + '_Mean_Subgenome_Table.png')
    plt.close()


def graph_table_median_TE_lengths(Camarosa, list_of_genome_data,
                                  selection,
                                  output_dir):
    # TODO Edit the arg list and try to reduce the number of ifs.
    """
    Args:
    """

    fig, ax = plt.subplots()
    ax.axis('off')
    ax.axis('tight')

    frames = []
    for GenomeData in list_of_genome_data:
        if selection == 'superfam':
            median_length_dict = GenomeData.median_superfam_length.to_dict()
        if selection == 'order':
            median_length_dict = GenomeData.median_order_length.to_dict()
        subset_df = \
        pd.DataFrame.from_dict(median_length_dict,
                               orient='index',
                               columns=[GenomeData.subgenome_identity])
        if selection == 'superfam':
            subset_df.index.name = 'SuperFamily'
        if selection == 'order':
            subset_df.index.name = 'Order'
        frames.append(subset_df)

    if selection == 'superfam':
        median_length_dict = Camarosa.median_superfam_length.to_dict()
    if selection == 'order':
        median_length_dict = Camarosa.median_order_length.to_dict()
    Cam_data = pd.DataFrame.from_dict(median_length_dict,
                                      orient='index',
                                      columns=[Camarosa.genome_id])
    if selection == 'superfam':
        Cam_data.index.name = 'SuperFamily'
    if selection == 'order':
        Cam_data.index.name = 'Order'
    frames.append(Cam_data)

    if selection == 'superfam':
        my_table = reduce(lambda x, y: pd.merge(x, y, on='SuperFamily'), frames)
    if selection == 'order':
        my_table = reduce(lambda x, y: pd.merge(x, y, on='Order'), frames)
    my_table = my_table.round(2)
    ax.table(cellText=my_table.values,
             colLabels=my_table.columns,
             loc='center', rowLabels=my_table.index)
    plt.title('Median TE Length (bp)')
    fig.tight_layout()
    plt.savefig(output_dir + '/' + Camarosa.genome_id +
                '_' + selection + '_Median_Subgenome_Table.png')
    plt.close()


if __name__ == '__main__':
    """Command line interface to generate pie graphs."""
    parser = argparse.ArgumentParser()
    path_main = os.path.abspath(__file__)
    parser.add_argument('--bar_output_dir', '-b', type=str,
                        default=os.path.join(path_main, '../../results/bar'),
                        help='Parent directory to output bar results')
    parser.add_argument('--pie_output_dir', '-p', type=str,
                        default=os.path.join(path_main, '../../results/pie'),
                        help='Parent directory to output pie results')
    parser.add_argument('--table_output_dir', '-ta', type=str,
                        default=os.path.join(path_main, '../../results/table'),
                        help='Parent directory to output table results')

    args = parser.parse_args()
    args.bar_output_dir = os.path.abspath(args.bar_output_dir)
    args.pie_output_dir = os.path.abspath(args.pie_output_dir)
    args.table_output_dir = os.path.abspath(args.table_output_dir)
    # Load Data

    data_place = '/home/scott/Documents/Uni/Research/Projects/TE_Density/filtered_input_data/'

    # NOTE Genome Sizes in MB:
    cam_size = 0.805489
    h4_size = 0.250
    fda_size = 0.29079
    fii_size = 0.26556
    fma_size = 0.26575
    fng_size = 0.30590
    fpe_size = 0.28233
    fvi_size = 0.22961

    # NOTE Camarosa Information, Import Camarosa Data:
    Gene_Data_Cam = pd.read_csv(os.path.join(data_place,
                                'Cleaned_Cam_Genes.tsv'),
                                header='infer', sep='\t',
                                index_col='Gene_Name')
    TE_Data_Cam = pd.read_csv(os.path.join(data_place,
                              'Cleaned_Cam_TEs.tsv'),
                              header='infer', sep='\t')
    Camarosa = GenomeData('Camarosa', Gene_Data_Cam, TE_Data_Cam, cam_size)
    # NOTE Camarosa subgenome chromosomes
    F_vesca = ['Fvb1-4', 'Fvb2-2', 'Fvb3-4', 'Fvb4-3', 'Fvb5-1', 'Fvb6-1', 'Fvb7-2']
    F_nipponica = ['Fvb1-3', 'Fvb2-1', 'Fvb3-3', 'Fvb4-2', 'Fvb5-4', 'Fvb6-2', 'Fvb7-1']
    F_iinumae = ['Fvb1-2', 'Fvb2-4', 'Fvb3-2', 'Fvb4-4', 'Fvb5-4', 'Fvb6-3', 'Fvb7-3']
    F_viridis = ['Fvb1-1', 'Fvb2-3', 'Fvb3-1', 'Fvb4-1', 'Fvb5-2', 'Fvb6-4', 'Fvb7-4']

    Vesca = SubgenomeData('Camarosa', Gene_Data_Cam,
                          TE_Data_Cam, F_vesca, 'F_vesca')
    Nipponica = SubgenomeData('Camarosa', Gene_Data_Cam,
                              TE_Data_Cam, F_nipponica, 'F_nipponica')
    Iinumae = SubgenomeData('Camarosa', Gene_Data_Cam,
                            TE_Data_Cam, F_iinumae, 'F_iinumae')
    Viridis = SubgenomeData('Camarosa', Gene_Data_Cam,
                            TE_Data_Cam, F_viridis, 'F_viridis')
    list_of_subgenomes = [Vesca, Nipponica, Iinumae, Viridis]

    # NOTE H4 Information, Import H4 Data:
    Gene_Data_H4 = pd.read_csv(os.path.join(data_place,
                               'Cleaned_H4_Genes.tsv'),
                               header='infer', sep='\t',
                               index_col='Gene_Name')
    TE_Data_H4 = pd.read_csv(os.path.join(data_place,
                             'Cleaned_H4_TEs.tsv'),
                             header='infer', sep='\t')
    H4 = GenomeData('H4', Gene_Data_H4, TE_Data_H4, h4_size)

    graph_general_genome_stats_pie(Camarosa, args.pie_output_dir, explode_bool=False)
    graph_order_sequences_as_percents(Camarosa, args.pie_output_dir)
    graph_superfam_sequences_as_percents(Camarosa, args.pie_output_dir)
    graph_order_numbers_as_percents(Camarosa, args.pie_output_dir)
    graph_superfam_numbers_as_percents(Camarosa, args.pie_output_dir)
    graph_cam_superfam_chromosomes(Camarosa, args.pie_output_dir,
                                   F_vesca, F_iinumae,
                                   F_viridis, F_nipponica)

    graph_general_genome_stats_pie(H4, args.pie_output_dir, explode_bool=False)
    graph_order_numbers_as_percents(H4, args.pie_output_dir)
    graph_order_sequences_as_percents(H4, args.pie_output_dir)
    graph_superfam_numbers_as_percents(H4, args.pie_output_dir)
    graph_superfam_sequences_as_percents(H4, args.pie_output_dir)

    graph_bar_superfam_lengths_genome_wide(list_of_subgenomes,
                                           'Camarosa_Subgenome_Content_Bar',
                                           args.bar_output_dir)

    graph_average_LTR_length_H4_Cam(H4, Camarosa, args.bar_output_dir)

    graph_table_superfam_lengths_genome_wide(Camarosa, list_of_subgenomes,
                                             args.table_output_dir)

    graph_table_mean_TE_lengths(Camarosa, list_of_subgenomes,
                                'order',
                                args.table_output_dir)

    graph_table_mean_TE_lengths(Camarosa, list_of_subgenomes,
                                'superfam',
                                args.table_output_dir)

    graph_table_median_TE_lengths(Camarosa, list_of_subgenomes,
                                  'order',
                                  args.table_output_dir)

    graph_table_median_TE_lengths(Camarosa, list_of_subgenomes,
                                  'superfam',
                                  args.table_output_dir)
