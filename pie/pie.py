#!/usr/bin/env python3

"""
Generate some pie graphs for showing TE distributions
"""

__author__ = "Scott Teresi"

import os
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import coloredlogs
import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)


from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.import_genes import import_genes
from transposon.import_transposons import import_transposons
from transposon.replace_names import te_annot_renamer
from transposon.density import verify_TE_cache, verify_gene_cache


def get_nulls(my_df):
    null_columns = my_df.columns[my_df.isnull().any()]
    count_of_null = my_df[null_columns].isnull().sum()
    print('Counts of null values per column: ' '\n', count_of_null, '\n')
    rows_where_null = my_df[my_df.isnull().any(axis=1)][null_columns].head()
    print('Rows where null exist: ', '\n', rows_where_null, '\n')


def drop_nulls(my_df, status=False):
    if status:
        print('DROPPING ROWS WITH AT LEAST ONE NULL VALUE!!!')
    my_df = my_df.dropna(axis=0, how='any')
    return my_df


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


def sum_genes_and_TEs(Gene_Data, TE_Data, selection):
    if selection == 'Camarosa':
        genome_size = 0.805489
    elif selection == 'H4':
        genome_size = 0.250
    elif selection == 'FDA':
        genome_size = 0.29079
    elif selection == 'FII':
        genome_size = 0.26556
    elif selection == 'FMA':
        genome_size = 0.26575
    elif selection == 'FNG':
        genome_size = 0.30590
    elif selection == 'FPE':
        genome_size = 0.28233
    elif selection == 'FVI':
        genome_size = 0.22961
    else:
        raise ValueError('What genome are you using?')
    gene_lengths = Gene_Data.Length.sum()
    TE_lengths = TE_Data.Length.sum()

    # Convert to Gbp
    gene_lengths = gene_lengths / 1000000000 / genome_size  # normalize
    TE_lengths = TE_lengths / 1000000000 / genome_size  # normalize
    other = 1 - gene_lengths - TE_lengths
    my_lengths = []
    my_lengths.append(other)
    my_lengths.append(gene_lengths)
    my_lengths.append(TE_lengths)
    return my_lengths


def general_genome_stats(Gene_Data, TE_Data, output_dir, selection):
    # order_labels = get_unique(TE_Data.Order)  # list
    # superfamily_labels = get_unique(TE_Data.SuperFamily)  # list
    sizes = sum_genes_and_TEs(Gene_Data, TE_Data, selection)
    # total_num_TEs = len(TE_Data.index)

    # Plotting
    # Genome Content Plot
    labels = ['Other', 'Genes', 'TEs']
    new_sizes = []
    for size in sizes:
        new_sizes.append(size * 100)
    colors = ['forestgreen', 'cornflowerblue', 'crimson']
    explode = []
    for i in new_sizes:
        explode.append(0)
    # explode[1] = 0.1
    explode = tuple(explode)
    patches, texts, autotext = plt.pie(new_sizes, colors=colors, shadow=False,
                                       startangle=120, explode=explode,
                                       autopct='%1.1f%%', labels=labels)
    plt.legend(patches, labels, loc='best')
    plt.title('Contribution to Genome Size')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output_dir + '/' + selection + '_General_Genome_Content.png')
    plt.close()


def te_order_content(TE_Data, output_dir, selection):
    # TODO get rid of the family notation
    my_family_counts = TE_Data.Order.value_counts()
    my_family_counts = my_family_counts.to_dict()
    family_sum = 0
    for key, val in my_family_counts.items():
        family_sum += val
    new_sizes = []
    for key, val in my_family_counts.items():
        my_family_counts[key] = val/family_sum*100  # normalize
        new_sizes.append(my_family_counts[key])

    labels = my_family_counts.keys()
    colors = ['seagreen', 'deepskyblue', 'lightcoral', 'darkviolet', 'crimson']
    explode = []
    for i in new_sizes:
        explode.append(0)
    # explode[0] = 0.1
    explode = tuple(explode)
    patches, texts, autotext = plt.pie(new_sizes, colors=colors, shadow=False,
                                       startangle=150, explode=explode,
                                       autopct='%1.1f%%', labels=labels)
    plt.legend(patches, labels, loc='best')
    plt.title('Relative %s of TEs Compared to Total Number of TEs (Order)')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output_dir + '/' + selection + '_TE_Average_Order_Content.png')
    plt.close()


def te_subfam_content(TE_Data, output_dir, selection):
    # TE Distro Plot
    # subfamily_labels = get_unique(TE_Data.SuperFamily)  # list
    my_subfamily_counts = TE_Data.SuperFamily.value_counts()
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
    plt.savefig(output_dir + '/' + selection + '_TE_Average_SuperFam_Content.png')
    plt.close()


def graphs_cam_chromosomes(TE_Data, output_dir):
    selection = 'Camarosa'
    # Loop
    # n = 2.9
    # fig.subplots_adjust(hspace = n, wspace = n)
    # n = 4
    colors = ['sienna', 'goldenrod', 'lightgreen', 'forestgreen',
              'steelblue', 'darkslateblue', 'orchid', 'crimson', 'mistyrose',
              'darkmagenta']

    # explode = (0, 0.1, 0.1, 0, 0, 0, 0, 0, 0)  # explode

    all_c = get_unique(TE_Data.Chromosome)  # list
    grouped_TEs = split(TE_Data, 'Chromosome')

    for xi in range(0, len(all_c), 4):
        fig, axes = plt.subplots(nrows=2, ncols=2)
        # print(f"XI is: {xi}")
        # raise ValueError

        for i, ax in enumerate(axes.ravel()):
            # print(f"i is: {i}")
            # list_subset = all_c[xi: xi + i + 1]
            titles = get_unique(grouped_TEs[xi + i].Chromosome)  # list
            # print(f"Title is: {titles} \n")
            size_vals = []
            my_subfamily_counts = grouped_TEs[xi+i].SuperFamily.value_counts()
            my_subfamily_counts = my_subfamily_counts.to_dict()
            subfamily_sum = 0

            for key, val in my_subfamily_counts.items():
                subfamily_sum += val
            for key, val in my_subfamily_counts.items():
                my_subfamily_counts[key] = val/subfamily_sum*100  # normalize
                size_vals.append(my_subfamily_counts[key])

            labels = my_subfamily_counts.keys()
            ax.pie(size_vals,
                   labels=None,
                   autopct='%1.1f%%',
                   shadow=False,
                   startangle=90,
                   radius=1.7,
                   pctdistance=0.84,
                   colors=colors,
                   textprops={'fontsize': 6.8})
            ax.set_title(titles[0], loc='left', pad=19)
            plt.tight_layout()
        # n = 5
        # plt.tight_layout(pad=n, w_pad = n, h_pad = n)
        fig.legend(labels=labels,
                   loc='center right',
                   title='TEs',
                   borderaxespad=-0.2)
        n = 0.9
        plt.subplots_adjust(right=0.75, hspace=n)
        plt.savefig(output_dir + '/' + selection +
                    '_' + all_c[xi] + '_All_Chromosomes.png')
        # plt.savefig(all_c[xi]+'_Cam_All.png')
        plt.close()


def graphs_diploid_chromosomes(TE_Data, output_dir, selection):
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


def validate_args(args, logger):
    """Raise if an input argument is invalid."""

    if not os.path.isfile(args.genes_input_file):
        logger.critical("argument 'genes_input_dir' is not a file")
        raise ValueError("%s is not a directory" % (args.genes_input_file))
    if not os.path.isfile(args.tes_input_file):
        logger.critical("argument 'tes_input_dir' is not a file")
        raise ValueError("%s is not a directory" % (args.tes_input_file))
    if not os.path.isdir(args.output_dir):
        logger.critical("argument 'output_dir' is not a directory")
        raise ValueError("%s is not a directory" % (args.output_dir))


if __name__ == '__main__':
    """Command line interface to generate pie graphs."""

    parser = argparse.ArgumentParser(description="""
                                     Generate pie graphs of TE distribtutions
                                     """)
    path_main = os.path.abspath(__file__)
    parser.add_argument('genes_input_file', type=str,
                        help='parent path of gene file')
    parser.add_argument('tes_input_file', type=str,
                        help='parent path of transposon file')
    parser.add_argument('--contig_del', default=True)
    parser.add_argument('--filtered_input_data', '-f', type=str,
                        default=os.path.join(path_main, '../..',
                                             'filtered_input_data'),
                        help='parent directory for cached input data')
    parser.add_argument('selection', type=str,
                        help='Genome')
    parser.add_argument('--output_dir', '-o', type=str,
                        default=os.path.join(path_main, '../../results/pie'),
                        help='parent directory to output results')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='set debugging level to DEBUG')

    args = parser.parse_args()
    args.genes_input_file = os.path.abspath(args.genes_input_file)
    args.tes_input_file = os.path.abspath(args.tes_input_file)
    args.selection = str(args.selection)
    args.filtered_input_data = os.path.abspath(args.filtered_input_data)
    args.output_dir = os.path.abspath(args.output_dir)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)

    for argname, argval in vars(args).items():
        logger.debug("%-12s: %s" % (argname, argval))
    validate_args(args, logger)

    # Load Data
    logger.info("Checking disk for previously filtered data...")
    g_fname = os.path.basename(os.path.splitext(args.genes_input_file)[0])
    t_fname = os.path.basename(os.path.splitext(args.tes_input_file)[0])
    cleaned_genes = os.path.join(args.filtered_input_data, str('Cleaned_' +
                                                               g_fname +
                                                               '.tsv'))
    cleaned_transposons = os.path.join(args.filtered_input_data, str('Cleaned_' +
                                                                     t_fname +
                                                                     '.tsv'))
    Gene_Data = verify_gene_cache(args.genes_input_file, cleaned_genes,
                                  args.contig_del, logger)
    TE_Data = verify_TE_cache(args.tes_input_file, cleaned_transposons,
                              te_annot_renamer, args.contig_del, logger)

    # Generate Graphs
    general_genome_stats(Gene_Data, TE_Data, args.output_dir, args.selection)
    te_order_content(TE_Data, args.output_dir, args.selection)
    te_subfam_content(TE_Data, args.output_dir, args.selection)
    if args.selection == 'Camarosa':
        graphs_cam_chromosomes(TE_Data, args.output_dir)
    else:
        graphs_diploid_chromosomes(TE_Data, args.output_dir, args.selection)
