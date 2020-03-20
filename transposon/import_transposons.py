import pandas as pd


def check_nulls(my_df):
    # TODO add logging event and possibly crash
    Bool = my_df.isnull().values.any()
    if Bool:
        print('You have Null values in your dataframe that were not caught!')


# TODO SCOTT add an input argument for filtering the data, so a client can
# provide a filtering function (TE_Renamer) instead of the import hard coded
# that way the client doesn't have to modify your code to get it to work
# also use lower case for function names and start only classes with capital
# letters
def import_transposons(tes_input_path, te_annot_renamer, contig_del):
    """Import TE File.
        Args: input_dir (command line argument) Specify the input directory of
        the TE annotation data, this is the same as the Gene annotation
        directory
        
        contig_drop (bool): logical whether to drop rows with a contig as the
        chromosome id
    """
    col_names = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop',
                 'Score', 'Strand', 'Frame', 'Attribute']

    col_to_use = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop',
                  'Strand']

    TE_Data = pd.read_csv(
        tes_input_path,
        sep='\t+',
        header=None,
        engine='python',
        names=col_names,
        usecols=col_to_use)

    TE_Data = TE_Data[~TE_Data.Chromosome.str.contains('#')]  # remove comment
    # rows in annotation
    TE_Data[['Order', 'SuperFamily']] = TE_Data.Feature.str.split('/', expand=True)

    TE_Data = TE_Data.drop(['Feature', 'Software'], axis=1)
    TE_Data = te_annot_renamer(TE_Data)  # NOTE call to the cleaner
    TE_Data.Strand = TE_Data.Strand.astype(str)
    TE_Data.Start = TE_Data.Start.astype('uint32')
    TE_Data.Stop = TE_Data.Stop.astype('uint32')
    TE_Data['Length'] = TE_Data.Stop - TE_Data.Start + 1
    if contig_del:
        TE_Data = TE_Data[~TE_Data.Chromosome.str.contains('contig')]
    check_nulls(TE_Data)
    return TE_Data
