import pandas as pd
import os
from transposon.replace_names import TE_Renamer

def check_nulls(my_df):
    # TODO add logging event and possibly crash
    Bool = my_df.isnull().values.any()
    if Bool:
        print('You have Null values in your dataframe that were not caught!')

def import_transposons(input_dir):
    """Import TE File
        Args: input_dir (command line argument) Specify the input directory of
        the TE annotation data, this is the same as the Gene annotation
        directory
    """

    # TODO remove MAGIC NUMBER (perhaps just search by extension (gtf)?)
    gff_filename = 'camarosa_gff_data.gff' # DECLARE YOUR DATA NAME

    col_names = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
        'Score', 'Strand', 'Frame', 'Attribute']

    col_to_use = ['Chromosome', 'Software', 'Feature', 'Start', 'Stop', \
                 'Strand']

    TE_Data = pd.read_csv(
            os.path.join(input_dir, gff_filename),
            sep='\t+',
            header=None,
            engine='python',
            names = col_names,
            usecols = col_to_use)

    TE_Data[['Order', 'SuperFamily']] = TE_Data.Feature.str.split('/', expand=True)

    to_drop = TE_Data[TE_Data.Chromosome.str.contains('##sequence-region', case=False)]
    TE_Data.drop(to_drop.index, inplace = True)
    TE_Data.SuperFamily.fillna(value='Unknown_SuperFam', inplace=True) # replace None w U
        # step to fix TE names

    TE_Data = TE_Data.drop(['Feature', 'Software'], axis=1)


    TE_Data.Strand = TE_Data.Strand.astype(str)
    TE_Data.Start = TE_Data.Start.astype('uint32')
    TE_Data.Stop = TE_Data.Stop.astype('uint32')
    TE_Data['Length'] = TE_Data.Stop - TE_Data.Start + 1
    TE_Data = TE_Data[TE_Data.Order != 'Simple_repeat'] # drop s repeat
    TE_Data = TE_Renamer(TE_Data)
    check_nulls(TE_Data)
    return TE_Data
