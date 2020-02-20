def TE_Renamer(TE_Data):
    U = 'Unknown_Order'
    master_order = {
        # Custom changes
        'unknown': 'Unknown',
        'MITE': 'DNA_(TIR)',
        'RC?': 'DNA',
        'RC': 'DNA',
        'SINE?': U,
        'tandem': 'Tandem',
        'No_hits': U,

        # EDTA/Wicker et al 2007 renames to common name:
        'LTR': 'LTR',
        'DIRS': 'DIRS',
        'PLE': 'PLE',
        'SINE': 'SINE',
        'LINE': 'LINE',
        'TIR': 'DNA_(TIR)',
        'Crypton': 'Crypton',
        'Helitron': 'Helitron',
        'Maverick': 'Maverick'

    }

    U = 'Unknown_SuperFam'
    master_superfamily = {
        # EDTA/Wicker et al 2007 renames to common name:
        'RLC': 'Copia',
        'RLG': 'Gypsy',
        'RLB': 'Bel_Pao',
        'RLR': 'Retrovirus',
        'RLE': 'ERV',
        'RYD': 'DIRS',
        'RYN': 'Ngaro',
        'RYV': 'VIPER',
        'RPP': 'Penelope',
        'RIR': 'R2',
        'RIT': 'RTE',
        'RIJ': 'Jockey',
        'RIL': 'L1',
        'RII': 'I',
        'RST': 'tRNA',
        'RSL': '7SL',
        'RSS': '5S',

        'DTT': 'Tc1_Mariner',
        'DTA': 'hAT',
        'DTM': 'Mutator',
        'DTE': 'Merlin',
        'DTR': 'Transib',
        'DTP': 'P',
        'DTB': 'PiggyBac',
        'DTH': 'PIF_Harbinger',
        'DTC': 'CACTA',
        'DYC': 'Crypton',
        'DHH': 'Helitron',
        'DMM': 'Maverick',


        # Custom changes
        'Uknown': U,
        'unknown': U,
        'Unknown': U,
        'MuDr': 'MULE',
        'MULE-MuDR': 'MULE',
        'Mutator|cleanup': 'MULE',
        'TcMar': U,
        'Pao': U,
        'Caulimovirus': U,
        'hAT-Tag1': 'hAT',
        'hAT-Tip100': 'hAT',
        'hAT-Charlie': 'hAT',
        'Helitron': U,
        'Maverick': U,
        'Harbinger': 'PIF_Harbinger',
        'TcMar-Pogo': U,
        'CR1': 'LINE',
        'hAT-Ac': 'hAT',
        'L2': 'LINE',
        'L1': 'LINE',
        'Jockey': 'LINE',
        'MuLE-MuDR': 'MULE',
        'MuDR': 'MULE',
        'Mutator': 'MULE',
        'Micro_like': U,
        'Micro-like-sequence': U,
        'Micro-like-sequence|cleanup': U,
        'Unclassified': U,
        'L1-Tx1': 'LINE',
        'CRE': 'LINE',
        'CACTA': 'CMC-EnSpm',
        'Tad1': U,
        'hAT|cleanup': 'hAT',
        '': U,
        'Line': 'LINE'
    }
    TE_Data.SuperFamily.fillna(value='Unknown_SuperFam', inplace=True)  # replace None w U
    # step to fix TE names
    TE_Data.Order.replace(master_order, inplace=True)
    TE_Data.SuperFamily.replace(master_superfamily, inplace=True)
    TE_Data.loc[TE_Data.Order == 'Tandem', 'SuperFamily'] = 'Tandem'

    to_drop = TE_Data.Chromosome.str.contains('##sequence-region')
    TE_Data = TE_Data[~to_drop]
    to_drop = TE_Data.Chromosome.str.contains('contig*')
    TE_Data = TE_Data[~to_drop]

    TE_Data = TE_Data[TE_Data.Order != 'Simple_repeat']  # drop s repeat
    return TE_Data
