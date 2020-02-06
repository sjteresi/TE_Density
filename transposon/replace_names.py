def TE_Renamer(TE_Dataframe):
    U = 'Unknown'
    master_order = {
        'RC?':'DNA',
        'RC':'DNA',
        'SINE?':U,
        'tandem':'Tandem',
        'No_hits':U
    }

    U = 'Unknown_SuperFam'
    master_superfamily = {
        'Uknown':U,
        'MuDr':'MULE',
        'MULE-MuDR':'MULE',
        'Mutator|cleanup':'MULE',
        'TcMar':U,
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
        'Micro_like':U,
        'Micro-like-sequence':U,
        'Micro-like-sequence|cleanup':U,
        'Unclassified':U,
        'L1-Tx1':'LINE',
        'CRE':'LINE',
        'CACTA':'CMC-EnSpm',
        'Tad1':U,
        'hAT|cleanup':'hAT',
        '':U,
        'Line':'LINE'
    }

    TE_Dataframe.Order.replace(master_order, inplace=True)
    TE_Dataframe.SuperFamily.replace(master_superfamily, inplace=True)
    return TE_Dataframe
