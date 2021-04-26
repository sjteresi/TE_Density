def te_annot_renamer(TE_Data):
    U = "Unknown_Order"
    master_order = {
        # Custom changes
        ## RepeatMasker-based Changes
        "unknown": U,
        "Unknown": U,
        "MITE": "DNA",
        "RC?": "DNA",
        "RC": "DNA",
        "SINE?": U,
        "tandem": "Tandem",
        "No_hits": U,
        ## EDTA-based Changes
        "pararetrovirus": "LTR",
        "mixture": "Mixture",
        "DNA": "TIR",
    }

    U = "Unknown_SuperFam"
    master_superfamily = {
        # EDTA/Wicker et al 2007 renames to common name:
        "RLC": "Copia",
        "RLG": "Gypsy",
        "RLB": "Bel_Pao",
        "RLR": "Retrovirus",
        "RLE": "ERV",
        "RYD": "DIRS",
        "RYN": "Ngaro",
        "RYV": "VIPER",
        "RPP": "Penelope",
        "RIR": "R2",
        "RIT": "RTE",
        "RIJ": "Jockey",
        "RIL": "L1",
        "RII": "I",
        "RST": "tRNA",
        "RSL": "7SL",
        "RSS": "5S",
        "DTT": "Tc1_Mariner",
        "DTA": "hAT",
        "DTM": "Mutator",
        "DTE": "Merlin",
        "DTR": "Transib",
        "DTP": "P",
        "DTB": "PiggyBac",
        "DTH": "PIF_Harbinger",
        "DTC": "CACTA",
        "DYC": "Crypton",
        "DHH": "Helitron",
        "DMM": "Maverick",
        # Custom changes
        "Uknown": U,
        "unknown": U,
        "Unknown": U,
        "EnSpm_CACTA": "CACTA"
        #'MuDr': 'MULE',
        #'MULE-MuDR': 'MULE',
        #'Mutator|cleanup': 'MULE',
        #'TcMar': U,
        #'Pao': U,
        #'Caulimovirus': U,
        #'hAT-Tag1': 'hAT',
        #'hAT-Tip100': 'hAT',
        #'hAT-Charlie': 'hAT',
        #'Helitron': U,
        #'Maverick': U,
        #'Harbinger': 'PIF_Harbinger',
        #'PIF-Harbinger': 'PIF_Harbinger',
        #'TcMar-Pogo': U,
        #'CR1': 'LINE',
        #'hAT-Ac': 'hAT',
        #'L2': 'LINE',
        #'L1': 'LINE',
        #'Jockey': 'LINE',
        #'MuLE-MuDR': 'MULE',
        #'MuDR': 'MULE',
        #'Mutator': 'MULE',
        #'Micro_like': U,
        #'Micro-like-sequence': U,
        #'Micro-like-sequence|cleanup': U,
        #'Unclassified': U,
        #'L1-Tx1': 'LINE',
        #'CRE': 'LINE',
        #'CACTA': 'CMC-EnSpm',
        #'Tad1': U,
        #'hAT|cleanup': 'hAT',
        #'': U,
        #'Line': 'LINE'
    }
    TE_Data.SuperFamily.fillna(
        value="Unknown_SuperFam", inplace=True
    )  # replace None w U
    # step to fix TE names
    TE_Data.Order.replace(master_order, inplace=True)
    TE_Data.SuperFamily.replace(master_superfamily, inplace=True)
    TE_Data.loc[TE_Data.Order == "Tandem", "SuperFamily"] = "Tandem"

    # to_drop = TE_Data.Chromosome.str.contains('##sequence-region')
    # TE_Data = TE_Data[~to_drop]
    # to_drop = TE_Data.Chromosome.str.contains('contig*')
    # TE_Data = TE_Data[~to_drop]

    TE_Data.loc[
        (TE_Data["Order"] == "Unknown_Order")
        & (TE_Data["SuperFamily"] == "Unknown_SuperFam"),
        ["Order", "SuperFamily"],
    ] = "Completely_Unknown"

    TE_Data.loc[
        (TE_Data["Order"] == "Helitron")
        & (TE_Data["SuperFamily"] == "Unknown_SuperFam"),
        ["SuperFamily"],
    ] = "Helitron"
    TE_Data.loc[(TE_Data["Order"] == "Helitron"), ["Order"]] = "Helitron"
    TE_Data.loc[(TE_Data["SuperFamily"] == "Helitron"), ["Order"]] = "Helitron"
    TE_Data.loc[(TE_Data["Order"] == "Mixture"), ["SuperFamily"]] = "Mixture"

    ltr_elements = ["Copia", "Gypsy"]
    TE_Data.loc[
        (TE_Data["Order"] == "LTR") & (~TE_Data["SuperFamily"].isin(ltr_elements)),
        ["SuperFamily"],
    ] = "Unknown_LTR_Superfam"

    TE_Data = TE_Data[TE_Data.Order != "Simple_repeat"]  # drop s repeat
    TE_Data = TE_Data[TE_Data.Order != "long_terminal_repeat"]  # drop
    TE_Data = TE_Data[TE_Data.Order != "Maverick"]  # drop if in Order category
    TE_Data = TE_Data[
        TE_Data.Order != "target_site_duplication"
    ]  # drop if in Order category
    return TE_Data
