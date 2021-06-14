def te_annot_renamer(TE_Data):
    U = "Unknown_Order"
    master_order = {
        "Unknown": U,
        "MITE": "TIR",
        "pararetrovirus": "pararetrovirus",
        "DNA": "TIR",
    }

    U = "Unknown_Superfam"
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
        "DTT": "Tc1-Mariner",
        "DTA": "hAT",
        "DTM": "Mutator",
        "DTE": "Merlin",
        "DTR": "Transib",
        "DTP": "P",
        "DTB": "PiggyBac",
        "DTH": "PIF-Harbinger",
        "DTC": "CACTA",
        "DYC": "Crypton",
        "DHH": "Helitron",
        "DMM": "Maverick",
        # Custom changes
        "unknown": U,
        "Unknown": U,
        "None": U,
        "EnSpm_CACTA": "CACTA",
        "MuDR_Mutator": "Mutator",
        "PIF_Harbinger": "PIF-Harbinger",
    }

    TE_Data.SuperFamily.fillna(
        value="Unknown_Superfam", inplace=True
    )  # replace None w U

    # Invoke dictionary to fix names
    TE_Data.Order.replace(master_order, inplace=True)
    TE_Data.SuperFamily.replace(master_superfamily, inplace=True)

    # Rename the superfamily value for pararetros as pararetrovirus
    TE_Data.loc[TE_Data.Order == "pararetrovirus", "SuperFamily"] = "pararetrovirus"

    # Rename unknown LINE element superfamilies to Unknown_LINE_Superfam to
    # distinguish between other unknowns
    TE_Data.loc[
        (TE_Data.Order == "LINE") & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        "SuperFamily",
    ] = "Unknown_LINE_Superfam"

    # Rename unknown LTR element superfamilies to Unknown_LTR_Superfam to
    # distinguish between other unknowns
    TE_Data.loc[
        (TE_Data["Order"] == "LTR") & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        "SuperFamily",
    ] = "Unknown_LTR_Superfam"

    # Rename unknown TIR element superfamilies to Unknown_TIR_Superfam to
    # distinguish between other unknowns
    TE_Data.loc[
        (TE_Data.Order == "TIR") & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        "SuperFamily",
    ] = "Unknown_TIR_Superfam"

    # Rename both values for Helitron elements, so that 'Helitron' is
    # both the Order and SuperFamily value
    # Some Helitron elements were labeled 'DNA' in the Order location, this is
    # technically correct but I prefer to differentiate the TIR DNA elements
    # from DNA elements as a whole
    TE_Data.loc[
        (TE_Data["Order"] == "TIR") & (TE_Data["SuperFamily"] == "Helitron"),
        ["Order", "SuperFamily"],
    ] = "Helitron"
    # If the Order is Helitron and the SuperFamily is unknown make the
    # superfamily 'Helitron'
    TE_Data.loc[
        (TE_Data["Order"] == "Helitron")
        & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        "SuperFamily",
    ] = "Helitron"

    # For TEs that are unknown for both Order AND SuperFamily we will call
    # those 'Completely_Unknown'
    TE_Data.loc[
        (TE_Data["Order"] == "Unknown_Order")
        & (TE_Data["SuperFamily"] == "Unknown_Superfam"),
        ["Order", "SuperFamily"],
    ] = "Completely_Unknown"

    return TE_Data


def diagnostic_cleaner_helper(TE_Data):
    print()
    print(TE_Data.Order.unique())
    print(TE_Data.SuperFamily.unique())
    print()

    # To see unique for a given type:
    # print(TE_Data.loc[TE_Data['Order'] == 'LINE'].SuperFamily.unique())
    return None
