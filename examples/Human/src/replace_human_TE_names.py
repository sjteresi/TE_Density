def te_annot_renamer(transposon_data):
    # TODO clean up all comments

    print()
    print(transposon_data.shape)
    print(transposon_data["Order"].unique())
    print(transposon_data["SuperFamily"].unique())
    print()
    # transposon_data.SuperFamily.fillna(
    # value="Unknown_SuperFam", inplace=True
    # )  # replace None w U
    # step to fix TE names
    # transposon_data.Order.replace(master_order, inplace=True)
    # transposon_data.SuperFamily.replace(master_superfamily, inplace=True)
    # transposon_data.loc[transposon_data.Order == "Tandem", "SuperFamily"] = "Tandem"

    #################################################
    # Make Penelope its own Order called PLE, make the SuperFamily Penelope
    # Corresponds to Wicker's grouping of Penelope elements
    transposon_data.loc[
        (transposon_data["Order"] == "LINE")
        & (transposon_data["SuperFamily"] == "Penelope"),
        ["Order"],
    ] = "PLE"

    # Drop CR1 elements from dataset, not in Wicker
    # transposon_data = transposon_data[transposon_data.SuperFamily != "CR1"]

    # Drop Dong-R4 elements from dataset, not in Wicker
    # transposon_data = transposon_data[transposon_data.SuperFamily != "Dong-R4"]

    # Drop all RNA related TE Orders
    transposon_data = transposon_data[transposon_data.Order != "snRNA"]
    transposon_data = transposon_data[transposon_data.Order != "tRNA"]
    transposon_data = transposon_data[transposon_data.Order != "rRNA"]
    transposon_data = transposon_data[transposon_data.Order != "rRNA"]
    transposon_data = transposon_data[transposon_data.Order != "srpRNA"]
    transposon_data = transposon_data[transposon_data.Order != "scRNA"]

    # Remove ERV
    transposon_data = transposon_data[transposon_data.SuperFamily != "ERVL-maLR"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "ERV1"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "ERVL"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "ERVK"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "ERV1?"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "ERVL?"]
    # Remove low confidence and other minor/misc things
    transposon_data = transposon_data[transposon_data.Order != "LTR?"]
    transposon_data = transposon_data[transposon_data.Order != "SINE?"]
    transposon_data = transposon_data[transposon_data.Order != "RC?"]
    transposon_data = transposon_data[transposon_data.Order != "DNA?"]
    transposon_data = transposon_data[transposon_data.Order != "RNA"]
    transposon_data = transposon_data[transposon_data.Order != "Low_complexity"]
    transposon_data = transposon_data[transposon_data.Order != "RC"]
    transposon_data = transposon_data[transposon_data.Order != "Satellite"]
    transposon_data = transposon_data[transposon_data.Order != "Simple_repeat"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "PiggyBac?"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "Gypsy?"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "TcMar?"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "hAT?"]
    transposon_data = transposon_data[transposon_data.SuperFamily != "hAT-Tip100?"]

    # Remove elements that are DNA as superfamily
    transposon_data = transposon_data[transposon_data.SuperFamily != "DNA"]

    # Rename DNA order to TIR
    transposon_data.loc[(transposon_data["Order"] == "DNA"), ["Order"]] = "TIR"

    ###################################################3
    # to_drop = transposon_data.Chromosome.str.contains('##sequence-region')
    # transposon_data = transposon_data[~to_drop]
    # to_drop = transposon_data.Chromosome.str.contains('contig*')
    # transposon_data = transposon_data[~to_drop]

    # transposon_data.loc[
    # (transposon_data["Order"] == "Unknown_Order")
    # & (transposon_data["SuperFamily"] == "Unknown_SuperFam"),
    # ["Order", "SuperFamily"],
    # ] = "Completely_Unknown"

    # transposon_data.loc[
    # (transposon_data["Order"] == "Helitron")
    # & (transposon_data["SuperFamily"] == "Unknown_SuperFam"),
    # ["SuperFamily"],
    # ] = "Helitron"
    # transposon_data.loc[(transposon_data["SuperFamily"] == "Helitron"), ["Order"]] = "Helitron"
    # transposon_data.loc[(transposon_data["Order"] == "Mixture"), ["SuperFamily"]] = "Mixture"

    # ltr_elements = ["Copia", "Gypsy"]
    # transposon_data.loc[
    # (transposon_data["Order"] == "LTR") & (~transposon_data["SuperFamily"].isin(ltr_elements)),
    # ["SuperFamily"],
    # ] = "Unknown_LTR_Superfam"

    # transposon_data = transposon_data[transposon_data.Order != "long_terminal_repeat"]  # drop
    # transposon_data = transposon_data[transposon_data.Order != "Maverick"]  # drop if in Order category
    # transposon_data = transposon_data[
    # transposon_data.Order != "target_site_duplication"
    # ]  # drop if in Order category
    # print(transposon_data[transposon_data["Order"] == "LINE"]["SuperFamily"].unique())

    # If TE Order is LINE and the SuperFamily is not in Wicker 2007
    # No CR1 superfamily in Fig 1
    # print(transposon_data[transposon_data["SuperFamily"] == "CR1"])
    print()
    print(transposon_data.shape)
    print(transposon_data["Order"].unique())
    print(transposon_data["SuperFamily"].unique())
    print()

    return transposon_data
