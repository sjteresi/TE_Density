import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os


def graph_barplot_density_differences(
    values,
    te_type,
    window_val,
    direction,
    number_of_zeros,
    output_dir,
    logger,
    display=False,
    align="left",
):
    """
    Plot a histogram of TE density differences between syntelog pairs

    Args:
        values (list): A list of values representing the TE density differences
        between syntelog pairs

        te_type (str): String representing the TE type being plotted

        window_val (int): Integer representing the current window of which the
        data is being plotted

        direction (str): string representing whether or not the graphs are
        coming from upstream or downstream TE density data

        number_of_zeros ():

        logger (logging.Logger): Object to log information to

        display (boolean): Defaults to False, if True shows the plot upon
        generation with the plt.show() command
    """

    # MAGIC, the bins to group density values for the histogram AND the values
    # for the xticks on the xaxis
    tick_bins = [
        -1.0,
        -0.9,
        -0.8,
        -0.7,
        -0.6,
        -0.5,
        -0.4,
        -0.3,
        -0.2,
        -0.1,
        0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
    ]

    plt.figure(figsize=(8, 6))
    n, bins, patches = plt.hist(
        values, bins=tick_bins, facecolor="blue", ec="black", alpha=0.5
    )
    plt.rcParams["xtick.labelsize"] = 7  # MAGIC set size of axis ticks
    plt.ylabel("Number of Genes")
    plt.xlabel("Difference in TE Density Values")
    plt.title("O. glaberrima vs O. sativa")  # MAGIC genome name order here
    N = mpatches.Patch(
        label="Total Plotted Genes: %s \nTE type: %s \nWindow: %s \nDirection: %s \nNo. 0 Differences: %s"
        % (len(values), te_type, window_val, direction, str(number_of_zeros))
    )
    plt.xticks(tick_bins)
    plt.legend(handles=[N])
    path = os.path.join(
        output_dir,
        (te_type + "_" + str(window_val) + "_" + direction + "_DensityDifferences.png"),
    )
    logger.info("Saving graph to: %s" % path)
    plt.savefig(path)
    if display:
        plt.show()
    plt.close()
