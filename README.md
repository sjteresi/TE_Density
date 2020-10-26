# OLD REPO
I am currently in the process of updating the pipeline to implement a number of speed and clarity improvements. I am also working on making this pipeline more universal. To view the old repository follow this [link](https://github.com/EdgerLab/TE_Density_Old)

# Purpose:
Calculate TE density for every gene in a genome, along sliding windows, and correctly allocate density values to a gene along TE Order and Superfamily classifications.

# Usage:
This code requires two input files. It requires an annotated gene file and an annotated TE file. Preferably, both files are in `gff` format, however the pipeline does work on `gtf` format, though support is not maintained for `gtf` formatted input data. TEs identities, which ought to be in the third column of the TE annotation file, should have a string name of the Order, followed by a `/`, and then the SuperFamily. For example, an entry might look like `LTR/Gypsey` or `DNA/DTH`.

Edit the `transposon/replace_names.py` dictionary to accomodate changes to the TE naming scheme if you want to re-write some annotations or collapse some annotations together.

# Command-line Arguments:
* `--reset_h5` Defaults to `False`. When `True`, resets the cache of wrapped GeneData and TransposonData input files. If the `--reset_h5` flag is provided as an argument to `density.py`, then the boolean `True` is supplied.
* `--contig_del` Defaults to `True`. When `True`, deletes any entries from both the gene annotation file and transposon annotation (GFF format) file that are labeled as "Contig" (case in-sensitive) for column 1. If the `--contig_del` flag is provided as an argument, then the boolean `False` is supplied.








## Configuration File:
A configuration file is required by the density algorithm to set the beginning window size, window delta, and max window size variables. Two configuration files are located in the `config` folder. `test_run_config.ini` contains a setup for quickly testing the code. `production_run_config.ini` contains the paradigm that we use for our implementation of the code. When running `density.py` the code defaults to `test_run_config.ini`, however with the `-c` option the path to an optional configuration file may be supplied.

## Temporary Overlap Data Location:
Gene and TE overlap values are an intermediate calculation in this pipeline. However the overlap values can become quite large, so they must be stored in a temporary location. The script defaults to storing the compressed overlap files at `/tmp` however with the `-s` option the user can specify a different directory to output the overlap files.

# Requirements
Please install Pip so that you may easily install Python packages.
Then use Pip to go over our **requirements.txt** and install the needed Python packages: `pip install -r requirements.txt`.


# Control-Flow:
1. Import and clean the gene and TE annotations.
	* Revise the TE annotation to remove nested TEs
	* Wrap the annotations using the `Gene_Data` and `Transposon_Data` classes.
	* Save te datastructures to disk in the form of altered annotation files and cache the wrapped data in H5 form per chromosome to disk.
2. Read the config parameters and feed them to the OverlapManager which is then started.
	* Call .calculate_overlap on the OverlapManager instance. Starts a progress bar and creates a pool for jobs.
	* `Overlap.calculate` is then run, results are stored to a file. `overlap.calculate` iterates over the genes, windows, and TEs, calculating the intra, left, and right overlap matrices. A file path is returned.
3. Jobs that were put on the queue to OverlapResult are then put on a result queue.
4. This is where I start to get confused. I do not know how the result queue gets finished? I assume it ends gracefully with the context manager. From `process.py` it seems get returned a `_ProgressBars` class object and its attribute `.result` which is a list of results.
5. Tie this list of results together with the normalization matrices? Currently, the normalization matrix doesn't have a lot of accessor functions, just the division functions, which shouldn't be too hard to interface with the overlap. We will just need to make sure the dataframes are in the same order and then they can be divided.
