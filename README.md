# OLD REPO
I am currently in the process of updating the pipeline to implement a number of speed and clarity improvements. I am also working on making this pipeline more universal. To view the old repository follow this [link](https://github.com/EdgerLab/TE_Density_Old)

# Purpose:
Calculate TE density for every gene in a genome, along sliding windows, and correctly allocate density values to a gene along TE Order and Superfamily classifications.

# Usage:
This code requires two input files. It requires an annotated gene file and an annotated TE file. Preferably, both files are in `gff` format, however the pipeline does work on `gtf` format, though support is not maintained for `gtf` formatted input data. TEs identities, which ought to be in the third column of the TE annotation file, should have a string name of the Order, followed by a `/`, and then the SuperFamily. For example, an entry might look like `LTR/Gypsey` or `DNA/DTH`.

`python transposon/density.py /path/to/gene_annotation_file /path/to/TE_annotation_file /path/to/config_file`

Edit the `transposon/replace_names.py` dictionary to accomodate changes to the TE naming scheme if you want to re-write some annotations or collapse some annotations together.

## Configuration File:
A configuration file is required by the density algorithm to set the beginning window size, window delta, and max window size variables. Two configuration files are located in the `config` folder. `test_run_config.ini` contains a setup for quickly testing the code. `production_run_config.ini` contains the paradigm that we use for our implementation of the code. When running `density.py` the code defaults to `test_run_config.ini`, however with the `-c` option the path to an optional configuration file may be supplied.

# Requirements
Please install Pip so that you may easily install Python packages.
Then use Pip to go over our **requirements.txt** and install the needed Python packages: `pip install -r requirements.txt`.

# To Do:
Talk with Michael about making `validate_args` its own separate file

Maybe set an option to use the EDTA TE naming system vs. Ning's naming system. Both are currently in the same file and there is the chance that they may overwrite one another. Also, the EDTA naming scheme contains some classifications that Ning specifically overwrote because she was not confident in those names, however that was for a hand-curated TE annotation file.
