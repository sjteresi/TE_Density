Calculate Transposable Element (TE) density for a genome.

Calculates TE density:

1. For every gene in a genome
2. Along sliding windows
3. With respect to TE Order and TE Superfamily classifications

# Contact Information:
|               | Name           | GitHub                                                  | Email              |
|---------------|----------------|---------------------------------------------------------|--------------------|
| Project Lead: | Scott Teresi   | [Personal GitHub](https://github.com/huckleberry-hound) | <teresisc@msu.edu> |
| Co-Author     | Michael Teresi | [Personal GitHub](https://github.com/teresi)            |                    |
| PI:           | Patrick Edger  | [Lab GitHub](https://github.com/EdgerLab)               | <edgerpat@msu.edu> |


# Usage:

```
./process_genome.py --help
```

This code requires two input files. It requires an annotated gene file and an annotated TE file. Preferably, both files are in `gff` format, however the pipeline does work on `gtf` format, though support is not maintained for `gtf` formatted input data. TEs identities, which ought to be in the third column of the TE annotation file, should have a string name of the Order, followed by a `/`, and then the SuperFamily. For example, an entry might look like `LTR/Gypsy` or `DNA/DTH`.

Edit the `transposon/replace_names.py` dictionary to accomodate changes to the TE naming scheme if you want to re-write some annotations or collapse some annotations together.


# Command-line Arguments:
* `--reset_h5` Defaults to `False`. When `True`, resets the cache of wrapped GeneData and TransposonData input files. If the `--reset_h5` flag is provided as an argument to `density.py`, then the boolean `True` is supplied.
* `--contig_del` Defaults to `True`. When `True`, deletes any entries from both the gene annotation file and transposon annotation (GFF format) file that are labeled as "Contig" (case in-sensitive) for column 1. If the `--contig_del` flag is provided as an argument, then the boolean `False` is supplied.


## Configuration File:
A configuration file is required by the density algorithm to set the beginning window size, window delta, and max window size variables. Two configuration files are located in the `config` folder. `test_run_config.ini` contains a setup for quickly testing the code. `production_run_config.ini` contains the paradigm that we use for our implementation of the code. When running `density.py` the code defaults to `test_run_config.ini`, however with the `-c` option the path to an optional configuration file may be supplied.


## Temporary Overlap Data Location:
Gene and TE overlap values are an intermediate calculation in this pipeline. However the overlap values can become quite large, so they must be stored in a temporary location. The script defaults to storing the compressed overlap files at `/tmp` however with the `-s` option the user can specify a different directory to output the overlap files.

# Requirements
SEE [requirements.txt](requirements.txt)
 
SEE [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/command_ref.html)
```bash
mkvirtualenv <mynewenv>  # or your preferred environment
pip install -r requirements.txt
```
