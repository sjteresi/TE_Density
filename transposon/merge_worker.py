#!/usr/bin/env python

"""
Orchestrates the summation of many OverlapData into many MergeData.
"""

__author__ = "Michael Teresi"

from collections import namedtuple, defaultdict
import logging
import os
from functools import partial
import multiprocessing
import h5py
import numpy as np

import transposon
from transposon import FILE_DNE, MAX_SYSTEM_RAM_GB
from transposon.overlap import OverlapData
from transposon.transposon_data import TransposonData

# FUTURE multiprocessing pool, for the initializer
# setup_merge_proc(overlap_paths)
#    global merge_element
#    merge_element = _MergeElement(overlap_paths)

def _process_merge_element(merge_config):
    """Helper function to execute a merge."""

    merge_element = _MergeElement(merge_config)
    return merge_element.process()

class _MergeElement():
    """Combines multiple OverlapData for one sub chromosome.

    Many OverlapData results are combined into one MergeData.
    This produces results for both superfamilies & order, SEE MergeData.
    This occurs after calculating overlap and prior to normalization.
    The OverlapData must be for the same sub-chromosome, e.g. Fvb1-1.
    """

    class Configuration:
        """Parameters for merging."""

        def __init__(self, overlap_paths, te_path, output_dir, log_callback=None, ram=2):
            """Initializer.

            Args:
                overlap_paths(list(str)): paths to OverlapData files
                te_path(str): path to TransposonData
                output_dir(str): path to write results
                log_callback(callable): configure local logger
                ram(float): gigabytes of RAM for output file
            """
            logging.info("paths: {}".format(overlap_paths))
            self.overlap_paths = overlap_paths
            self.te_path = str(te_path)
            self.output_dir = str(output_dir)
            self.log_callback = log_callback
            self.ram = float(ram)

            self.validate()

        def validate(self):
            """Raise if bad inputs."""

            if not self.overlap_paths:
                raise ValueError("overlap paths are empty")
            for overlap_path in self.overlap_paths:
                if not os.path.exists(overlap_path):
                    raise FILE_DNE(overlap_path)
            if not os.path.exists(self.te_path):
                raise FILE_DNE(self.te_path)
            if not os.path.exists(self.output_dir):
                raise FILE_DNE(self.output_dir)
            if self.ram <= 0 or self.ram > MAX_SYSTEM_RAM_GB:
                raise ValueError("requested GB {} > system RAM of {}"
                                 .format(self._ram, system_gb))

    def __init__(self, config):
        """Initializer.

        """

        self._logger = logging.getLogger(self.__class__.__name__)
        self._logger.info("overlap files: {}".format(config.overlap_paths))
        self._overlaps = [OverlapData.from_file(o) for o in config.overlap_paths]
        self.acquire()
        self.validate_overlaps(self._overlaps)
        self._windows = self._scrape_windows(self._overlaps)
        self._gene_names = self._scrape_genes(self._overlaps)
        self.release()
        self._te_data = TransposonData.read(config.te_path)
        self._output_dir = os.path.abspath(config.output_dir)
        self._ram = float(config.ram) if config.ram > 0 else 1  # gigabytes

        if self._max_workers <= 0:
            raise ValueError("worker count '{}' <= 0".format(self._max_workers))

    def acquire(self):
        [overlap.start() for overlap in self._overlaps]
        # TODO acquire output file

    def release(self):
        [overlap.stop() for overlap in self._overlaps]
        # TODO release output file

    def process(self):
        self.acquire()
        merge_data = MergeData.from_param(
            self._te_data,
            self._gene_names,
            self._windows,
            self._output_dir,
            self._ram
        )
        output_path = None
        with merge_data as merging:
            output_path = merging.filepath;
            for o in self._overlaps:
                merging.sum(o)

        self.release()
        return output_path

    @staticmethod
    def validate_overlaps(overlaps):
        """Raise if overlap list cannot be used."""

        if not overlaps:
            raise ValueError("input overlap list is empty")

        ids = set()
        for overlap in overlaps:
            ids.add(overlap.chromosome_id)
            if overlap.chromosome_id is None:
                msg = ("OverlapData chromosome ID is None {}"
                       .format(overlap.filepath))
                raise ValueError(msg)
        if len(ids) != 1:
            msg = ("OverlapData list has multiple chromosome IDs: {} for {}"
                   .format(ids, [o.filepath for o in overlaps]))
            raise ValueError(msg)

    @staticmethod
    def _scrape_windows(overlaps):
        """Find the set of windows the overlap files contain."""

        window_set = set()
        for overlap in overlaps:  # requires an open overlap instance
            if overlap.windows is None:
                raise ValueError("missing window list in {}"
                                 .format(overlap.filepath))
            for w in overlap.windows:
                window_set.add(w)
        return list(window_set)

    @staticmethod
    def _scrape_genes(overlaps):
        """Find the set of gene names the overlap files contain."""

        gene_set = set()
        for overlap in overlaps:  # requires an open overlap instance
            if overlap.gene_names is None:
                raise ValueError("missing gene names in {}"
                                 .format(overlap.filepath))
            for g in overlap.gene_names:
                gene_set.add(g)
        return list(gene_set)


class MergeManager():
    """Orchestrates merging multiple OverlapData."""

    def __init__(self, overlaps, transposons, genome, output_dir, ram, log_level, workers=None):
        """Initialzer.

        Args:
            overlaps (iterable(str)): filepaths to OverlapData.
            transposons (str): path to TransposonData file.
            genes (str): path to GeneData file.
            windows (iterable(int)): input window values.
            output_dir (str): path to output results.
        """

        self._logger = logging.getLogger(self.__class__.__name__)
        self._logger.setLevel(log_level)
        self._overlaps = [OverlapData.from_file(str(o)) for o in overlaps]
        self._te_path = str(transposons)
        self._output_dir = os.path.abspath(output_dir)
        self._max_workers = workers or os.cpu_count()
        self._pool = multiprocessing.Pool(self._max_workers)
        self._ram = int(ram)
        self._proc_manager = multiprocessing.Manager()
        self._log_queue = self._proc_manager.Queue(-1)

        [o.start() for o in self._overlaps]  # is `map` more pythonic?
        self._chromosomes = set(o.chromosome_id for o in self._overlaps)
        self._chrome2files = self._chrome_ids_to_path(self._overlaps)
        self._chrome2windows = self._chrome_ids_to_windows(self._overlaps)
        [o.stop() for o in self._overlaps]
        self._logger.info("merging {} overlaps files".format(len(self._overlaps)))
        self._logger.info("merging {} sub chromosomes".format(len(self._chromosomes)))
        self._merge_configs = self._index_overlaps(self._chromosomes, self._chrome2files)

        max_ram_per_worker = MAX_SYSTEM_RAM_GB / self._max_workers
        if self._ram > max_ram_per_worker:
            raise ValueError("ram per worker {} is > {}"
                             .format(self._ram, max_ram_per_worker))

    def merge(self):

        try:
            output_paths = self._pool.map(
                _process_merge_element, self._merge_configs)
        except KeyboardInterrupt:
            return []
        else:
            return output_paths


    def __enter__(self):
        """Context manager start."""

        self._stop()
        self._start()
        return self

    def __exit__(self, exc_type, exc_val, exc_traceback):
        """Context manager stop."""

        self._stop()

    def _start(self):

        pass

    def _stop(self):

        pass

    def _index_overlaps(self, chromosomes, chrome2file):

        merge_configs = []
        for chrome in chromosomes:
            cfg = _MergeElement.Configuration(
                chrome2file[chrome],
                self._te_path,
                self._output_dir,
                ram=self._ram
            )
            merge_configs.append(cfg)
        return merge_configs

    @staticmethod
    def _chrome_ids_to_path(overlaps):
        """Scrape the chromosome IDs to paths."""

        chrome_2_path = defaultdict(list)
        for overlap in overlaps:  # requires an open overlap instance
            if overlap.chromosome_id is None:
                raise ValueError("OverlapData  has no path, is the file active?:  {}"
                                 .format(overlap))
            chrome_2_path[overlap.chromosome_id].append(overlap.filepath)
        return chrome_2_path

    @staticmethod
    def _chrome_ids_to_windows(overlaps):
        """Scrape the chromosome IDs to windows"""

        chrome_2_windows = defaultdict(set)
        for overlap in overlaps:  # requires an open overlap instance
            if overlap.windows is None:
                raise ValueError("missing window list in {}".format(overlap.filepath))
            current_set = chrome_2_windows[overlap.chromosome_id]
            [current_set.add(win) for win in overlap.windows]
        return chrome_2_windows

    # scrape the windows from the overlap files, pass to MergeData
    # scrape genes from the overlap files, pass to MergeData


if __name__ == "__main__":
    import argparse
    import glob
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser(description='temporary integration test')
    parser.add_argument('input_dir')
    args = parser.parse_args()
    cache_dir = os.path.abspath(args.input_dir)
    data_dir = os.path.abspath(os.path.join(cache_dir, "../../TE_Density/filtered_input_data"))
    logging.info("gene data? {}".format(data_dir))
    o_files = glob.glob(cache_dir + "/*.h5")
    # BUG the TSV files are not what we want here
    # we need the GeneData files, but TSV is the unprocessed data
    te_file = glob.glob(data_dir + "/*gtf_data.tsv")[0]  # MAGIC there can be only one
    gene_file = glob.glob(data_dir + "/*gff_data.tsv")[0]  # MAGIC there can be only one
    logging.info("transposon file  {}".format(te_file))
    logging.info("gene_file        {}".format(gene_file))
    manager = MergeManager(o_files, te_file, gene_file, '/tmp/', 2, logging.DEBUG)


    with manager as boss:
        logging.info("merging...")
        merged_files = boss.merge()
        logging.info("complete: {}".format(merged_files))
