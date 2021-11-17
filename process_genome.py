#!/usr/bin/env python

"""
Calculate transposable element density.
"""

__author__ = "Scott Teresi, Michael Teresi"

import argparse
import os
import cProfile, pstats, io
import logging
import coloredlogs
import numpy as np
import pandas as pd
from tqdm import tqdm
from configparser import ConfigParser
import sys
import time

from transposon import FILE_DNE, set_numexpr_threads
from transposon import raise_if_no_file, raise_if_no_dir
from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.preprocess import PreProcessor
from transposon.overlap_manager import OverlapManager
from transposon.overlap import OverlapData
from transposon.merge_data import MergeData


def validate_args(args, logger):
    """Raise if an input argument is invalid."""

    raise_if_no_file(
        args.genes_input_file,
        logger=logger,
        msg_fmt="arg 'genes_input_file' not a file: %s",

    )
    raise_if_no_file(
        args.tes_input_file,
        logger=logger,
        msg_fmt="arg 'tes_input_file' not a file: %s",
    )
    raise_if_no_dir(
        args.output_dir,
        logger=logger,
        msg_fmt="""arg 'output_dir' not a dir: %s \n Please make sure this
        directory exists"""
    )


def parse_algorithm_config(config_path):
    """Return parameters for running density calculations."""

    raise_if_no_file(config_path)
    parser = ConfigParser()
    parser.read(config_path)
    window_start = parser.getint("density_parameters", "first_window_size")
    window_step = parser.getint("density_parameters", "window_delta")
    window_stop = parser.getint("density_parameters", "last_window_size")
    alg_param = {
        "window_range": range(window_start, window_stop + 1, window_step)
    }  # MAGIC we want inclusive
    # range of windows
    return alg_param


# FUTURE move to class
from collections import namedtuple
from functools import partial
from multiprocessing import Pool, Manager
from threading import Thread, Event
from queue import Empty


MergeJob = namedtuple(
    "MergeJob",
    ["overlap_file", "te_file", "gene_file", "windows", "output_dir", "progress_bar"],
)


def job_2_merge_and_overlap(job):
    """Return an instance of MergeData, OverlapData given the job."""

    transposons = TransposonData.read(job.te_file)
    windows = list(job.windows)
    output_dir = str(job.output_dir)
    gene_data = GeneData.read(job.gene_file)
    merge_data = MergeData.from_param(transposons, gene_data, windows, output_dir)
    overlap_data = OverlapData.from_file(job.overlap_file)
    return merge_data, overlap_data


def calc_merge_number_operations(job):
    """Number of expected progress updates for processing a MergeJob.

    NOTE this is a candidate for refactoring alongside the way density is calculated.

    Args:
        job(MergeJob): the job
    """

    merge_data, overlap_data = job_2_merge_and_overlap(job)
    with merge_data as merge_output:
        with overlap_data as overlap_input:
            return merge_data.n_updates(overlap_data)


def calc_merge(job):
    """Target for a process to calculate density given a job."""

    merge_data, overlap_data = job_2_merge_and_overlap(job)
    gene_data = GeneData.read(job.gene_file)
    with merge_data as merge_output:
        with overlap_data as overlap_input:
            merge_output.sum(overlap_input, gene_data, job.progress_bar)


class MergeProgress:
    """Show progress of the merge."""

    def __init__(self, queue, progress):

        self.stop_event = Event()
        self.queue = queue
        self.progress = progress
        self._thread = None

    def __enter__(self):

        self.stop_event.clear()
        if self._thread is not None:
            self._thread.join()
        self._thread = Thread(target=self.exec)
        self._thread.start()
        return self

    def __exit__(self, type, val, trace):

        self.stop_event.set()

    def exec(self):
        while not self.stop_event.is_set():
            try:
                result = self.queue.get(timeout=0.2)
            except Empty:
                continue

            self.progress.update()


def result_to_job(result, windows, output_dir, pbar_callback):
    """Convert an overlap result to a merge job."""
    job = MergeJob(
        str(result.overlap_file),
        str(result.te_file),
        str(result.gene_file),
        list(windows),
        str(output_dir),
        pbar_callback,
    )
    return job


if __name__ == "__main__":
    """Command line interface to calculate density."""

    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    output_default = os.path.join(dir_main, "..", "TE_Data")
    parser = argparse.ArgumentParser(description="calculate TE density")

    parser.add_argument("genes_input_file", type=str, help="parent path of gene file")

    parser.add_argument(
        "tes_input_file", type=str, help="parent path of transposon file"
    )

    parser.add_argument("genome_id", type=str, help="name of genome")

    parser.add_argument(
        "--num_threads",
        "-n",
        default=None,
        help="number of threads for code, defaults to machine max",
    )

    parser.add_argument(
        "--config_file",
        "-c",
        type=str,
        default=os.path.join(path_main, "..", "config/test_run_config.ini"),
        help="parent path of config file",
    )

    parser.add_argument(
        "--reset_h5",
        action="store_true",
        help="Rewrite h5 intermediate files for gene & TEs",
    )

    parser.add_argument(
        "--revise_anno",
        action="store_true",
        help="""Forces the
                        recreation of a revised TE annotation file. Desirable if
                        you have previously created a revised TE annotation but
                        you want the pipeline to create a new one from scratch
                        and overwrite the cache. This is especially useful if
                        you have modified the input TE annotation but have not
                        changed the filename.""",
    )

    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default=output_default,
        help="parent directory to output results",
    )

    parser.add_argument(
        "--single_process",
        action="store_true",
        help="""Run without multiprocessing; useful for profiling.""",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.genes_input_file = os.path.abspath(args.genes_input_file)
    args.tes_input_file = os.path.abspath(args.tes_input_file)
    args.config_file = os.path.abspath(args.config_file)
    args.output_dir = os.path.abspath(args.output_dir)
    if args.num_threads is not None:
        args.num_threads = int(args.num_threads)

    filtered_input_data_loc = os.path.abspath(
        os.path.join(args.output_dir, "filtered_input_data")
    )
    input_h5_cache_loc = os.path.abspath(
        os.path.join(args.output_dir, filtered_input_data_loc, "input_h5_cache")
    )

    revised_input_data_loc = os.path.abspath(
        os.path.join(args.output_dir, filtered_input_data_loc, "revised_input_data")
    )
    tmp_overlap_loc = os.path.abspath(os.path.join(args.output_dir, "tmp", "overlap"))

    # NOTE make directories for intermediate and final output data
    os.makedirs(filtered_input_data_loc, exist_ok=True)
    os.makedirs(input_h5_cache_loc, exist_ok=True)
    os.makedirs(revised_input_data_loc, exist_ok=True)
    os.makedirs(tmp_overlap_loc, exist_ok=True)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    for argname, argval in vars(args).items():
        logger.debug("%-18s: %s" % (argname, argval))
    validate_args(args, logger)
    alg_parameters = parse_algorithm_config(args.config_file)

    set_numexpr_threads(
        args.num_threads
    )  # prevents an unenecessary log call from numexpr

    logger.info("preprocessing...")
    preprocessor = PreProcessor(
        args.genes_input_file,
        args.tes_input_file,
        filtered_input_data_loc,
        input_h5_cache_loc,
        revised_input_data_loc,
        args.reset_h5,
        args.genome_id,
        args.revise_anno,
    )
    preprocessor.process()
    n_data_files = sum(1 for _ in preprocessor.data_filepaths())
    rel_preproc = os.path.relpath(input_h5_cache_loc)
    logger.info("preprocessed %d files to %s" % (n_data_files, rel_preproc))
    logger.info("preprocessing... complete")

    logger.info("process overlap...")

    gene_te_filepaths = list(preprocessor.data_filepaths())
    overlap_mgr = OverlapManager(
        gene_te_filepaths, tmp_overlap_loc, alg_parameters["window_range"]
    )
    overlap_results = overlap_mgr.calculate_overlap()
    logger.info("processed %d overlap jobs" % len(overlap_results))
    logger.info("process overlap... complete")

    logger.info("process density")

    win = alg_parameters["window_range"]
    pbar_update_mgr = Manager()
    pbar_update_queue = pbar_update_mgr.Queue()
    pbar_update = partial(pbar_update_queue.put_nowait, None)
    jobs = [
        result_to_job(res, win, args.output_dir, pbar_update) for res in overlap_results
    ]
    n_subsets = sum(calc_merge_number_operations(job) for job in jobs)
    pbar_subsets = tqdm(
        total=n_subsets,
        desc="subsets",
        position=0,
        ncols=80,
    )
    with MergeProgress(pbar_update_queue, pbar_subsets) as my_progress:
        if not args.single_process:
            with Pool(processes=args.num_threads) as my_pool:
                my_pool.map(calc_merge, jobs)
        else:
            for job in jobs:
                try:
                    pr = cProfile.Profile()
                    pr.enable()
                    calc_merge(job)
                except KeyboardInterrupt as keybr:
                    pr.disable()
                    stream = io.StringIO()
                    sortby = "cumulative"
                    ps = pstats.Stats(pr, stream=stream).sort_stats(sortby)
                    ps.print_stats(0.1)  # MAGIC percent to print
                    print(stream.getvalue())
                    raise keybr

    logger.info("process density... complete")
