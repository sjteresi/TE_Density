#!/usr/bin/env python3

"""
Produce OverlapData with multiprocessing.
"""

from collections import namedtuple
from functools import partial
import logging
import multiprocessing
import multiprocessing.managers
from threading import Thread
import queue
import os

from tqdm import tqdm

from transposon import FILE_DNE, raise_if_no_file, raise_if_no_dir
from transposon.worker import WorkerProcess, Sentinel
from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.overlap import OverlapWorker, OverlapResult


_OverlapJob = namedtuple(
    "_OverlapJob",
    ['gene_path', 'te_path', 'output_dir', 'window_range', 'gene_names', 'progress_queue', 'result_queue', 'stop_event']
)


# monkeypatch multiprocessing in py 3.6
# SEE https://stackoverflow.com/questions/46779860/multiprocessing-managers-and-custom-classes
# Backup original AutoProxy function
backup_autoproxy = multiprocessing.managers.AutoProxy

# Defining a new AutoProxy that handles unwanted key argument 'manager_owned'
def redefined_autoproxy(token, serializer, manager=None, authkey=None,
          exposed=None, incref=True, manager_owned=True):
    # Calling original AutoProxy without the unwanted key argument
    return backup_autoproxy(token, serializer, manager, authkey,
                     exposed, incref)

# Updating AutoProxy definition in multiprocessing.managers package
multiprocessing.managers.AutoProxy = redefined_autoproxy


def _overlap_job(job, update_cb=None):

    genes = GeneData.read(job.gene_path)
    transposons = TransposonData.read(job.te_path)
    overlap = OverlapWorker(job.output_dir)
    progress_cb = update_cb or partial(job.progress_queue.put_nowait, 1)
    file = overlap.calculate(
        genes,
        transposons,
        job.window_range,
        job.gene_names,
        stop=job.stop_event,
        progress=progress_cb,
    )
    return file

def process_overlap_job(job):

    try:
        file = _overlap_job(job)
    except Exception as err:
        result = OverlapResult(exception=err)
    else:
        result = OverlapResult(filepath=file)
    finally:
        job.result_queue.put(result)

class OverlapProcess(WorkerProcess):
    """Computes overlap in a process."""

    def execute_job(self, job):
        """Calculate the overlap."""

        try:
            self.logger.info("start calculation")
            file = process_overlap_job(job)
        except Exception as err:
            result = OverlapResult(exception=err)
        else:
            result = OverlapResult(filepath=file)
        finally:
            return result


class _ProgressBars(Thread):
    """Status of OverlapData production."""

    COLUMNS = 79  # MAGIC arbirtary col limit (no one has physical terminals anymore)

    def __init__(self, n_gene_names, n_chrome, result_queue, stop_event, logger=None):
        """Initializer.

        Args:
            sub_genes(list(GeneData)): the chunk of genes for each worker
            logger(logging.Logger): logger to use, None to create
        """

        super().__init__()
        self._logger = logger or logging.getLogger(self.__class__.__name__)
        self.chromosomes = tqdm(
            total=n_chrome, desc="chromosome  ", position=0, ncols=self.COLUMNS)
        self.gene_names = tqdm(
            total=n_gene_names, desc="  genes     ", position=1, ncols=self.COLUMNS)
        self.result_queue = result_queue
        self.stop = stop_event

    def _handle_result(self, result):

        self.gene_names.update(result.genes_processed)
        if result.filepath:
            self.chromosomes.update(1)
        else:
            self.chromosomes.refresh()

        if not result.exception:
            self._logger.error("failed to process gene '{}':  {}"
                               .format(result.genome_id, result.exception))

    def _pop_result(self):
        """Return an OverlapResult, or None if timed out."""

        try:
            # MAGIC NUMBER mainly delays shutdown
            result = self.result_queue.get(timeout=0.2)
        except queue.Empty:
            return None
        else:
            return result

    def run(self):
        """Target to execute."""

        while not self.stop.is_set():
            result = self._pop_result()
            if result:
                self._handle_result(result)
        self.chromosomes.close()
        self.gene_names.close()


class OverlapManager:
    """Orchestrate multiple OverlapWorkers."""

    def __init__(self,
                 data_paths,
                 output_dir,
                 window_range,
                 gene_names=None,
                 n_workers=None
                 ):
        """Initializer.

        Args:
            data_paths(list(str, str)): paths to GeneData, TransposonData file pairs
            window_range(range(int)): window sizes for calculating overlap
            gene_names(list(str)): gene names to calculate overlap for, default all
            n_workers(int): process count
        """


        self._logger = logging.getLogger(self.__class__.__name__)
        self.validate_files(data_paths, self._logger)
        self.output_dir = output_dir
        raise_if_no_dir(self.output_dir, self._logger)

        self.gene_transposon_paths = list(data_paths)
        if not self.gene_transposon_paths:
            msg = "input paths are empty: {}".format(self.gene_transposon_paths)
            raise ValueError(msg)

        self.window_range = window_range
        self._proc_mgr = multiprocessing.Manager()
        self.n_workers = n_workers or multiprocessing.cpu_count()
        self.n_workers = 1  # TODO remove when done testing
        self._stop_event = multiprocessing.Event()
        self._workers = None
        self._result_queue = None  # multiprocessing.Queue
        self._job_queue = None  # multiprocessing.Queue
        self._progress = None  # _ProgressBars

    def __enter__(self):
        """Acquire resources and begin processing."""

        self.start()
        return self

    def start(self):
        """Begin processing jobs."""

        self.stop()
        if self._workers is not None:
            msg = ("worker list not empty, cannot restart: {}".format(self._workers))
            self._logger.critical(msg)
            raise RuntimeError(msg)
        self._stop_event.clear()
        self._job_queue = self._proc_mgr.Queue()
        self._result_queue = self._proc_mgr.Queue()
        self._fill_jobs()
        self._progress = self._new_progress_bars()
        self._progress.start()
        self._workers = [self._new_worker() for i in range(self.n_workers)]
        for i, w in enumerate(self._workers):
            self._logger.debug("starting worker %d" % i)
            w.start()

    def __exit__(self, exc_type, exc_val, traceback):
        """End and release resources."""

        self.stop()

    def stop(self):
        """End processing and join."""

        self._signal_end()
        self.join_workers()
        if self._job_queue:
            self._job_queue = None
        if self._result_queue:
            self._result_queue = None
        if self._progress:
            self._progress.join()
            self._progress = None

    def join_workers(self):
        """Join worker processes."""

        if self._workers:
            for w in self._workers:
                w.join()
        self._workers = None

    def _new_worker(self):
        """Return a initialized OverlapProcess."""

        return OverlapProcess(self._job_queue, self._result_queue, self._stop_event)

    def _produce_jobs(self):
        """Yield _OverlapJob for each input."""

        for gene_path, te_path in self._yield_gene_trans_paths():
            # NB one can reduce the names used per worker and concat later
            all_names = GeneData.read(gene_path).names
            job = _OverlapJob(
                    gene_path=gene_path,
                    te_path=te_path,
                    output_dir=self.output_dir,
                    window_range=self.window_range,
                    gene_names=list(all_names),
                    progress_queue=self._result_queue,
                    stop_event=None,
                    )
            yield job

    def _fill_jobs(self):

        self._logger.debug("filling jobs...")
        for i, job in enumerate(self._produce_jobs()):
            self._job_queue.put(job)
        for i in range(self.n_workers):
            self._job_queue.put_nowait(Sentinel())

    def _yield_gene_trans_paths(self):

        for file_pair in self.gene_transposon_paths:
            gene_data_path, transposon_data_path = file_pair
            yield (gene_data_path, transposon_data_path)

    def _new_progress_bars(self):

        prog = _ProgressBars(
                self.n_gene_names,
                self.n_chrome,
                self._result_queue,
                self._stop_event
                )
        return prog

    @property
    def n_chrome(self):

        return len(self.gene_transposon_paths)

    @property
    def n_gene_names(self):

        n_gene_names_ = 0
        for gene_path, te_path in self._yield_gene_trans_paths():
            gene_data = GeneData.read(gene_path)
            sub_count = sum(1 for i in gene_data.names)
            n_gene_names_ += sub_count
        return n_gene_names_

    def _signal_end(self):
        """Set the stop event and enqueue sentinels for each worker."""

        self._stop_event.set()
        if self._job_queue:
            [self._job_queue.put_nowait(Sentinel()) for i in range(self.n_workers)]

    @staticmethod
    def validate_files(file_pairs, logger):
        """Raise FileNotFoundError if the files do not exist.

        Args:
            file_pairs(list(str, str)): filepath pairs
            logger(logging.Logger): logging instance
        """

        genes = [g for g, t in file_pairs]
        transposons = [t for g, t in file_pairs]
        paths = genes + transposons
        validate = partial(raise_if_no_file, logger=logger)
        map(validate, paths)

