#!/usr/bin/env python3

"""
Produce OverlapData with multiprocessing.
"""

from collections import namedtuple
from functools import partial
import logging
import multiprocessing
import multiprocessing.managers
import threading
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


def _calculate_overlap_job(job):
    """Returns an OverlapResult given the job.

    Args:
        job(_OverlapJob): input data for the job
    Returns:
        OverlapResult: output data for the job
    """

    genes = GeneData.read(job.gene_path)
    transposons = TransposonData.read(job.te_path)
    overlap = OverlapWorker(job.output_dir)
    progress_cb = partial(job.progress_queue.put_nowait, 1)
    file = overlap.calculate(
        genes,
        transposons,
        job.window_range,
        job.gene_names,
        stop=job.stop_event,
        progress=progress_cb,
    )
    result = OverlapResult(filepath=file, gene_file=job.gene_path)
    return result

def _process_overlap_job(job):
    """Makes the call to calculate overlap and enqueues the result.

    Args:
        job(_OverlapJob): input data for the job
    """

    try:
        result = _calculate_overlap_job(job)
    except Exception as err:
        result = OverlapResult(exception=err, gene_file=job.gene_path)
    finally:
        job.result_queue.put(result)


class _ProgressBars:
    """Status of OverlapData production."""

    COLUMNS = 79  # MAGIC arbirtary col limit (no one has physical terminals anymore)

    def __init__(self, n_gene_names, n_chrome, result_queue, progress_queue, logger=None):
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
        self.progress_queue = progress_queue
        self.stop_event = threading.Event()
        self._chrome_thread = None
        self._gene_thread = None
        self.results = list()

    def start(self):
        """Begin execution."""

        self.results = list()
        self.stop_event.clear()
        self._chrome_thread = threading.Thread(target=self.handle_chrome)
        self._chrome_thread.start()
        self._gene_thread = threading.Thread(target=self.handle_gene)
        self._gene_thread.start()

    def stop(self):
        """End execution."""

        self.stop_event.set()
        self._chrome_thread.join()
        self._gene_thread.join()
        self.chromosomes.close()
        self.gene_names.close()

    def __enter__(self):
        """Context manager begin."""

        self.start()
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        """Context manager end."""

        self.stop()

    def handle_chrome(self):
        """Target to consume chromosome results."""

        while not self.stop_event.is_set():
            result = self._pop(self.result_queue)
            if result is not None:
                self.chromosomes.update()
                self.gene_names.refresh()
                self._log_result(result)
                self.results.append(result)

    def _log_result(self, result):
        """Log details on result if necessary."""

        if not isinstance(result, OverlapResult):
            self._logger.warning("expecting %s but got %s" %
                                 (OverlapResult.__name__, result))
            return
        if result.exception is not None:
            self._logger.error("failed to process gene '{}':  {}"
                               .format(result.gene_file, result.exception))

    def handle_gene(self):
        """Target to consume gene results."""

        while not self.stop_event.is_set():
            result = self._pop(self.progress_queue)
            if result is not None:
                self.gene_names.update()
                self.chromosomes.refresh()

    @staticmethod
    def _pop(my_queue):
        """Pop an item, return None if timeout."""

        try:
            # MAGIC NUMBER mainly delays shutdown
            result = my_queue.get(timeout=0.2)
        except queue.Empty:
            return None
        else:
            return result


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
        # TODO filter jobs yielded based on whether the file needs to be updated
        # don't recalculate overlap if already exists (add overwrite flag or just have caller delete?)
        # requires OverlapData output filename to be known, and not randomized

        self._logger = logging.getLogger(self.__class__.__name__)
        self.validate_files(data_paths, self._logger)
        self.output_dir = output_dir
        raise_if_no_dir(self.output_dir, self._logger)

        self.gene_transposon_paths = list(data_paths)
        if not self.gene_transposon_paths:
            msg = "input paths are empty: {}".format(self.gene_transposon_paths)
            raise ValueError(msg)

        self.window_range = window_range
        self.n_workers = n_workers or multiprocessing.cpu_count()
        self._stop_event = multiprocessing.Event()

        self._proc_mgr = multiprocessing.Manager()
        self._result_queue = self._proc_mgr.Queue()
        self._progress_queue = self._proc_mgr.Queue()

    def calculate_overlap(self):
        """Calculate and write OverlapData to disk."""

        self._clear()
        progress = self._new_progress_bars()
        with progress:
            with multiprocessing.Pool(processes=self.n_workers) as pool:
                jobs = self._produce_jobs()
                pool.map(_process_overlap_job, jobs)  # blocks
        return progress.results

    def _clear(self):
        """Remove items in queues and signal old jobs to stop."""

        self._stop_event.set()
        map(self._clear_queue, [self._result_queue, self._progress_queue])

    @staticmethod
    def _clear_queue(my_queue):
        """Remove all items from the queue."""

        with my_queue.mutex:
            my_queue.queue.clear()

    def _produce_jobs(self):
        """Yield _OverlapJob for each input."""

        for gene_path, te_path in self._yield_gene_trans_paths():
            # NB one can provide a partial list of gene names to process
            # and concat the results to reduce the amount of work per job

            # TODO create some sort of function to get a more appropriate name
            # for the overlap data (on a chromosome by chromosome basis)
            # TODO replace 'output_dir' w/ new filepath
            #      and remove the random file name generation in the OverlapData.from_param
            all_names = GeneData.read(gene_path).names
            job = _OverlapJob(
                    gene_path=gene_path,
                    te_path=te_path,
                    output_dir=self.output_dir,
                    window_range=self.window_range,
                    gene_names=list(all_names),
                    progress_queue=self._progress_queue,
                    result_queue=self._result_queue,
                    stop_event=None,
                    )
            yield job

    def _yield_gene_trans_paths(self):
        """Yield tuple of paths to input data files.

        Yields:
            (str, str): GeneData, TransposonData paths
        """

        for file_pair in self.gene_transposon_paths:
            gene_data_path, transposon_data_path = file_pair
            yield (gene_data_path, transposon_data_path)

    def _new_progress_bars(self):
        """An instance of the progress bars to print status."""

        prog = _ProgressBars(
                self.n_gene_names,
                self.n_chrome,
                self._result_queue,
                self._progress_queue,
                )
        return prog

    @property
    def n_chrome(self):
        """Number of chromosomes to process."""

        return len(self.gene_transposon_paths)

    @property
    def n_gene_names(self):
        """Number of named genes to process."""

        n_gene_names_ = 0
        for gene_path, te_path in self._yield_gene_trans_paths():
            gene_data = GeneData.read(gene_path)
            sub_count = sum(1 for i in gene_data.names)
            n_gene_names_ += sub_count
        return n_gene_names_

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


