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
from transposon.overlap import OverlapWorker, OverlapResult, OverlapData


_OverlapJob = namedtuple(
    "_OverlapJob",
    [
        "gene_uid",
        "gene_path",
        "te_path",
        "output_filepath",
        "window_range",
        "gene_names",
        "progress_queue",
        "result_queue",
        "stop_event",
    ],
)

# monkeypatch multiprocessing in py 3.6
# SEE https://stackoverflow.com/questions/46779860/multiprocessing-managers-and-custom-classes
# Backup original AutoProxy function
backup_autoproxy = multiprocessing.managers.AutoProxy

# Defining a new AutoProxy that handles unwanted key argument 'manager_owned'
def redefined_autoproxy(
    token,
    serializer,
    manager=None,
    authkey=None,
    exposed=None,
    incref=True,
    manager_owned=True,
):
    # Calling original AutoProxy without the unwanted key argument
    return backup_autoproxy(token, serializer, manager, authkey, exposed, incref)


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
    overlap = OverlapWorker(job.output_filepath)
    progress_cb = partial(job.progress_queue.put_nowait, 1)
    file = overlap.calculate(
        genes,
        transposons,
        job.window_range,
        job.gene_names,
        stop=job.stop_event,
        progress=progress_cb,
    )
    result = OverlapResult(
        overlap_file=file, gene_file=job.gene_path, te_file=job.te_path
    )
    return result


def _process_overlap_job(job):
    """Makes the call to calculate overlap and enqueues the result.

    Args:
        job(_OverlapJob): input data for the job
    """

    try:
        result = _calculate_overlap_job(job)
    except Exception as err:  # BUG SIGINT not caught here during testing?
        result = OverlapResult(exception=err, gene_file=job.gene_path)
        os.path.remove(result.overlap_file)
    finally:
        job.result_queue.put(result)


class _ProgressBars:
    """Status of OverlapData production."""

    COLUMNS = 79  # MAGIC arbirtary col limit (no one has physical terminals anymore)

    def __init__(
        self, n_gene_names, n_chrome, result_queue, progress_queue, logger=None
    ):
        """Initializer.

        Args:
            sub_genes(list(GeneData)): the chunk of genes for each worker
            logger(logging.Logger): logger to use, None to create
        """

        super().__init__()
        self._logger = logger or logging.getLogger(self.__class__.__name__)

        self.chromosomes = tqdm(
            total=n_chrome,
            desc="process".ljust(12, " "),
            position=0,
            ncols=self.COLUMNS,
        )
        self.gene_names = tqdm(
            total=n_gene_names,
            desc="genes".ljust(12, " "),
            position=1,
            ncols=self.COLUMNS,
        )
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
            self._logger.warning(
                "expecting %s but got %s" % (OverlapResult.__name__, result)
            )
            return
        if result.exception is not None:
            self._logger.error(
                "failed to process gene '{}':  {}".format(
                    result.gene_file, result.exception
                )
            )

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

    def __init__(self, data_paths, results_dir, window_range, n_workers=None):
        """Initializer.

        Args:
            data_paths(list(str, str)): paths to GeneData, TransposonData file pairs
            results_dir(str): top level output directory, create a sub folder here
            window_range(range(int)): window sizes for calculating overlap
            n_workers(int): process count, default cores available (or hyperthreads)
        """
        # TODO filter jobs yielded based on whether the file needs to be updated
        # don't recalculate overlap if already exists (add overwrite flag or just have caller delete?)
        # requires OverlapData output filename to be known, and not randomized

        self._logger = logging.getLogger(self.__class__.__name__)
        self.validate_files(data_paths, self._logger)
        self.output_dir = self._validate_output_dir(results_dir)

        self.gene_transposon_paths = list(data_paths)
        if not self.gene_transposon_paths:
            msg = "input paths are empty: {}".format(self.gene_transposon_paths)
            raise ValueError(msg)

        # FUTURE OPTIMIZE accept list of gene names to process, process subset per core
        self.window_range = window_range
        self.n_workers = n_workers or multiprocessing.cpu_count()
        self._stop_event = multiprocessing.Event()

        self._proc_mgr = multiprocessing.Manager()
        self._result_queue = self._proc_mgr.Queue()
        self._progress_queue = self._proc_mgr.Queue()

    def _validate_output_dir(self, directory):
        """Return a sub folder for the results, creating if necessary.

        Args:
            dir(str): top level directory for results.
        """

        raise_if_no_dir(directory, self._logger)
        sub_dir = os.path.abspath(directory)
        self._logger.info("output overlap data to %s" % sub_dir)
        return sub_dir

    def calculate_overlap(self):
        """Calculate and write OverlapData to disk.

        Returns:
            list(OverlapResult): the result from completing the job
        """

        self._clear()
        jobs = list(self._produce_jobs())
        completed, todo = self._filter_jobs(jobs)
        completed_results = [self._completed_job_2_result(job) for job in completed]
        with self._new_progress_bars(todo) as progress:
            with multiprocessing.Pool(processes=self.n_workers) as pool:
                pool.map(_process_overlap_job, todo)  # blocks execution
        return progress.results + completed_results

    def _clear(self):
        """Remove items in queues and signal old jobs to stop."""

        self._stop_event.set()
        self._result_queue = self._proc_mgr.Queue()
        self._progress_queue = self._proc_mgr.Queue()

    @staticmethod
    def _clear_queue(my_queue):
        """Remove all items from the queue."""

        with my_queue.mutex:
            my_queue.queue.clear()

    def _produce_jobs(self):
        """Yield _OverlapJob for each input."""

        for gene_path, te_path in self._yield_gene_trans_paths():
            # FUTURE OPTIMIZE process a subset of gene names per worker
            # one can provide a partial list of gene names to process
            # however, then you will need to deconflict the output filenames
            # and concat the output OverlapData files or merge them later
            gene_data = GeneData.read(gene_path)
            filepath = self._overlap_filepath(gene_data)
            yield self._overlap_job(gene_data, gene_path, te_path, filepath)

    def _filter_jobs(self, jobs):
        """Split completed jobs if output file exists."""

        completed = []
        todo = []
        for job in jobs:
            if os.path.isfile(job.output_filepath):
                self._logger.debug(
                    "skipping '%s', output exists '%s'",
                    job.gene_uid,
                    job.output_filepath,
                )
                completed.append(job)
            else:
                todo.append(job)
        return completed, todo

    def _overlap_filepath(self, gene_data):
        """Filename for the overlap temporary file."""

        # NB assuming there is only one worker for this chromosome
        filename = gene_data.chromosome_unique_id + "_overlap." + OverlapData.EXT
        filepath = os.path.join(self.output_dir, filename)
        return filepath

    def _overlap_job(self, gene_data, gene_path, te_path, filepath):
        """Create job for processing an overlap."""

        job = _OverlapJob(
            gene_uid=gene_data.chromosome_unique_id,
            gene_path=gene_path,
            te_path=te_path,
            output_filepath=filepath,
            window_range=self.window_range,
            gene_names=list(gene_data.names),  # NB process all genes
            progress_queue=self._progress_queue,
            result_queue=self._result_queue,
            stop_event=None,
        )
        return job

    def _yield_gene_trans_paths(self):
        """Yield tuple of paths to input data files.

        Yields:
            (str, str): GeneData, TransposonData paths
        """

        for file_pair in self.gene_transposon_paths:
            gene_data_path, transposon_data_path = file_pair
            yield (gene_data_path, transposon_data_path)

    def _new_progress_bars(self, jobs):
        """An instance of the progress bars to print status."""

        n_genes = sum(len(job.gene_names) for job in jobs)
        n_jobs = sum(1 for j in jobs)
        prog = _ProgressBars(n_genes, n_jobs, self._result_queue, self._progress_queue,)
        return prog

    @property
    def n_chrome(self):
        """Number of chromosomes to process."""

        return len(self.gene_transposon_paths)

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

    def _completed_job_2_result(self, job):
        """Convert an _OverlapJob to OverlapResult.

        Assumes the job has completed prior to this call.
        """

        result = OverlapResult(
            genes_processed=len(job.gene_names),
            exception=None,
            overlap_file=job.output_filepath,
            gene_file=job.gene_path,
            te_file=job.te_path,
        )
        return result
