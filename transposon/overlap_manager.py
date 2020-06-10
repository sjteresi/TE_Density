#!/usr/bin/env python3

"""
Produce OverlapData with multiprocessing.
"""

from collections import namedtuple
from functools import partial
import logging
import multiprocessing
import os

from tqdm import tqdm

from transposon import FILE_DNE
from transposon.worker import WorkerProcess, Sentinel
from transposon.gene_data import GeneData
from transposon.transposon_data import TransposonData
from transposon.overlap import OverlapWorker, OverlapResult


_OverlapJob = namedtuple(
    "_OverlapJob",
    ['gene_path', 'te_path', 'output_dir', 'window_slice', 'gene_names', 'progress_queue']
)


class OverlapProcess(WorkerProcess):
    """Computes overlap in a process."""

    def execute_job(self, job):
        """Calculate the overlap."""

        try:
            genes = GeneData.read(job.gene_path)
            transposons = TransposonData.read(job.te_path)
            overlap = OverlapWorker(job.output_dir)
            # TODO refactor the progress update into OverlapWorker
            progress_cb = partial(job.progress_queue.put_nowait, 1)
            file = overlap.calculate(
                genes,
                transposons,
                job.window_slice,
                job.gene_names,
                stop=self._stop_event,
                progress=progress_cb,
            )
        except Exception as err:
            result = OverlapResult(exception=err)
        else:
            result = OverlapResult(filepath=file)
        finally:
            return result


class _ProgressBars:
    """Status of OverlapData production."""

    COLUMNS = 79  # MAGIC arbirtary col limit (no one has physical terminals anymore)

    def __init__(self, sub_genes, logger=None):
        """Initializer.

        Args:
            sub_genes(list(GeneData)): the chunk of genes for each worker
            logger(logging.Logger): logger to use, None to create
        """

        super().__init__()
        self._logger = logger or logging.getLogger(self.__name__.__class__)
        self.chromosomes = tqdm(
            total=len(sub_genes), desc="chromosome  ", position=0, ncols=self.COLUMNS)
        n_gene_names = sum(sum(1 for _ in gene.names) for gene in sub_genes)
        self.gene_names = tqdm(
            total=n_gene_names, desc="  genes     ", position=1, ncols=self.COLUMNS)

    def handle_result(self, result):

        self.gene_names.update(result.genes_processed)
        if result.filepath:
            self.chromosomes.update(1)
        else:
            self.chromosomes.refresh()

        if not result.exception:
            self._logger.error("failed to process gene '{}':  {}"
                               .format(result.genome_id, result.exception))


class OverlapManager:
    """Orchestrate multiple OverlapWorkers."""

    def __init__(self,
                 sub_genes,
                 sub_tes,
                 output_dir,
                 window_slice,
                 gene_names=None,
                 n_workers=None
                 ):
        """Initializer.

        Args:
            sub_genes(list(str)): paths to GeneData files
            sub_tes(list(str)): paths to TransposonData files
            window_slice(slice): window sizes for calculating overlap
            gene_names(list(str)): gene names to calculate overlap for, default all
            n_workers(int): process count
        """

        self.validate_files(sub_genes)
        self.validate_files(sub_tes)

        self._logger = logging.getLogger(self.__name__.__class__)
        self._proc_mgr = multiprocessing.Manager()
        self.sub_genes = sub_genes
        self.sub_tes = sub_tes
        self.output_dir = output_dir
        self.n_workers = n_workers or max(min(n_workers, multiprocessing.cpu_count), 1)
        self._stop_event = multiprocessing.Event()
        self._workers = None
        self._consumer = None
        self._producer = None
        self._result_queue = None
        self._job_queue = None

        if not os.path.isdir(self.output_dir):
            # TODO change to FileNotFoundError
            raise ValueError("output dir DNE {}".format(self.output_dir))

    def __enter__(self):
        self.start()
        return self

    def start(self):

        self.stop()
        # TODO start each one
        # if not self._workers:

    def __exit__(self, exc_type, exc_val, traceback):

        self.stop()

    def stop(self):

        self._signal_end()
        self.join_workers()
        if self._consumer:
            self._consumer.join()
            self._consumer = None
        if self._producer:
            self._producer.join()
            self._producer = None
        if self._job_queue:
            self._job_queue = None
        if self._result_queue:
            self._result_queue = None

    def join_workers(self):

        if self._workers:
            for w in self._workers:
                w.join()

    def _signal_end(self):
        """Set the stop event and enqueue sentinels."""

        self._stop_event.set()
        if self._job_queue:
            [self._job_queue.put_nowait(Sentinel()) for _ in self._n_workers]

    @staticmethod
    def validate_files(file_list):
        """Raise FileNotFoundError if the file does not exist."""

        for file in file_list:
            if not os.path.isfile(file):
                raise FILE_DNE(file)
