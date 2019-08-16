#!/usr/bin/env python3

"""
Consumes gene/transposon data and produces density values.
"""

__author__ = "Michael Teresi, Scott Teresi"


from multiprocessing import Process, Queue, Event
from queue import Empty, Full

from transposon.density import

class DensityInput(object):
    """Contains data to request a density calculation."""

    __slots__ = ['gene_id', 'window']
    def __init__(self, gene_id, window):
        """Initializer.

        Args:
            gene_id (str): the gene to calculate density on.
            window (int): the window size to calculate density on.
        """
        self.gene_id = gene_id
        self.window = window

class DensityOutput(object):
    """Contains data on a density result."""


class DensityWorker(Process):
    """Worker to produce density values provided gene / transposon data."""

    def __init__(self, input_queue, output_queue, genes, transposons, **kwargs):
        """Initializer.

        Args:
            input_queue (multiprocessing.Queue): queue of DensityInput.
            output_queue (multiprocessing.Queue): queue of DensityOutput.
            genes (): gene data.
            transposons (): transposon data.
            kwargs: keyword arguments for multiprocessing.Process.
        """

        super(DensityWorker, self).__init__(**kwargs)
        self.input_queue = input_queue
        self.output_queue = output_queue
        self.genes = genes
        self.transposons = transposons
        self.stop_event = Event()
        # NOTE we should use this: https://gist.github.com/JesseBuesking/10674086
        self._logger = logging.getLogger(__name__)

    def stop(self):
        """End calculations."""

        self.stop_event.set()

    def run(self):
        """Begin calculations; called by Process.start"""

        density_output = None
        while not self.stop_event.is_set():
            if density_output is None:
                density_input = self._pop_request()
                density_output = self._calc_density(density_input)
            else:
                density_output = self._enqueue_result(density_output)

    def _pop_request(self):
        """Consume and return, None if timed out."""

        try:
            # MAGIC NUMBER mostly affects time to exit in our case
            return self.input_queue.get(timeout=0.5)
        except Empty:
            return None

    def _enqueue_result(self, density_output):
        """Produce, return object to enqueue if timed out."""

        try:
            # MAGIC NUMBER mostly affects time to exit in our case
            self.output_queue.put(timeout=0.5)
            return None
        except Full:
            return density_output


    def _calc_density(self, density_input):

        if density_input is None:
            return None

        # TODO left off here
