#!/usr/bin/env python3

"""
Produces density values provided gene/transposon data.
"""

__author__ = "Michael Teresi, Scott Teresi"


from multiprocessing import Process, Queue, Event
from queue import Empty, Full

from transposon.density import rho_left_window, rho_intra, rho_right_window

from transposon.density import

class DensityInput(object):
    """Contains data to request a density calculation."""

    def __init__(self, gene_name, window):
        """Initializer.

        Args:
            gene_name (str): the gene to calculate density on.
            window (int): the window size to calculate density on.
        """

        assert window >= 0
        self.gene_name = str(gene_name)
        self.window = int(window)

class OverlapOutput(object):
    """Contains a density result, prior to division by relevant area."""

    def __init__(self):
        """Initializer.

        Args:

        """
        self.gene_name = str(gene_name)
        self.window = int(window)
        self.te_family = str(te_family)
        self.n_basepairs = int(relevant_basepairs)



class DensityOutput(object):
    """Contains data on a density result."""

    def __init__(self):
        """Initializer.

        Args:
        """
        self.input = density_input
        self.is_family = is_family  # true if family, false if subfamily
        self.family_or_subfamily_name = name
        self.densities = densities
        self.exc_type = None
        self.exc_value = None
        self.exc_traceback = None


class DensityWorker(Process):
    """Calculates TE densities provided gene / transposon data.

    Density values are summed with respect to the TE sub/family.
    The final output is the accumulated densities.
    If a density cannot be completed, the request is put on the error queue.

    Design:
        The worker accumulates densities as it progresses.
        This reduces interprocess communication compared to the case where
            another process(es) would sum each individual density result.
        In said case, the accumulation would need to lock the entry for each
            gene, subfamily combo prior to making the sum.
        The increase in memory is preferred over simplifyin
    """

    def __init__(self, genes, transposons, request_queue, error_queue, **kwargs):
        """Initializer.

        Args:
            genes (data.GeneData): genes.
            transposons (data.TransposableElementsData): transposons.
            input_queue (multiprocessing.Queue): queue of DensityInput.
            output_queue (multiprocessing.Queue): queue of DensityOutput.
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

        self._families = self._extract_families()
        self._sub_families = self._extract_sub_families()

    def stop(self):
        """End calculations."""

        self.stop_event.set()

    def run(self):
        """Begin calculations; called by Process.start"""

        density_output = None
        # maybe break it up like: pop, check, calc, enqueue?
        while not self.stop_event.is_set():
            if density_output is None:
                density_input = self._pop_request()
                density_output = self._calc_density(
                    density_input,
                    self.genes,
                    self.transposons,
                    self.logger)
            else:
                density_output = self._enqueue_result(density_output)

    def _pop_request(self):
        """Consume and return, None if timed out."""

        try:
            # MAGIC NUMBER mostly affects time to exit in our case
            return self.input_queue.get(timeout=0.5)
        except Empty:
            return None

    def _enqueue_result(self, rho_output):
        """Produce, return object to enqueue if timed out."""

        try:
            # MAGIC NUMBER mostly affects time to exit in our case
            self.output_queue.put(timeout=0.5)
            return None
        except Full:
            return density_output

    def _check_gene_name(self, rho_input):

        return rho_input.name in self.genes.names

    def _validate_input(self, rho_input):
        raise NotImplementedError()

    @static_method
    def _extract_families(self, transposons):
        """Return set of the TE families."""
        raise NotImplementedError()

    @static_method
    def _extract_sub_families(self, transposons):
        """Return set of the TE sub-families."""
        raise NotImplementedError()

    def _extract_output(self, left, intra, right, key_list, key):
        """Return an output for the densities wrt TE family | sub-family."""
        raise NotImplementedError()


    @static_method
    def _calc_density(rho_input, genes, transposons logger):
        """Caculate the transposable element densities.

        This is the target function of the worker.

        Args:
            rho_input (DensityOutput): input request.
            logger (logging.Logger): logger to use.
        Returns:
            (DensityOutput): the TE densities.

        """

        # TODO need to try/catch, but might need to define our own exceptions

        if density_input is None:
            return None

        name = rho_input.gene_name
        left = rho_left_window(
            genes.start(name),
            genes.stop(name),
            rho_input.window,
            transposons.starts
            transposons.stops
            transposons.lengths
        )
        intra = rho_intra(
            genes.start(name),
            genes.stop(name),
            rho_input.window,
            transposons.starts,
            transposons.stops,
            transposons.lengths,
        )
        left = rho_right_window(
            genes.start(name),
            genes.stop(name),
            rho_input.window,
            transposons.starts
            transposons.stops
            transposons.lengths
        )
        densities = DensityOutput()
        densities.input = rho_input
        densities.
