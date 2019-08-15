#!/usr/bin/env python3

"""
Consumes gene/transposon data and produces density values.
"""

__author__ = "Michael Teresi, Scott Teresi"


from multiprocessing import Process, Queue


class DensityInput(object):
    """Contains config data to request a density calculation."""

    def __init__(self,):

class DensityOutput(object):
    """Contains data on a density result."""


class DensityWorker(object):
    """Worker to produce density values provided gene / transposon data."""

    def __init__(self, genes, transposons, input_queue, output_queue, logger):
        """Initialize."""

    def __enter__(self):

        self.start()

    def __exit__(self, exc_type, exc_value, exc_traceback):

        self.stop()

    def start(self):
        """Begin calculating densities."""


    def stop(self):
        """End calculations."""
