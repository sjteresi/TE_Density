#!/usr/bin/env python3

"""
Collates density results.

Each density is independent, so the processing order is not important.
Densities are summed for one gene/window pair wrt one TE type.
DensityAccumulator receives the densities and produces
"""

__author__ = "Michael Teresi, Scott Teresi"

from multiprocessing import Process, Queue


class LoggerAccumulator(object):
    """Collates logs from the workers."""

    # SEE https://gist.github.com/JesseBuesking/10674086
    pass

class SubDensityAccumulator(Process):
    """Sums density results for a subset of transposable elements."""

    pass


class DensityAccumulatorMap(dict):
    """Maps a gene name to the accumulated density results."""

    pass


class DensityAccumulator(object):
    """Sums density values provided results."""

    # NOTE maybe use a multiprocessing manager to do the merging?
    # https://stackoverflow.com/questions/8640367/python-manager-dict-in-multiprocessing

    def __init__(self, )
