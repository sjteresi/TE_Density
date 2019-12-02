#!/usr/bin/env python3

"""
Produces requests for density calculations.
"""


from transposon.worker import DensityInput


from multiprocessing import Process, Queue, Event
from queue import Empty, Full

from transposon.worker import DensityInput


def input_sweep(window_list, gene_names):
    for window in window_list:
        for gene in gene_names:
            rho_input = DensityInput(gene, window)
            yield rho_input

class DensityRequestor(object):
    """Produces DensityInput objects and tracks failures."""

    def __init_(self, windows, gene_

