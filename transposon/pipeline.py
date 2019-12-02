#!/usr/bin/env python3

"""
Setup a sequence of objects to calculate density using multiprocessing.
"""


from multiprocessing import Process, Queue, Event

from transposon.producer import 
from transposon.worker import DensityInput, DensityOutput, DensityWorker
from transposon.accumulator import DensityAccumulator


class PipelineStage(object):
    """Base object containing a process that runs a worker."""

    def __init__(self):

        self._stop_event = Event()
        self.input_queue = Queue.Queue()
        self.output_queue = Queue.Queue()

    def __enter__(self):
        self.start()

    def __exit__(self, exc_type, exc_val, exc_trace):
        self._stop_event.set()

        


class DensityPipeline(object):
    """Contains the producer / consumers for calculating TE density."""

    def __init__(self):
