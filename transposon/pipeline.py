#!/usr/bin/env python3

"""
Setup a sequence of objects to calculate density using multiprocessing.
"""


from transposon.producer import 
from transposon.worker import DensityInput, DensityOutput, DensityWorker
from transposon.accumulator import DensityAccumulator

class DensityPipeline(object):
    """Contains the producer / consumers for calculating TE density."""

    def __init__(self):
