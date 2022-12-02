#!/usr/bin/env python3

"""
Unit test WorkerProcess.
"""

__author__ = "Michael Teresi"


from collections import namedtuple
import coloredlogs
import logging
import os
import pytest
import multiprocessing

from transposon.worker import WorkerProcess


class FakeWorker(WorkerProcess):
    """Implementation of WorkerProcess for testing."""

    def execute_job(self, job):
        """Arbitrary job to square an input."""

        logging.debug("got {}".format(job))
        return (job * job, os.getpid())


@pytest.fixture()
def worker():
    mgr = multiprocessing.Manager
    input = multiprocessing.Queue()
    output = multiprocessing.Queue()
    stop = multiprocessing.Event()
    yield FakeWorker(input, output, stop)

@pytest.fixture()
def worker_running(worker):
    """Yield a running process."""
    worker.start()
    yield worker
    worker.stop_event.set()
    worker.join()

def test_start_stop(worker):
    """Can we start / stop the worker?"""

    worker.start()
    worker.stop_event.set()
    worker.join()

def test_answer(worker_running):
    """Do we get feedback?"""
    number = 4
    worker_running.input.put(4)
    answer, worker_pid = worker_running.output.get(timeout=1)
    assert answer == number**2

def test_newprocess(worker_running):
    """Is the worker in a different process?"""

    worker_running.input.put(4)
    answer, worker_pid = worker_running.output.get(timeout=1)
    assert os.getpid() != worker_pid


if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=logging.INFO)
    pytest.main(['-s', __file__])  # for convenience
