#!/usr/bin/env python3

"""
Generic consumer / producer with Process.

Inherit and implement `Worker.execute_job`.
Implement your 'job' container, enqueue to the job queue.
Consumes jobs from the job queue, passes them to `execute_job`.
Stop processing using the Sentinel and / or stop event.
"""

from abc import ABC, abstractmethod
from collections import namedtuple
import logging
from multiprocessing import Process
import queue

Sentinel = namedtuple("Sentinel", [])


class WorkerProcess(Process, ABC):
    """Consumes jobs, produces results."""

    TIMEOUT = 0.2  # MAGIC arbitrary, affects time to shutdown if not using sentinel

    def __init__(self, job_queue, result_queue, stop_event):
        """Initializer.

        Args:
            job_queue(Queue): input queue
            result_queue(Queue): output queue
            stop_event(Event): end if set
        """

        super().__init__()
        self._logger = logging.getLogger(self.__class__.__name__)
        self.input = job_queue
        self.output = result_queue
        self.stop_event = stop_event

    @abstractmethod
    def execute_job(self, job):
        """Target function for the worker.

        Args:
            job (tuple): caller's container pulled from the job queue
        Returns:
            tuple: caller's container to be enqueued to the result queue
        """

        pass

    def run(self):
        """Get jobs, put results.

        Returns if the stop event is set or the sentinel is received.
        """

        job = None
        result = None
        while not self.stop_event.is_set():

            if not self._send_result(result):
                continue
            else:
                result = None

            try:
                job = job or self.input.get(timeout=self.TIMEOUT)
            except queue.Empty:
                continue
            else:
                if isinstance(job, Sentinel):
                    self.input.put_nowait(job)
                    break

            result = self.execute_job(job)
            job = None

    def _send_result(self, result):
        """True if the result was enqueued, True if 'result' is None."""

        if result is None:
            return True

        success = False
        try:
            self.output.put(result, timeout=self.TIMEOUT)
        except queue.Full:
            self._logger.warning("output queue is full!")
        else:
            sucess = True

        return True

