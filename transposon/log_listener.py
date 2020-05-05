#!/usr/bin/env python3

"""
Threaded logger for multiprocessing logging.

Receives log records through a queue using the QueueHandler.
Handles the log on receipt.
Provides helper functions to setup the client.
Use logger calls as normal after 'connecting' the process to the QueueHandler.
"""

__author__ = "Michael Teresi"


from functools import partial
import logging
import logging.handlers
import multiprocessing
import threading
from threading import Event, Thread
import queue


class LogListener(Thread):
    """Threaded server to handle logs from sub processes.

    Pass the listener's message queue to the worker processes.
    Setup the worker's logger using self.connect helper function.
    """

    def __init__(self, message_queue=None):
        """Initializer.

        Args:
            message_queue(multiprocessing.Queue): communication pathway for log records
        """

        super().__init__()
        self.log_queue = message_queue or multiprocessing.Queue(-1)
        self.stop_event = Event()

        # if not isinstance(self.log_queue, multiprocessing.queues.Queue):
        #     raise TypeError("expecting {} but got {}"
        #                     .format(multiprocessing.queues.Queue,
        #                             type(self.log_queue)))

    def run(self):
        """Process logs."""

        while not self.stop_event.is_set():
            try:
                record = self.log_queue.get(timeout=1)
                if record is None:
                    break  # poison pill
            except queue.Empty:
                continue
            except (KeyboardInterrupt, SystemExit):
                raise
            else:
                logger = logging.getLogger(record.name)
                logger.handle(record)

    def stop(self):
        """End target function."""

        self.log_queue.put_nowait(None)
        self.stop_event.set()

    def connect_callback(self, level):
        """Callback to configure a client.

        Intended for use in the caller's process.

        Args:
            level(int): log level filter of the process
        """

        return partial(self.connect, self.log_queue, level)

    @staticmethod
    def connect(log_queue, level):
        """Configure local logger to send to the provided queue.

        Intended for use in the caller's process.

        Args:
            log_queue(Queue): log message queue from the listener
            level(logging.Level): filter level, e.g. DEBUG
        """

        q_handler = logging.handlers.QueueHandler(log_queue)
        root_logger = logging.getLogger()
        root_logger.addHandler(q_handler)
        root_logger.setLevel(level)
