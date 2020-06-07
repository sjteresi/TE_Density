#!/usr/bin/env python3

"""
Unit test LogListener.
"""

__author__ = "Michael Teresi"


import logging
import logging.handlers
import pytest
import multiprocessing
from multiprocessing import Process, Pool
import threading
import queue

from transposon.log_listener import LogListener

LOG_LEVELS = [logging.DEBUG,
              logging.INFO,
              logging.WARNING,
              logging.ERROR,
              logging.CRITICAL]


def log_in_process(message, level, listener_queue, level_filter=logging.DEBUG):
    """Log a message in a process.

    Args:
        message(str): log message
        level(int): log verbosity level
        listener_queue(multiprocessing.Queue): queue to the listener
        level_filter(int): verbosity level for the process
    """

    args = (message, level, level_filter, listener_queue)
    proc = Process(target=log, args=args)
    proc.start()
    proc.join()

def log(message, level, level_filter, log_queue):
    """Log a message to the log listener queue.

    Args:
        message(str): log message
        level(int): log verbosity level
        level_filter(int): caller verbosity limit
        log_queue(multiprocessing.Queue): log listener Queue
    """

    LogListener.connect(log_queue, level_filter)
    logging.log(level, message)


# NOTE couldn't get this to work with module scope b/c of race conditions?
@pytest.yield_fixture()
def listener():
    """A running log listener."""

    listener_queue = multiprocessing.Queue(-1)
    listener = LogListener(message_queue=listener_queue)
    listener.start()
    yield listener
    listener.stop()
    listener.join()


@pytest.yield_fixture()
def feedback():
    """Messages received by the root logger."""

    # setup root logger to put records in a queue for feedback
    logger = logging.getLogger()
    feedback_queue = multiprocessing.Queue(-1)
    q_handler = logging.handlers.QueueHandler(feedback_queue)
    logger.addHandler(q_handler)
    yield feedback_queue
    logger.removeHandler(q_handler)

def test_instance_default():
    """Can we instantiate the listener?"""

    listener = LogListener()


def test_instance_queue():
    """Can we instantiate the listener w/ and input queue?"""

    log_queue = multiprocessing.Queue()
    listener = LogListener(log_queue)


def test_instance_fail():
    """Does the initializer raise if given the wrong queue?"""

    with pytest.raises(TypeError) as excinfo:
        listener = LogListener(queue.Queue())


@pytest.mark.parametrize("level", LOG_LEVELS)
def test_message(listener, feedback, level):

    message = "I'll be back " + str(level)
    log_in_process(message, level, listener.log_queue, level_filter=logging.DEBUG)
    answer = feedback.get(timeout=0.1)  # MAGIC experimental
    assert answer.levelno == level
    assert answer.msg == message


if __name__ == "__main__":
    pytest.main(['-s', __file__])  # for convenience
