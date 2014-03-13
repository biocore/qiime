#!/usr/bin/env python
"""Tests for functions in utils."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

import signal
import os
from os import remove, rmdir, environ, mkdir
from time import sleep, time
from os.path import exists
from StringIO import StringIO
from socket import error

from unittest import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.flowgram import Flowgram
from cogent import Sequence
from cogent.app.util import ApplicationNotFoundError
from cogent.util.misc import remove_files
from qiime.util import get_tmp_filename

from qiime.util import load_qiime_config
from qiime.denoiser.cluster_utils import submit_jobs, setup_server,\
    setup_workers, adjust_workers, stop_workers, check_workers

# timeout handling taken from test_workflow.py


class TimeExceededError(Exception):
    pass

allowed_seconds_per_test = 60


def timeout(signum, frame):
    raise TimeExceededError("Test failed to run in allowed time (%d seconds)."
                            % allowed_seconds_per_test)


class TestUtils(TestCase):

    def setUp(self):

        self.home = environ['HOME']
        self.server_socket = None
        self.tmp_result_file = get_tmp_filename(tmp_dir=self.home,
                                                prefix="/test_hello_",
                                                suffix=".txt")
        self.tmp_dir = get_tmp_filename(tmp_dir=self.home,
                                        prefix="test_cluster_util",
                                        suffix="/")

        self.files_to_remove = [self.tmp_result_file]
        self.command = "echo hello > %s\n" % self.tmp_result_file

        signal.signal(signal.SIGALRM, timeout)
        # set the 'alarm' to go off in allowed_seconds seconds
        signal.alarm(allowed_seconds_per_test)

    def tearDown(self):
        """Clean up tmp files."""

        # turn off the alarm
        signal.alarm(0)

        remove_files(self.files_to_remove, False)
        if self.server_socket:
            self.server_socket.close()
        # give clients time to clean up
        sleep(1)
        if exists(self.tmp_dir):
            try:
                rmdir(self.tmp_dir)
            except OSError:
                # give clients some more time, fail if still error
                sleep(5)
                rmdir(self.tmp_dir)

    def _setup_server_and_clients(self):

        self.server_socket = setup_server()
        mkdir(self.tmp_dir)
        workers, client_sockets = setup_workers(4, self.tmp_dir,
                                                self.server_socket,
                                                verbose=False)

        return workers, [sock for sock, addr in client_sockets]

    def test_submit_jobs(self):
        """submit_jobs executes commands on the cluster"""

        submit_jobs([self.command], prefix="test_job")
        # Try and wait ten times, could be made nicer with alarm()
        for i in range(10):
            if exists(self.tmp_result_file):
                observed_text = "".join(list(open(self.tmp_result_file)))
                self.assertEqual(observed_text, "hello\n")
                return
            else:
                sleep(10)
        # if we get here we failed
        self.fail("The test job apparently never finished.\n"
                  + "check the jobs error log and check the queue status\n.")

    def test_setup_workers(self):
        """setup_workers starts clients"""

        mkdir(self.tmp_dir)
        self.server_socket = setup_server()
        workers, client_sockets = setup_workers(4, self.tmp_dir,
                                                self.server_socket,
                                                verbose=False)
        self.assertEqual(len(workers), 4)
        # Try sending some data, return of send should be length of message
        for client_sock, addr in client_sockets:
            self.assertEqual(client_sock.send("Hello"), 5)
        # workers die, once server_socket closes in tearDown()

    def test_adjust_workers(self):
        """adjust_workers stops clients"""

        workers, client_sockets = self._setup_server_and_clients()
        last_sock = client_sockets[-1]

        qiime_config = load_qiime_config()
        min_per_core = int(qiime_config['denoiser_min_per_core'])

        # no sockets get stopped
        self.assertEqual(
            adjust_workers(
                4 *
                min_per_core -
                1,
                4,
                client_sockets),
            4)
        # if we can send something the socket is still alive
        self.assertEqual(last_sock.send("Hello"), 5)

        # now, kill one client
        self.assertEqual(
            adjust_workers(
                3 *
                min_per_core -
                1,
                4,
                client_sockets),
            3)
        # socket should be closed
        self.assertRaises(error, last_sock.send, "Hello")

    def test_stop_workers(self):
        """stop_workers terminates all clients"""

        workers, client_sockets = self._setup_server_and_clients()

        stop_workers(client_sockets)
        for client_socket in client_sockets:
            self.assertRaises(error, client_socket.send, "hello")

    def test_stop_workers_on_closed_socket(self):
        """stop_workers terminates all clients"""

        # Repeat test but this time close one of the sockets early.
        # simulates crashed client
        workers, client_sockets = self._setup_server_and_clients()
        client_sockets[-1].close()
        fake_fh = StringIO()
        stop_workers(client_sockets, fake_fh)
        self.assertEqual(fake_fh.getvalue(),
                         "Worker 3 seems to be dead already. Check for runaways!\n")
        for client_socket in client_sockets:
            self.assertRaises(error, client_socket.send, "hello")

    def test_check_workers(self):
        """check_workers checks for dead workers"""

        workers, client_sockets = self._setup_server_and_clients()

        self.assertTrue(check_workers(workers, client_sockets))

        # Now close and terminate a client, wait and check again
        client_sockets[0].close()
        self.server_socket.close()
        sleep(1)
        self.assertFalse(check_workers(workers, client_sockets))

    def test_setup_server(self):
        """setup_server opens a port and listens"""

        self.server_socket = setup_server()
        host, port = self.server_socket.getsockname()
        # Not much to test here, if we get something back it should work
        self.assertGreaterThan(port, 1023)
        # binding to a known port should fail
        self.assertRaises(error, setup_server, 80)

    def test_save_send(self):
        """save_send reliably send data to a socket"""
        # Don't really know how to test this effectively...
        # Would require to simulate a blocking socket on the recipient side...
        pass

if __name__ == "__main__":
    main()
