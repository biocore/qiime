#!/usr/bin/env python

"""A simple client waiting for data to clean up 454 sequencing data"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight", "Nigel Cook"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os import remove, makedirs
from os.path import exists, split, dirname
from time import sleep, time
from subprocess import Popen, PIPE, STDOUT
from asyncore import dispatcher, loop
from asynchat import async_chat
from socket import socket, AF_INET, SOCK_STREAM, gethostname, error

from qiime.util import get_qiime_project_dir
from qiime.denoiser.utils import get_denoiser_data_dir,\
    init_flowgram_file, get_flowgram_ali_exe


def setup_worker(fp, server_addr, port, counter=0, verbose=False,
                 error_profile=None):
    """ This routine starts the worker.

    fp: fp that should be used to store intermediate data

    server_addr: ip address of server

    port: port on server to connect to

    counter: counts each round of denoising

    verbose: verbose flag

    error_profile: path to error profile .dat file
    """
    if fp is None:
        raise ValueError("setup_worker needs file path for worker")
    log_fh = None
    if verbose:
        dir = dirname(fp + ".log")
        if not exists(dir):
            makedirs(dir)
        log_fh = open(fp + ".log", "a", 0)

    # use local tmp if possible
    new_fp = fp
    if exists("/tmp"):
        new_fp = "/tmp/" + split(fp)[1]

    # set up the workers and start the loop
    worker = DenoiseWorker(new_fp, server_addr, port, counter=counter,
                           log_fh=log_fh, error_profile=error_profile)
    # this asyncore loop will run until the server closes the connection
    loop()
    # we're done


def _process_data(this_round_fp, log_fh=None, error_profile=None):
    """compute alignment scores for flowgrams in this_round_fp.

    this_round_fp: fp to input data
    log_fh: fh to log file
    error_profile: path to error profile
    """

    # we have data!
    cmd = "%s -relscore_pairid %s %s.dat"\
        % (get_flowgram_ali_exe(),
           error_profile, this_round_fp)
    proc = Popen(cmd, shell=True, universal_newlines=True,
                 stdout=PIPE, stderr=PIPE)
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        host = gethostname()
        if log_fh:
            log_fh.write(
                "An error occured on %s at %f\n%s" %
                (host, time(), cmd))
            log_fh.write(stderr)
            log_fh.close()
        raise RuntimeError("Worker process crashed. Aborting...!\n" +
                           "Note: You need to kill the other jobs yourself")
    if log_fh:
        log_fh.write(this_round_fp + "... done!\n")
    return stdout


class DenoiseWorker(async_chat):

    def __init__(self, fp, address, port=0, log_fh=None,
                 counter=0, error_profile=None):
        """connect to server and init object.

        fp: file path base name to use for intermediate files
        address: server address
        port: server port to connect to
        log_fh: if not None log some messages to log_fh
        counter: counts each round of denoising
        """

        async_chat.__init__(self)
        self.log_fh = log_fh
        self.counter = counter
        self.fp = fp
        self.error_profile = error_profile
        self.in_buffer = []

        self.create_socket(AF_INET, SOCK_STREAM)
        if self.log_fh:
            self.log_fh.write(
                "Client %s trying to connect to %s:%s\n" %
                (fp, address, port))
        # loop until connection to server established
        while True:
            try:
                self.connect((address, port))
                if self.log_fh:
                    self.log_fh.write(
                        "client connected to %s\n" %
                        str(self.getpeername()))
                break
            except error as msg:
                if self.log_fh:
                    self.log_fh.write("Server not ready: %s\n" % msg)
                sleep(1)

        self.set_terminator("--END--")

    def handle_accept(self):
        """called when the other side connects."""
        pass

    def handle_connect(self):
        if self.log_fh:
            self.log_fh.write("Server connection established\n")

    def handle_close(self):
        """Called when the other side closes connection"""
        if self.log_fh:
            self.log_fh.write(
                "Server closed connection at %s. Shutting down.\n" %
                time())
        self.close()

    def collect_incoming_data(self, data):
        """Buffer the data"""
        self.in_buffer.append(data)

    def found_terminator(self):
        """Action performed when the terminator is found."""

        # Note this function is event-triggered
        # This means we received all necessary data
        if self.log_fh is not None:
            self.log_fh.write("Data for round %d received: %s\n"
                              % (self.counter, time()))

        this_round_fp = "%s_%d" % (self.fp, self.counter)

        (fh, filename) = init_flowgram_file(
            filename=this_round_fp + ".dat", n=0)
        for chunk in self.in_buffer:
            fh.write(chunk)
        fh.close()
        self.in_buffer = []

        result = _process_data(this_round_fp, self.log_fh, self.error_profile)

        remove(this_round_fp + ".dat")

        # return results to server
        # Do we need buffering here?
        # No, push() does the buffering for us
        self.push(result)
        self.push("--END--")
        self.counter += 1
