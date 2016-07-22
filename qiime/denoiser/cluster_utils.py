#!/usr/bin/env python

"""Some utility functions for operating on a cluster or MP machine."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight", "Nigel Cook", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os import remove, system
from string import join, lowercase
from os.path import exists, join
from time import sleep, time
from random import sample

from asynchat import async_chat
from socket import socket, AF_INET, SOCK_STREAM, gethostname, error
from burrito.util import ApplicationNotFoundError
from burrito.util import which
from qiime.util import load_qiime_config, get_qiime_temp_dir


def submit_jobs(commands, prefix):
    """submit jobs using exe pointed to by cluster_jobs_fp.

    commands: List of commands (strings) that should be executed

    prefix: A uniq prefix used to name submit script
"""
    qiime_config = load_qiime_config()
    CLUSTER_JOBS_SCRIPT = qiime_config['cluster_jobs_fp']

    if not CLUSTER_JOBS_SCRIPT:
        raise ApplicationNotFoundError(
            "cluster_jobs_fp not set in config file!")
    if not (exists(CLUSTER_JOBS_SCRIPT) or which(CLUSTER_JOBS_SCRIPT)):
        raise ApplicationNotFoundError(
            "cluster_jobs_fp not in $PATH or provided as full path!")

    outfilename = join(get_qiime_temp_dir(), "%s_commands.txt" % prefix)
    fh = open(outfilename, "w")
    fh.write("\n".join(commands))
    fh.close()
    cmd = '%s -ms %s %s' % (CLUSTER_JOBS_SCRIPT, outfilename, prefix)
    system(cmd)
    remove(outfilename)


def setup_workers(num_cpus, outdir, server_socket, verbose=True,
                  error_profile=None):
    """Start workers waiting for data.

    num_cpus: number of cores

    outdir: directory were the workers will work in

    server_socket: an open socket to the server

    verbose: verbose flag passed to the workers

    error_profile: filepath to the error profiles, passed to workers

"""
    DENOISE_WORKER = "denoiser_worker.py"
    workers = []
    client_sockets = []
    # somewhat unique id for cluster job
    tmpname = "".join(sample(list(lowercase), 8))

    host, port = server_socket.getsockname()

    # TODO: this should be set to a defined wait time using alarm()
    for i in range(num_cpus):
        name = outdir + ("/%sworker%d" % (tmpname, i))
        workers.append(name)
        cmd = "%s -f %s -s %s -p %s" % (DENOISE_WORKER, name, host, port)

        if verbose:
            cmd += " -v"
        if error_profile:
            cmd += " -e %s" % error_profile

        submit_jobs([cmd], tmpname)
        # wait until the client connects
        # This might be a race condition -> make the client robust
        client_socket, client_address = server_socket.accept()
        client_sockets.append((client_socket, client_address))

    return workers, client_sockets


def adjust_workers(num_flows, num_cpus, worker_sockets, log_fh=None):
    """Stop workers no longer needed.

    num_flows: number of flowgrams

    num_cpus: number of CPUs currently used

    worker_sockets: list of connected sockets

    log_fh: open fh to log file

    Returns new number of CPUs
    """

    qiime_config = load_qiime_config()
    min_per_core = int(qiime_config['denoiser_min_per_core'])
    if(num_flows < (num_cpus - 1) * min_per_core):
        if log_fh:
            log_fh.write("Adjusting number of workers:\n")
            log_fh.write("flows: %d   cpus:%d\n" % (num_flows, num_cpus))
        # TODO: make sure this works with __future__ division
        per_core = max(min_per_core, (num_flows / num_cpus) + 1)
        for i in range(num_cpus):
            if(i * per_core > num_flows):
                worker_sock = worker_sockets.pop()
                worker_sock.close()
                num_cpus = num_cpus - 1
                if log_fh:
                    log_fh.write("released worker %d\n" % i)
        if log_fh:
            log_fh.write("New number of cpus:%d\n" % num_cpus)
    if (num_cpus == 0 or num_cpus != len(worker_sockets)):
        raise ValueError("Adjust_workers screwed up!")
    return num_cpus


def stop_workers(worker_sockets, log_fh=None):
    """Stop all workers.

    worker_sockets:  list of connected sockets

    log_fh: open fh to log file
    """
    for i, worker in enumerate(worker_sockets):
        try:
            worker.send("Server shutting down all clients")
        except error:
            # socket already closed, client dead
            if log_fh:
                log_fh.write(
                    "Worker %s seems to be dead already. Check for runaways!\n" %
                    i)
        worker.close()


def check_workers(workers, worker_sockets, log_fh=None):
    """Check if all workers are still alive. Exit otherwise.

    workers: list of worker names

    worker_sockets: list of connected sockets

    log_fh: open fh to log file

    """

    # Do a dummy send and see if it fails
    for worker, sock in zip(workers, worker_sockets):
        try:
            sock.send("")
        except error:
            if log_fh:
                log_fh.write(
                    "FATAL ERROR\nWorker %s not alive. Aborting\n" %
                    worker)
            stop_workers(worker_sockets, log_fh)
            return False
    return True


def setup_server(port=0, verbose=False):
    """Open a port on the server for workers to connect to.

    port: the port number to use, 0 means let OS decide

    verbose: a verbose flag
    """

    host = gethostname()
    sock = socket(AF_INET, SOCK_STREAM)
    try:
        sock.bind((host, port))
    except error as msg:
        raise error("Could not open Socket on server: " + str(msg))
    sock.listen(5)  # max num of queued connections usually [1..5]
    if verbose:
        print "Server listening on %s" % str(sock.getsockname())
    return sock


def setup_cluster(num_cpus, outdir, verbose, error_profile):
    """Setup server and clients"""

    server_socket = setup_server()
    workers, client_socks_and_adrs = setup_workers(
        num_cpus, outdir, server_socket,
        verbose=verbose,
        error_profile=error_profile)
    # we don't need the client adresses anywhere, so get rid of them
    client_sockets = [sock for sock, addr in client_socks_and_adrs]

    return client_sockets, workers, server_socket


def save_send(socket, data):
    """send data to a socket.

    socket: a connected socket object

    data: string to send over the socket

    """

    # We have no control about how much data the clients accepts,
    # thus we send in chunks until done
    while len(data) > 0:
        try:
            send_data_size = socket.send(data)
            # remove sent portion form data
            data = data[send_data_size:]
        except error as msg:
            # most likely socket busy, buffer full or not yet ready
            sleep(0.01)


def send_flowgram_to_socket(identifier, flowgram, socket, trim=False):
    """send one flowgram over a socket.

    id: identifier of this flowgram

    flowgram: the flowgram itself

    socket: socket to write to

    trim: Boolean flag for quality trimming flowgrams
    """

    if trim:
        flowgram = flowgram.getQualityTrimmedFlowgram()

    # store space separated string representation of flowgram
    # storing this is much quicker than re-generating everyt we send it
    if (not hasattr(flowgram, "spaced_flowgram")):
        spaced_flowgram_seq = " ".join(map(str, flowgram.flowgram))
        flowgram.spaced_flowgram = spaced_flowgram_seq
    else:
        spaced_flowgram_seq = flowgram.spaced_flowgram

    data = "%s %d %s\n" % (identifier, len(flowgram), spaced_flowgram_seq)
    save_send(socket, data)


class ClientHandler(async_chat):

    """A convenience wrapper around a socket to collect incoming data"""

    # This handler is called from the main routine with an open socket.
    # It waits for the client to return its data on the socket and stores it
    # in a global variable result_array. Afterwards the handler is deleted
    # by removing it from the global asyncore map
    # Note: the incomgn socket is expected to be connected upon initialization
    #      and remains connected after this handler is destroyed

    def __init__(self, sock, worker_number, result_array, timing):
        async_chat.__init__(self, sock)
        self.in_buffer = []
        self.set_terminator("--END--")
        self.number = worker_number
        self.results = result_array
        self.timing = timing

    def collect_incoming_data(self, data):
        """Buffer the data"""

        # Note data might come in chunks of arbitrary size
        self.in_buffer.append(data)

    def found_terminator(self):
        """Action performed when the terminator is found."""
        # Note this function is event-triggered

        # Data on sockets comes in chunks of strings. Cat first then split on
        # \n
        data = "".join(self.in_buffer)
        self.results[self.number] = [map(float, (s.split())) for s in data.split("\n")
                                     if s != ""]
        self.in_buffer = []
        self.timing[self.number] = time()
        # delete this channel from the global map, but don't close the socket
        # as we will use it again in the next round.
        # Once global map is empty, asynchronous loop in server will finish
        self.del_channel()
