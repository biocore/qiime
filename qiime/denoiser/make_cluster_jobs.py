#!/usr/bin/env python

"""A simple qsub based cluster submission script."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os.path import exists
from os import remove, rename, rmdir, makedirs, close
from subprocess import Popen, PIPE, STDOUT
from tempfile import mkstemp

from skbio.util import create_dir
from burrito.util import ApplicationNotFoundError
from burrito.util import which


# qsub template
# requires format string (walltime, ncpus, nodes, queue, job_name,
# keep_output, command)
QSUB_TEXT = """# Walltime Limit: hh:nn:ss
#PBS -l walltime=%s

# Node Specification:
#PBS -l ncpus=%d -l nodes=%d

# Queue: Defaults to friendlyq
#PBS -q %s

# Mail: options are (a) aborted, (b) begins execution, (e) ends execution
# use -M <email> for additional recipients
# supress email notification
#PBS -m n

# Job Name:
#PBS -N %s

# Keep output
#PBS -k %s

echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------
cd $PBS_O_WORKDIR
%s
"""


def make_jobs(commands, job_prefix, queue, jobs_dir="jobs/",
              walltime="72:00:00", ncpus=1, nodes=1, keep_output="oe"):
    """prepare qsub text files.

    command: list of commands

    job_prefix: a short, descriptive name for the job.

    queue: name of the queue to submit to

    jobs_dir: path to directory where job submision scripts are written

    walltime: the maximal walltime

    ncpus: number of cpus

    nodes: number of nodes

    keep_output: keep standard error, standard out, both, or neither
                 o=std out, e=std err, oe=both, n=neither
    """

    filenames = []
    create_dir(jobs_dir)
    for command in commands:
        fd, job_name = mkstemp(dir=jobs_dir, prefix=job_prefix + "_",
                              suffix=".txt")
        close(fd)
        out_fh = open(job_name, "w")

        out_fh.write(QSUB_TEXT % (walltime, ncpus, nodes, queue, job_prefix,
                                  keep_output, command))
        out_fh.close()
        filenames.append(job_name)
    return filenames


def submit_jobs(filenames, verbose=False):
    """Submit jobs in filenames.

    filenames: list of prepared qsub job scripts, ready to be submitted

    verbose: a binary verbose flag
    """
    if not which("qsub"):
        raise ApplicationNotFoundError("qsub not found. Can't submit jobs.")

    for file in filenames:
        command = 'qsub %s' % file
        result = Popen(command, shell=True, universal_newlines=True,
                       stdout=PIPE, stderr=STDOUT).stdout.read()
        if verbose:
            print result
