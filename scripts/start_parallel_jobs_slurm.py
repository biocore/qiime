#!/usr/bin/env python

"""A simple slurm based cluster submission script."""

__author__ = "Simon Jacobs"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder", "Rob Knight", "Greg Caporaso", "Jai Ram Rideout", "Evan Bolyen", "Simon Jacobs"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Simon Jacobs"
__email__ = "sdjacobs@uchicago.edu"

from optparse import OptionParser
from os.path import exists, normpath, sep
from os import remove, rename, rmdir, makedirs, close
from tempfile import mkstemp

from qiime.util import  make_option,\
    parse_command_line_parameters, load_qiime_config, qiime_system_call

from qiime.denoiser.make_cluster_jobs import make_jobs, submit_jobs

qiime_config = load_qiime_config()

script_info = {}
script_info[
    'brief_description'] = "Starts multiple jobs in parallel on slurm based multiprocessor systems."
script_info[
    'script_description'] = "This script is designed to start multiple jobs in parallel on cluster systems with a slurm based scheduling system."
script_info['script_usage'] = [
    ("Job submission example",
     "Start each command listed in test_jobs.txt in parallel. The run ID for these jobs will be RUNID.",
     "%prog -ms test_jobs.txt RUNID"),
    ("Queue specification example",
     "Submit the commands listed in test_jobs.txt to the specified queue.",
     "%prog -ms test_jobs.txt -q himem RUNID"),
    ("Jobs output directory specification example",
     "Submit the commands listed in test_jobs.txt, with the jobs put under the "
     "my_jobs/ directory.",
     "%prog -ms test_jobs.txt -j my_jobs/ RUNID")
]

default_slurm_queue_desc = qiime_config['slurm_queue'] or "slurm's default"
default_slurm_memory_desc = qiime_config['slurm_memory'] or "slurm's default"
default_slurm_time_desc = qiime_config['slurm_time'] or "slurm's default"

script_info['output_description'] = "No output is created."
script_info['required_options'] = []
script_info['optional_options'] = [
    make_option('-m', '--make_jobs', action='store_true',
                help='make the job files [default: %default]',
                default=False),

    make_option('-s', '--submit_jobs', action='store_true',
                help='submit the job files [default: %default]',
                default=False),

    make_option('-q', '--queue',
                help=('name of queue to submit to '
                      '[default: %s]' % default_slurm_queue_desc),
                default=qiime_config['slurm_queue']),

    make_option('-K', '--mem_per_cpu',
                help=('megabytes of memory to request per '
                      'CPU [default: %s]' % default_slurm_memory_desc),
                default=qiime_config['slurm_memory']),

    make_option('-j', '--job_dir',
                help='directory to store the jobs [default: %default]',
                default="jobs/"),

    make_option('-t', '--time',
                help=('run time limit of the jobs in dd-hh:mm:ss format '
                      '[default: %s]' % default_slurm_time_desc),
                default=qiime_config['slurm_time']),
]

script_info['version'] = __version__
script_info['disallow_positional_arguments'] = False

def make_jobs(commands, job_prefix, queue, jobs_dir="jobs/"):
    filenames = []
    for command in commands:
        fd, job_name = mkstemp(dir=jobs_dir, prefix=job_prefix + "_",
                              suffix=".txt")
        close(fd)
        out_fh = open(job_name, "w")
        out_fh.write("#!/bin/sh\n")
        out_fh.write(command)
        out_fh.close()
        filenames.append(job_name)
    return filenames


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    if opts.submit_jobs and not opts.make_jobs:
        option_parser.error('Must pass -m if passing -s. (Sorry about this, '
                            'it\'s for backwards-compatibility.)')

    min_args = 2
    if len(args) != min_args:
        option_parser.error('Program requires <commands file> and '
                            '<job prefix>')

    if (len(args[1]) > 10 or len(args[1]) == 0):
        option_parser.error('job prefix must be 1-10 characters long')

    if(not exists(opts.job_dir)):
        try:
            makedirs(opts.job_dir)
        except OSError:
            exit(" Jobs directory can not be created. "
                 "Check for permissions or file with the same name: %s\n"
                 % opts.job_dir)

    commands = list(open(args[0]))
    job_prefix = args[1]

    if opts.mem_per_cpu is not None:
        mem_per_cpu = " --mem-per-cpu=" + opts.mem_per_cpu
    else:
        mem_per_cpu = ""

    if opts.queue is not None:
        queue = " -p " + opts.queue
    else:
        queue = ""

    if (opts.make_jobs):
        filenames = make_jobs(
            commands,
            job_prefix,
            opts.queue,
            opts.job_dir)
    else:
        exit("Should we ever get here???")

    if opts.time is not None:
        time = " --time=" + opts.time
    else:
        time = ""

    if (opts.submit_jobs):
        for f in filenames:
            cmd = "".join([
                    "sbatch",
                    queue,
                    " -J ", job_prefix,
                    mem_per_cpu,
                    time,
                    " -o ", normpath(opts.job_dir), sep, job_prefix, "_%j.out",
                    " ", f
                ])
            qiime_system_call(cmd, shell=True)

if __name__ == "__main__":
    main()
