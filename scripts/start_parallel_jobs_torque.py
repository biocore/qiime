#!/usr/bin/env python

"""A simple qsub based cluster submission script."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder", "Rob Knight", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from optparse import OptionParser
from os.path import exists
from os import remove, rename, rmdir, makedirs

from qiime.util import  make_option,\
    parse_command_line_parameters, load_qiime_config

from qiime.denoiser.make_cluster_jobs import make_jobs, submit_jobs

qiime_config = load_qiime_config()

script_info = {}
script_info[
    'brief_description'] = "Starts multiple jobs in parallel on torque/qsub based multiprocessor systems."
script_info[
    'script_description'] = "This script is designed to start multiple jobs in parallel on cluster systems with a torque/qsub based scheduling system."
script_info['script_usage'] = [
    ("Job submission example",
     "Start each command listed in test_jobs.txt in parallel. The run ID for these jobs will be RUNID.",
     "%prog -ms test_jobs.txt RUNID"),
    ("Queue specification example",
     "Submit the commands listed in test_jobs.txt to the specified queue.",
     "%prog -ms test_jobs.txt -q friendlyq RUNID"),
    ("Jobs output directory specification example",
     "Submit the commands listed in test_jobs.txt, with the jobs put under the "
     "my_jobs/ directory.",
     "%prog -ms test_jobs.txt -j my_jobs/ RUNID")
]
script_info['output_description'] = "No output is created."
script_info['required_options'] = []
script_info['optional_options'] = [
    make_option('-m', '--make_jobs', action='store_true',
                help='make the job files [default: %default]'),

    make_option('-s', '--submit_jobs', action='store_true',
                help='submit the job files [default: %default]'),

    make_option('-q', '--queue',
                help='name of queue to submit to [default: %default]',
                default=qiime_config['torque_queue']),

    make_option('-j', '--job_dir',
                help='directory to store the jobs [default: %default]',
                default="jobs/"),

    make_option('-w', '--max_walltime', type='int',
                help='maximum time in hours the job will run for [default: %default]',
                default=72),

    make_option('-c', '--cpus', type='int',
                help='number of CPUs to use [default:%default]',
                default=1),

    make_option('-n', '--nodes', type='int',
                help='number of nodes to use [default:%default]',
                default=1)

]

script_info['version'] = __version__
script_info['disallow_positional_arguments'] = False


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

    commands = list(open(args[0]))
    job_prefix = args[1]

    if(not exists(opts.job_dir)):
        try:
            makedirs(opts.job_dir)
        except OSError:
            exit(" Jobs directory can not be created. "
                 "Check for permissions or file with the same name: %s\n"
                 % opts.job_dir)

    if (opts.make_jobs):
        filenames = make_jobs(
            commands,
            job_prefix,
            opts.queue,
            opts.job_dir,
            (str(opts.max_walltime) + ":00:00"),
            opts.cpus,
            opts.nodes)
    else:
        exit("Should we ever get here???")
    if (opts.submit_jobs):
        submit_jobs(filenames, opts.verbose)

if __name__ == "__main__":
    main()
