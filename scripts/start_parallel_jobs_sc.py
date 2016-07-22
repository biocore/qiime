#!/usr/bin/env python
# File created on 10 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os import makedirs, chmod, getenv, remove
from os.path import exists
from shutil import rmtree
from stat import S_IRWXU
from qiime.util import make_option
from qiime.util import parse_command_line_parameters
from qiime.util import load_qiime_config, qiime_system_call, get_qiime_temp_dir

qiime_config = load_qiime_config()

script_info = {}
script_info[
    'brief_description'] = "Starts parallel jobs on Sun GridEngine queueing systems."
script_info[
    'script_description'] = "Starts multiple jobs in parallel on Sun GridEngine systems. This is designed to work with StarCluster EC2 instances, but may be applicable beyond there."
script_info['script_usage'] = [
    ("Job submission example",
     "Start each command listed in test_jobs.txt in parallel. The run ID for these jobs will be RUNID.",
     "%prog -ms test_jobs.txt RUNID"),
    ("Queue specification example",
     "Submit the commands listed in test_jobs.txt to the specified queue.",
     "%prog -ms test_jobs.txt -q all.q RUNID")
]
script_info['output_description'] = "No output is created."
script_info['required_options'] = []
script_info['optional_options'] = [
    make_option('-m', '--make_jobs', action='store_true',
                help='make the job files [default: %default]'),
    make_option('-s', '--submit_jobs', action='store_true',
                help='submit the job files [default: %default]'),
    make_option('-q', '--queue_name', default=qiime_config['sc_queue'],
                help='the queue to submit jobs to [default: %default]')
]
script_info['version'] = __version__
script_info['disallow_positional_arguments'] = False

# qsub template
QSUB_TEXT = """#!/bin/bash
#$ -V
#$ -q %s
#$ -N %s
#$ -cwd

#######
%s
"""


def write_job_files(output_dir, commands, run_id, queue_name):
    jobs_dir = '%s/jobs/' % output_dir
    job_fps = []
    if not exists(jobs_dir):
        try:
            makedirs(jobs_dir)
        except OSError as e:
            raise OSError("Error creating jobs directory in temp dir: %s" % output_dir +
                          " (specified in qiime_config). Do you have write access?. " +
                          "Original error message follows:\n%s" % str(e))
        paths_to_remove = [jobs_dir]
    else:
        # point paths_to_remove at job_fps
        paths_to_remove = job_fps

    for i, command in enumerate(commands):
        job_fp = '%s/%s%d' % (jobs_dir, run_id, i)
        f = open(job_fp, 'w')
        f.write(QSUB_TEXT %
                (queue_name, run_id + str(i), '\n'.join(command.split(';'))))
        f.close()
        job_fps.append(job_fp)

    return job_fps, paths_to_remove


def run_commands(output_dir, commands, run_id,
                 submit_jobs, keep_temp, queue_name):
    """
    """
    job_fps, paths_to_remove = write_job_files(
        output_dir, commands, run_id, queue_name)

    # Call the jobs
    if submit_jobs:
        for job_fp in job_fps:
            qiime_system_call(' '.join(['qsub', job_fp]))

    # clean up the shell scripts that were created
    if not keep_temp:
        for p in paths_to_remove:
            try:
                # p is file
                remove(p)
            except OSError:
                # p is directory
                rmtree(p)
    return


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.submit_jobs and not opts.make_jobs:
        option_parser.error('Must pass -m if passing -s. (Sorry about this, '
                            'it\'s for backwards-compatibility.)')

    min_args = 2
    if len(args) < min_args:
        option_parser.error('Exactly two arguments are required.')

    output_dir = get_qiime_temp_dir()
    run_commands(output_dir,
                 open(args[0]),
                 args[1],
                 submit_jobs=opts.submit_jobs,
                 keep_temp=True,
                 queue_name=opts.queue_name)

if __name__ == "__main__":
    main()
