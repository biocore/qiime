#!/usr/bin/env python

from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Daniel McDonald", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

from qiime.util import make_option
from qiime.util import (parse_command_line_parameters, load_qiime_config,
                        get_options_lookup)
from os import popen, system, makedirs, mkdir
from os.path import split, splitext, join
from subprocess import check_call, CalledProcessError
from cogent.core.tree import TreeNode
from random import choice
from time import sleep, time
from qiime.parallel.merge_otus import start_job, local_job, torque_job, \
    job_complete, initial_has_dependencies, initial_nodes_to_merge, \
    mergeorder, mergetree
qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Parallel merge BIOM tables"""
script_info[
    'script_description'] = """This script works like the merge_otu_tables.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example""",
     """Merge the OTU tables $PWD/t1.biom,$PWD/t2.biom,$PWD/t3.biom,$PWD/t4.biom and write the resulting output table to the $PWD/merged/ directory.""",
     """%prog -i $PWD/t1.biom,$PWD/t2.biom,$PWD/t3.biom,$PWD/t4.biom -o $PWD/merged/"""))
script_info[
    'output_description'] = """The output consists of many files (i.e. merged_table.biom, merged_table.log and all intermediate merge tables). The .biom file contains the result of merging the individual BIOM tables. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file, the dependency hierarchy and runtime information for each individual merge."""

script_info['required_options'] = [
    make_option('-i', '--input_fps', type='existing_filepaths',
                help='the otu tables in biom format (comma-separated)'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='the output otu table directory path')]

script_info['optional_options'] = [
    make_option('-C', '--cluster', action='store_true', default=False,
                help="Submit to a torque cluster"),
    options_lookup['seconds_to_sleep'],
    options_lookup['job_prefix']]
script_info['version'] = __version__

RANDOM_JOB_PREFIX_CHARS = "abcdefghigklmnopqrstuvwxyz"
RANDOM_JOB_PREFIX_CHARS += RANDOM_JOB_PREFIX_CHARS.upper()
RANDOM_JOB_PREFIX_CHARS += "0123456790"


def get_random_job_prefix(fixed_prefix='',
                          max_job_prefix_len=10,
                          leading_trailing_underscores=True):
    """ Return a string to use as job prefix """

    length = max_job_prefix_len - len(fixed_prefix)
    if leading_trailing_underscores:
        length -= 2

    result = [choice(RANDOM_JOB_PREFIX_CHARS) for i in range(length)]
    if leading_trailing_underscores:
        return fixed_prefix + '_' + ''.join(result) + '_'
    else:
        return fixed_prefix + ''.join(result)


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    input_fps = opts.input_fps
    output_dir = opts.output_dir
    seconds_to_sleep = opts.seconds_to_sleep
    verbose = opts.verbose

    merge_otus_serial_script = 'merge_otu_tables.py'
    created_temp_paths = []

    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with MOTU
    job_prefix = opts.job_prefix or get_random_job_prefix('MOTU')

    # A temporary output directory is created in output_dir named
    # job_prefix. Output files are then moved from the temporary
    # directory to the output directory when they are complete, allowing
    # a poller to detect when runs complete by the presence of their
    # output files.
    working_dir = '%s/%s' % (output_dir, job_prefix)
    try:
        makedirs(working_dir)
    except OSError:
    # working dir already exists
        pass

    import os.path
    # wrapper log output contains run details
    log_fp = os.path.join(output_dir, 'parallel_merge_otus.log')

    wrapper_log_output = open(log_fp, 'w')
    wrapper_log_output.write("Parallel merge output\n\n")

    # construct the dependency tree
    import os

    for f in input_fps:
        if not os.path.exists(f):
            raise IOError("%f does not exist!" % f)

    tree = mergeorder(input_fps, working_dir)

    if verbose:
        print tree.asciiArt()

    wrapper_log_output.write('Dependency tree:\n')
    wrapper_log_output.write(tree.asciiArt())
    wrapper_log_output.write('\n\n')
    wrapper_log_output.flush()

    to_process = initial_nodes_to_merge(tree)
    has_dependencies = initial_has_dependencies(tree, to_process)

    # loop until the whole shabang is done
    pending = []  # jobs that are currently running

    while not tree.Processed:
        # check if we have nodes to process, if so, shoot them off
        for node in to_process:
            if opts.cluster:
                start_job(node, merge_otus_serial_script,
                          qiime_config['torque_queue'], wrap_call=torque_job)
            else:
                start_job(node, merge_otus_serial_script,
                          qiime_config['torque_queue'], wrap_call=local_job)

            wrapper_log_output.write(node.FullCommand)
            wrapper_log_output.write('\n')
            wrapper_log_output.flush()

            pending.append(node)
        to_process = set([])

        # check running jobs
        current_pending = []
        for pending_node in pending:
            # if we're complete, update state
            if job_complete(pending_node):
                wrapper_log_output.write("Node %s completed in %f seconds" %
                                         (pending_node.Name, pending_node.TotalTime))
                wrapper_log_output.write('\n')
                wrapper_log_output.flush()
            else:
                current_pending.append(pending_node)
        pending = current_pending

        # check for new jobs to add
        current_dependencies = []
        for dep_node in has_dependencies:
            # if children are satisfied, then allow for processing
            # the logic here is odd to handle the case where an internal node
            # has both a tip that is a child and child that is an internal node
            children_are_complete = [(c.Processed or c.istip())
                                     for c in dep_node.Children]
            if all(children_are_complete):
                to_process.add(dep_node)
            else:
                current_dependencies.append(dep_node)
        has_dependencies = current_dependencies

        sleep(seconds_to_sleep)
    os.rename(tree.FilePath, "%s/%s" % (output_dir, "merged.biom"))
if __name__ == '__main__':
    main()
