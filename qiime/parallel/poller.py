#!/usr/bin/env python

from __future__ import division
from time import sleep
from optparse import OptionParser
from os import getenv, remove
from os.path import exists, isdir
from shutil import rmtree
from skbio.util import remove_files
from qiime.parse import parse_tmp_to_final_filepath_map_file

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


def get_function_handle(s):
    last_dot = s.rindex('.')
    module_name = s[:last_dot]
    function_name = s[last_dot + 1:]
    module = __import__(module_name, globals(), locals(), [function_name])
    function = eval('module.%s' % function_name)
    return function


def remove_all(paths_to_remove):
    for path in paths_to_remove:
        if isdir(path):
            rmtree(path)
        else:
            try:
                remove(path)
            except OSError:
                # File doesn't exist
                pass
    return


def basic_check_run_complete_f(f):
    """ Return True if all filepaths exist

        f: file containing list of filepaths

        example f:
         f1.txt
         f2.txt
         f3.txt

        If f contains the three lines above, this function would return
         False if any of these three files did not exist, and True otherwise.

    """
    filepaths = [l.strip() for l in f]
    for fp in filepaths:
        if not exists(fp):
            return False
    return True


def basic_process_run_results_f(f):
    """ Copy each list of infiles to each outfile and delete infiles

        f: file containing one set of mapping instructions per line

        example f:
         f1.txt f2.txt f3.txt f_combined.txt
         f1.log f2.log f3.log f_combined.log

        If f contained the two lines above, this function would
         concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
         and f1.log, f2.log, and f3.log into f_combined.log
    """
    infiles_lists, out_filepaths = parse_tmp_to_final_filepath_map_file(f)
    for infiles_list, out_filepath in zip(infiles_lists, out_filepaths):
        try:
            of = open(out_filepath, 'w')
        except IOError:
            raise IOError("Poller can't open final output file: %s" % out_filepath +
                          "\nLeaving individual jobs output.\n Do you have write access?")

        for fp in infiles_list:
            for line in open(fp):
                of.write('%s\n' % line.strip('\n'))
        of.close()
    # It is a good idea to have your clean_up_callback return True.
    # That way, if you get mixed up and pass it as check_run_complete_callback,
    # you'll get an error right away rather than going into an infinite loop
    return True


def basic_clean_up_f(f):
    """ Removes list of files in f

        f: file containing list of filepaths

        example f:
         f1.txt
         f2.txt
         f3.txt
         temp_dir

        If f contains the four lines above, this function would delete
         those three files/directories.

    """
    deletion_list = [l.strip() for l in f]
    remove_all(deletion_list)
    return True


def verbose_check_run_complete_f(f):
    """ Return True if all filepaths exist

        f: file containing list of filepaths

        example f:
         f1.txt
         f2.txt
         f3.txt

        If f contains the three lines above, this function would return
         False if any of these three files did not exist, and True otherwise.

    """
    filepaths = [l.strip() for l in f]
    for fp in filepaths:
        if not exists(fp):
            print "At least one fp doesn't exist: %s" % fp
            return False
    print "All filepaths exist."
    return True


def verbose_process_run_results_f(f):
    """ Copy each list of infiles to each outfile and delete infiles (verbose)

        f: file containing one set of mapping instructions per line

        example f:
         f1.txt f2.txt f3.txt f_combined.txt
         f1.log f2.log f3.log f_combined.log

        If f contained the two lines above, this function would
         concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
         and f1.log, f2.log, and f3.log into f_combined.log
    """
    infiles_lists, out_filepaths = parse_tmp_to_final_filepath_map_file(f)
    for infiles_list, out_filepath in zip(infiles_lists, out_filepaths):
        try:
            of = open(out_filepath, 'w')
            print 'Final result file (%s) contains temp files:'\
                % out_filepath
        except IOError:
            raise IOError("Poller can't open final output file: %s" % out_filepath +
                          "\nLeaving individual jobs output.\n Do you have write access?")

        for fp in infiles_list:
            print '\t%s' % fp
            for line in open(fp):
                of.write(line)
        of.close()
    return True


def verbose_clean_up_f(f):
    """ Removes list of files in f (verbose)

        f: file containing list of filepaths

        example f:
         f1.txt
         f2.txt
         f3.txt
         temp_dir

        If f contains the four lines above, this function would delete
         those three files/directories.

    """
    deletion_list = [l.strip() for l in f]
    remove_all(deletion_list)
    print "Post-run clean-up complete."
    return True


def poller(check_run_complete_f,
           process_run_results_f,
           clean_up_f,
           check_run_complete_file,
           process_run_results_file,
           clean_up_file,
           seconds_to_sleep):
    """ Polls for completion of job(s) and then processes/cleans up results

        check_run_complete_f: function which returns True when polled
         job(s) complete and False otherwise
        process_run_results_f: function applied to process the results
         of the polled job(s) -- run only after check_run_complete_f => True
        clean_up_f: function applied to clean up after the polled
         job(s) -- run after process_run_results_f
        check_run_complete_file: file passed to check_run_complete_f
         on each call
        process_run_results_file: file passed to process_run_results_f
        clean_up_file: file passed to clean_up_f
        seconds_to_sleep: number of seconds to sleep between calls
         to check_run_complete_f

    """
    number_of_loops = 0
    while(not check_run_complete_f(check_run_complete_file)):
        sleep(seconds_to_sleep)
        number_of_loops += 1
    process_run_results_f(process_run_results_file)
    clean_up_f(clean_up_file)
    est_per_proc_run_time = number_of_loops * seconds_to_sleep
    return est_per_proc_run_time
