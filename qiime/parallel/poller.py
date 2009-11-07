#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# poller.py

""" Description
File created on 30 Sep 2009.

"""
from __future__ import division
from time import sleep
from optparse import OptionParser
from os import getenv, remove
from os.path import exists, isdir
from shutil import rmtree
from cogent.util.misc import remove_files

usage_str = """usage: %prog [options] {-f RUN_OUTPUT_FILES}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
 See Qiime/scripts/poller_examples.py
 
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog 0.1'
    parser = OptionParser(usage=usage, version=version)

    # A binary 'verbose' flag
    parser.add_option('-v','--verbose',action='store_true',\
        help='Print information during execution -- '+\
        'useful for debugging [default: %default]')

    parser.add_option('-r','--check_run_complete_f',\
           help='function which returns True when run is completed '+\
           '[default: %default]',\
           default='qiime.parallel.poller.basic_check_run_complete_f')
    parser.add_option('-f','--check_run_complete_file',\
           help='path to file containing a list of files that must exist to' +\
           ' declare a run complete [REQUIRED]')
           
    parser.add_option('-p','--process_run_results_f',\
           help='function to be called when runs complete [default: %default]',\
           default='qiime.parallel.poller.basic_process_run_results_f')
    parser.add_option('-m','--process_run_results_file',\
           help='path to file containing a map of tmp filepaths which should' +\
           ' be written to final output filepaths [default: %default]')
           
    parser.add_option('-c','--clean_up_f',\
           help='function called after processing result [default: %default]',\
           default='qiime.parallel.poller.basic_clean_up_f')
    parser.add_option('-d','--clean_up_file',
           help='List of files and directories to remove after run'+\
           ' [default: %default]')
    
    parser.add_option('-t','--time_to_sleep',type='int',\
           help='time to wait between calls to status_callback_f'+\
           ' (in seconds) [default: %default]')

    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False,time_to_sleep=60)

    opts,args = parser.parse_args()
    
    if not opts.check_run_complete_f and not opts.run_output_file:
        parser.error('Must supply -r and/or -f.')

    return opts,args

def get_function_handle(s):
    last_dot = s.rindex('.')
    module_name = s[:last_dot]
    function_name = s[last_dot+1:]
    module = __import__(module_name,globals(),locals(),[function_name])
    function = eval('module.%s' % function_name)
    return function

def parse_filepath_list_file(lines):
    return [l.strip() for l in lines]

def parse_tmp_to_final_filepath_map_file(lines):
    infiles_lists = []
    out_filepaths = []
    for line in lines:
        fields = line.split()
        infiles_lists.append(fields[:-1])
        out_filepaths.append(fields[-1])
    return infiles_lists, out_filepaths

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
 
def basic_check_run_complete_f(f,verbose=False):
    """ Return a function which returns True if all filepaths exist
    
        f: file containing list of filepaths
        verbose: if True, the resulting function will print information to
         stdout (default:False)
    """
    filepaths = parse_filepath_list_file(f)
    for fp in filepaths:
        if not exists(fp):
            if verbose:
                print "At least one fp doesn't exist: %s" % fp
            return False
    if verbose: print "All filepaths exist."
    return True
    
def basic_process_run_results_f(f,verbose=False):
    """ Copy each list of infiles to each outfile and delete infiles
    """
    infiles_lists,out_filepaths = parse_tmp_to_final_filepath_map_file(f)
    for infiles_list, out_filepath in zip(infiles_lists,out_filepaths):
        try:
            of = open(out_filepath,'w')
        except IOError:
            raise IOError,\
             "Poller can't open final output file: %s" % out_filepath  +\
             "\nLeaving individual jobs output.\n Do you have write access?"

        for fp in infiles_list:
            for line in open(fp):
               of.write(line)
        of.close()

    if verbose: print "Post-run processing complete."
    # It is a good idea to have your clean_up_callback return True.
    # That way, if you get mixed up and pass it as check_run_complete_callback, 
    # you'll get an error right away rather than going into an infinite loop
    return True
    
def basic_clean_up_f(f,verbose=False):
    deletion_list = parse_filepath_list_file(f)
    remove_all(deletion_list)
    if verbose: 
        print "Post-run clean-up complete."
    return True
 
    
def poller(check_run_complete_f,\
            process_run_results_f,\
            clean_up_f,\
            check_run_complete_file,\
            process_run_results_file,\
            clean_up_file,\
            seconds_to_sleep,\
            verbose=False,
            log_filepath=None):
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
        verbose: turn on verbose output [default: False]
        log_filepath: path to store log file [default: None]
        
    """
    number_of_loops = 0
    while(not check_run_complete_f(check_run_complete_file,verbose)):
        sleep(seconds_to_sleep)
        number_of_loops += 1
    process_run_results_f(process_run_results_file,verbose)
    clean_up_f(clean_up_file,verbose)
    est_per_proc_run_time = number_of_loops * seconds_to_sleep
    return est_per_proc_run_time



if __name__ == "__main__":
    opts,args = parse_command_line_parameters()

    poller(get_function_handle(opts.check_run_complete_f),\
            get_function_handle(opts.process_run_results_f),\
            get_function_handle(opts.clean_up_f),\
            list(open(opts.check_run_complete_file)),\
            list(open(opts.process_run_results_file)),\
            list(open(opts.clean_up_file)),\
            seconds_to_sleep=opts.time_to_sleep,\
            verbose=opts.verbose,\
            log_filepath=None)


  