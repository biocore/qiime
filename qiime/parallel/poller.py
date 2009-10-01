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

usage_str = """usage: %prog [options] {(-r CHECK_RUN_COMPLETE_F | -f RUN_OUTPUT_FILES) (-c CLEAN_UP_CALLBACK_F | -m TMP_TO_FINAL_FILEPATH_MAP)}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
 Call the poller with the example_check_run_complete_callback (-r) function
 which waits for three files to be created:
  $HOME/poller_test/poller_test_0.txt
  $HOME/poller_test/poller_test_1.txt
  $HOME/poller_test/poller_test_2.txt
  
 Existence of these three files is checked every 5 seconds (-t), and when
 they all exist qiime.parallel.poller.example_clean_up_callback (-c) is called.
 This clean-up function cats all the files into a single file:
  $HOME/poller_test/poller_test_completed.txt
  
 and removes the original three files it was waiting on.
  
 python Qiime/qiime/parallel/poller.py -r qiime.parallel.poller.example_check_run_complete_callback -p qiime.parallel.poller.example_clean_up_callback -t 5
 
 poller.py is designed for use on the cluster to clean up parallel jobs. However,
 you can see it in action by running the above command and in a new terminal 
 creating the three files above in $HOME/poller_test. Put some text in each file
 and once all three exist, the poller will join them into a single new file and 
 remove the originals.
 
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
           '[default: %default]')
    parser.add_option('-f','--run_output_file',\
           help='path to file containing a list of files that must exist to' +\
           ' declare a run complete [default: %default]')
           
    parser.add_option('-p','--process_run_results_f',\
           help='function to be called when runs complete [default: %default]')
    parser.add_option('-m','--tmp_to_final_filepath_map',\
           help='path to file containing a map of tmp filepaths which should' +\
           ' be written to final output filepaths [default: %default]')
           
    parser.add_option('-c','--clean_up_f',\
           help='function to be called after processing result [default: %default]')
    parser.add_option('-d','--deletion_list',
           help='List of files and directories to remove after run'+\
           ' [default: %default]')
    
    parser.add_option('-t','--time_to_sleep',type='int',\
           help='time to wait between calls to status_callback_f'+\
           ' (in seconds) [default: %default]')

    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False,time_to_sleep=60)

    opts,args = parser.parse_args()
    
    if not opts.check_run_complete_f and not opts.run_output_file:
        parser.error('Must supply either -r or -f.')
    if not opts.process_run_results_f and not opts.tmp_to_final_filepath_map:
        parser.error('Must supply either -p or -m.')
    # if not opts.clean_up_f and not opts.deletion_list:
    #     parser.error('Must supply either -c or -d.')

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

# # Defining an example check_run_complete_f which returns True
# # if the following files all exist:
# #   $HOME/poller_test/poller_test_0.txt
# #   $HOME/poller_test/poller_test_1.txt
# #   $HOME/poller_test/poller_test_2.txt
# home_dir = getenv('HOME')
# example_filepaths = [home_dir + '/poller_test/poller_test_%d.txt'\
#       % i for i in range(3)]
# example_check_run_complete_callback = \
#  get_basic_check_run_complete_callback(example_filepaths,True)
# 
# ## Defining an example process_run_results_f
# # Merges the text from the three 'original' files:
# #   $HOME/poller_test/poller_test_0.txt
# #   $HOME/poller_test/poller_test_1.txt
# #   $HOME/poller_test/poller_test_2.txt
# # into a single output file:
# #   $HOME/poller_test/poller_test_completed.txt      
# # and removes the original files.
# example_out_filepath = home_dir + '/poller_test/poller_test_completed.txt'
# example_clean_up_callback = \
#  get_basic_clean_up_callback([example_filepaths],[example_out_filepath],True)

def get_basic_check_run_complete_f(filepaths,verbose=False):
    """ Return a function which returns True if all filepaths exist
    
        filepaths: list of file paths which should be checked for existence
        verbose: if True, the resulting function will print information to
         stdout (default:False)
    """
    def result():
        for fp in filepaths:
            if not exists(fp):
                if verbose:
                    print "At least one fp doesn't exist: %s" % fp
                return False
        if verbose: print "All filepaths exist."
        return True
    
    return result

def get_basic_process_run_results_f(infiles_lists,out_filepaths,verbose=False):
    """ Copy each list of infiles to each outfile and delete infiles
    """
    def result():
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
    
        if verbose: print "Clean up completed."
        # It is a good idea to have your clean_up_callback return True.
        # That way, if you get mixed up and pass it as check_run_complete_callback, 
        # you'll get an error right away rather than going into an infinite loop
        return True
    return result
    
def get_basic_clean_up_f(deletion_list):
    def result():
        remove_all(deletion_list)
        return True
    return result
    
def poller(check_run_complete_f,process_run_results_f,clean_up_f,\
           seconds_to_sleep=60,log_filepath=None):
    """ The core poller function
    """
    number_of_loops = 0
    while(not check_run_complete_f()):
        sleep(seconds_to_sleep)
        number_of_loops += 1
    process_run_results_f()
    clean_up_f()
    est_per_proc_run_time = number_of_loops * seconds_to_sleep
    return est_per_proc_run_time

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    seconds_to_sleep = opts.time_to_sleep
    
    if opts.check_run_complete_f:
        check_run_complete_f = \
         get_function_handle(opts.check_run_complete_f)
    else:
        output_filepaths = \
         parse_filepath_list_file(open(opts.run_output_file))
        check_run_complete_f = \
         get_basic_check_run_complete_f(output_filepaths,verbose)
    
    if opts.process_run_results_f:
        process_run_results_f = \
         get_function_handle(opts.process_run_results_f)
    else:
        infiles_lists, out_filepaths = \
         parse_tmp_to_final_filepath_map_file(\
         open(opts.tmp_to_final_filepath_map))
        process_run_results_f = \
         get_basic_process_run_results_f(infiles_lists,out_filepaths,verbose)
    
    if opts.clean_up_f:
        clean_up_f = \
         get_function_handle(opts.clean_up_f)
    elif opts.deletion_list:
        deletion_list = \
         parse_filepath_list_file(open(opts.deletion_list))
        clean_up_f = \
         get_basic_clean_up_f(deletion_list)
    else:
        clean_up_f = \
         get_basic_clean_up_f([])
    
    poller(check_run_complete_f,process_run_results_f,clean_up_f,\
     seconds_to_sleep=seconds_to_sleep,log_filepath=None)
    
    
    
    
    