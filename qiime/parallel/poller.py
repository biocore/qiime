#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# poller.py

""" Description
File created on 30 Sep 2009.

"""
from __future__ import division
from time import sleep
from optparse import OptionParser
from os import getenv
from os.path import exists
from cogent.util.misc import remove_files

usage_str = """usage: %prog [options] {-r CHECK_RUN_COMPLETE_F -c CLEAN_UP_CALLBACK_F}

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
  
 python Qiime/qiime/parallel/poller.py -r qiime.parallel.poller.example_check_run_complete_callback -c qiime.parallel.poller.example_clean_up_callback -t 5
 
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

    parser.add_option('-r','--check_run_complete_callback_f',\
           help='function which returns True when run is completed [REQUIRED]')
    parser.add_option('-c','--clean_up_callback_f',\
           help='function to be called when runs complete [REQUIRED]')
    parser.add_option('-t','--time_to_sleep',type='int',\
           help='time to wait between calls to status_callback_f'+\
           ' (in seconds) [default: %default]')

    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False,time_to_sleep=60)

    opts,args = parser.parse_args()
    required_options = ['check_run_complete_callback_f','clean_up_callback_f']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args

def example_check_run_complete_callback():
    """ An example check_run_complete_callback_f
        
        Returns True when three files all exist:
          $HOME/poller_test/poller_test_0.txt
          $HOME/poller_test/poller_test_1.txt
          $HOME/poller_test/poller_test_2.txt
    """
    home_dir = getenv('HOME')
    filepaths = [home_dir + '/poller_test/poller_test_%d.txt'\
      % i for i in range(3)]
    for fp in filepaths:
        if not exists(fp):
            print "At least one file does not yet exist: %s" % fp
            return False
    print "All files exist."
    return True
    
def example_clean_up_callback():
    """ An example clean_up_callback_f
    
        Merges the text from the three 'original' files:
          $HOME/poller_test/poller_test_0.txt
          $HOME/poller_test/poller_test_1.txt
          $HOME/poller_test/poller_test_2.txt
          
        into a single output file:
          $HOME/poller_test/poller_test_completed.txt
          
        and removes the original files.
        
    """
    home_dir = getenv('HOME')
    filepaths = [home_dir + '/poller_test/poller_test_%d.txt'\
     % i for i in range(3)]
    out_filepath = home_dir + '/poller_test/poller_test_completed.txt'
    
    try:
        of = open(out_filepath,'w')
    except IOError:
        raise IOError,\
         "Poller can't open final output file: %s" % out_filepath  +\
         "\nLeaving individual jobs output.\n Do you have write access?"
    
    for fp in filepaths:
        for line in open(fp):
           of.write(line)
    # clean up the individual job files
    remove_files(filepaths,error_on_missing=False)
    
    print "Clean up completed."
    # It is a good idea to have your clean_up_callback return True.
    # That way, if you get mixed up and pass it as check_run_complete_callback, 
    # you'll get an error right away rather than going into an infinite loop
    return True
    
def poller(check_run_complete_callback,clean_up_callback,\
           seconds_to_sleep=60,log_filepath=None):
    """ The core poller function
    """
    number_of_loops = 0
    while(not check_run_complete_callback()):
        sleep(seconds_to_sleep)
        number_of_loops += 1
    clean_up_callback()
    est_per_proc_run_time = number_of_loops * seconds_to_sleep
    return est_per_proc_run_time

def get_function_handle(s):
    last_dot = s.rindex('.')
    module_name = s[:last_dot]
    function_name = s[last_dot+1:]
    module = __import__(module_name,globals(),locals(),[function_name])
    function = eval('module.%s' % function_name)
    return function

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    check_run_complete_callback_f = \
     get_function_handle(opts.check_run_complete_callback_f)
    clean_up_callback_f = \
     get_function_handle(opts.clean_up_callback_f)
    seconds_to_sleep = opts.time_to_sleep
    
    poller(check_run_complete_callback_f,clean_up_callback_f,\
     seconds_to_sleep=seconds_to_sleep,log_filepath=None)
    
    
    
    
    