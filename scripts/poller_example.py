#!/usr/bin/env python
# File created on 08 Nov 2009.
from __future__ import division
from os import system, makedirs
from qiime.parallel.poller import remove_all
from qiime.util import load_qiime_config
from os.path import exists
from optparse import OptionParser

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

qiime_config = load_qiime_config()

usage_str = """usage: %prog [options] {-d POLLED_DIR}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

poller.py is designed for use with parallel jobs to wait for their
 completion, and subsequently process the results and clean up. This 
 script allows users to see it in action, and also to allow manual testing
 as this is a difficult process to unit test.
 
To test, call the example command below. The poller will begin running,
 at which time you can create the three polled files in POLLED_DIR. When
 all three are created, the poller will process the results, clean up, and exit.

Example usage:
 The poller waits for three files to be created:
  <POLLED_DIR>/poller_test_0.txt
  <POLLED_DIR>/poller_test_1.txt
  <POLLED_DIR>/poller_test_2.txt
  
  <POLLED_DIR> is defined via -d.
  
 Existence of these three files is checked every 5 seconds with
 verbose_check_run_complete_f. When all three exist verbose_process_run_results_f 
 is called, which cats all the files into a single file:
  <POLLED_DIR>/poller_test_completed.txt
 Finally, verbose_clean_up_f is called which removes the original three 
 files the poller was waiting on.
  
  python qiime/scripts/poller_example.py -d /Users/caporaso/poller_test/
  
 The actual call to the polling command is printed for reference just prior 
  to calling it. This illustrates how to pass both functions and filepaths to
  the poller. For an example where the default (non-verbose) check_run_complete_f,
  process_run_results_f, and clean_up_f are used, pass -c. Again, the polling
  command will be printed just prior to calling:
  
  python qiime/scripts/poller_example.py -d /Users/caporaso/poller_test/ -c
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)
    
    parser.add_option('-d','--polled_dir',\
           type='string',help='path to directory to poll '+\
           ' [REQUIRED]')
    
    parser.add_option('-c','--suppress_custom_functions',\
           action='store_true',help='use the default functions for '+\
           'checking run completion, processing results, and '+\
           'cleaning up (these are quiet) [defualt: %default]')
    
    parser.add_option('-P','--poller_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/parallel/poller.py [default: %default]')

    parser.add_option('-Y','--python_exe_fp',
           help='full path to python executable [default: %default]',\
           default=qiime_config['python_exe_fp'])

    # Set default values here if they should be other than None
    parser.set_defaults(poller_fp=qiime_config['poller_fp'],\
     suppress_custom_functions=False)

    opts,args = parser.parse_args()
    required_options = ['polled_dir']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args

def write_poller_files(polled_dir):
    """ write files to support the poller
    """
    try:
        makedirs(polled_dir)
    except OSError:
        pass
    
    check_run_complete_fp = '%s/check_run_complete.txt' % polled_dir
    process_run_results_fp = '%s/process_run_results.txt' % polled_dir
    clean_up_fp = '%s/clean_up.txt' % polled_dir
    
    # Define files to be polled for existence:
    #   <polled_dir>/poller_test_0.txt
    #   <polled_dir>/poller_test_1.txt
    #   <polled_dir>/poller_test_2.txt
    example_filepaths = [polled_dir + '/poller_test_%d.txt'\
          % i for i in range(3)]
    open(check_run_complete_fp,'w').write('\n'.join(example_filepaths))

    # Defining map for merging of original files:
    #   <polled_dir>/poller_test_0.txt
    #   <polled_dir>/poller_test_1.txt
    #   <polled_dir>/poller_test_2.txt
    # into a single output file:
    #   $HOME/poller_test/poller_test_completed.txt
    m = example_filepaths + ['%s/poller_test_completed.txt' % polled_dir]
    open(process_run_results_fp,'w').write('\t'.join(m))
 
    # Defining files to be removed:
    #   <polled_dir>/poller_test_0.txt
    #   <polled_dir>/poller_test_1.txt
    #   <polled_dir>/poller_test_2.txt
    open(clean_up_fp,'w').write('\n'.join(example_filepaths))
    
    return check_run_complete_fp, process_run_results_fp, clean_up_fp

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    poller_fp = opts.poller_fp
    python_exe_fp = opts.python_exe_fp
    polled_dir = opts.polled_dir
    suppress_custom_functions = opts.suppress_custom_functions
    
    check_run_complete_fp, process_run_results_fp, clean_up_fp = \
     write_poller_files(polled_dir)
    
    if not suppress_custom_functions:
        print 'Polling directory:\n %s' % polled_dir
        command = '%s %s -r %s -f %s -p %s -m %s -c %s -d %s -t 5' %\
         (python_exe_fp, 
          poller_fp,
          'qiime.parallel.poller.verbose_check_run_complete_f',
          check_run_complete_fp,
          'qiime.parallel.poller.verbose_process_run_results_f',
          process_run_results_fp,
          'qiime.parallel.poller.verbose_clean_up_f',
          clean_up_fp)
        print 'Polling command:\n %s' % command
        system(command)
    else:
        print 'Polling directory:\n %s' % polled_dir
        command = '%s %s -f %s -m %s -d %s -t 5' %\
         (python_exe_fp, 
          poller_fp,
          check_run_complete_fp,
          process_run_results_fp,
          clean_up_fp)
        print 'Polling command:\n %s' % command
        system(command)
          
          
          