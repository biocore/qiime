#!/usr/bin/env python 

"""A simple qsub based cluster submission script."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Jens Reeder", "Rob Knight"]#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.2.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"

from optparse import OptionParser
from os.path import exists
from os import remove, rename, rmdir, makedirs
from subprocess import Popen, PIPE, STDOUT

from cogent.util.misc import app_path
from cogent.app.util import get_tmp_filename, ApplicationNotFoundError

from qiime.denoiser.make_cluster_jobs import make_jobs, submit_jobs

# Does not use qiime option parsers as we might want to keep it
# independent from qiime
def parse_command_line_parameters(commandline_args=None):
    """ Parses command line arguments """
    usage   = 'usage: %prog [options] <commands file> <job prefix>'
    version = 'Version: %prog 0.1'
    parser  = OptionParser(usage=usage, version=version)

    # A binary 'verbose' flag
    parser.add_option('-v','--verbose',action='store_true', dest='verbose',\
                          help='Print information during execution -- '+\
                          'useful for debugging [default: %default]')

    parser.add_option('-m','--make_jobs',action='store_true',\
                          dest='make_jobs', help='Prepare qsub text file'+\
                          ' [default: %default]')

    parser.add_option('-s','--submit_jobs',action='store_true',\
                          dest='submit_jobs', help='Prepare qsub text file'+\
                          ' [default: %default]')

    parser.add_option('-q','--queue',action='store',\
                          type='string',dest='queue', \
                          help='name of Queue to submit to '+\
                          ' [default: %default]')

    parser.add_option('-j','--job_dir', action='store',\
                          type='string',dest='job_dir',\
                          help='directory to store the jobs '+\
                          '[default: %default]')

    # Define defaults
    parser.set_defaults(verbose=False, make_jobs=True, submit_jobs=True,
                        job_dir="jobs/", queue='friendlyq')

    opts,args = parser.parse_args(commandline_args)
    
    if len(args)!= 2:
        parser.error('Program requires <commands file> and  <job prefix>')

    if (len(args[1])>10 or len(args[1])==0):
        parser.error('job prefix must be 1-10 characters long')

    return opts,args

def main(commandline_args=None):
    opts, args = parse_command_line_parameters(commandline_args)
 
    commands = list(open(args[0]))
    job_prefix = args[1]

    if(not exists(opts.job_dir)):
        try:
            makedirs(opts.job_dir)
        except OSError:
            exit(" Jobs directory can not be created. "
                 +"Check for permissions or file with the same name: %s\n"
                 % opts.job_dir)
                 
    if (opts.make_jobs):
        filenames = make_jobs(commands, job_prefix, opts.queue, opts.job_dir)
    else:
        exit("Should we ever get here???")
    if (opts.submit_jobs):
        submit_jobs(filenames, opts.verbose)

if __name__ == "__main__":
    main()
