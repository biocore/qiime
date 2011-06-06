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

# qsub template
#requires format string (walltime, ncpus, nodes, queue, job_name, keep_output, command)
QSUB_TEXT = """# Walltime Limit: hh:nn:ss 
#PBS -l walltime=%s 

# Node Specification:
#PBS -l ncpus=%d -l nodes=%d

# Queue: Defaults to friendlyq 
#PBS -q %s 

# Mail: options are (a) aborted, (b) begins execution, (e) ends execution
# use -M <email> for additional recipients
# supress email notification
#PBS -m n

# Job Name:
#PBS -N %s 

# Keep output
#PBS -k %s

echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------
cd $PBS_O_WORKDIR 
%s
""" 

def make_jobs(commands, job_prefix, queue, jobs_dir="jobs/",
              walltime="72:00:00", ncpus=1, nodes=1, keep_output="oe"):
    """prepare qsub text files.
    
    command: list of commands
    
    job_prefix: a short, descriptive name for the job.

    queue: name of the queue to submit to
    
    jobs_dir: path to directory where job submision scripts are written

    walltime: the maximal walltime 
    
    ncpus: number of cpus
    
    nodes: number of nodes
    
    keep_output: keep standard error, standard out, both, or neither
                 o=std out, e=std err, oe=both, n=neither
    """

    filenames=[]
    for command in commands:
        job_name = get_tmp_filename(tmp_dir=jobs_dir, prefix=job_prefix+"_",
                                    suffix = ".txt")
        out_fh = open(job_name,"w")

        out_fh.write(QSUB_TEXT % (walltime, ncpus, nodes, queue, job_prefix,
                                  keep_output, command))        
        out_fh.close()
        filenames.append(job_name)
    return filenames

def submit_jobs(filenames, verbose=False):
    """Submit jobs in filenames.

    filenames: list of prepared qsub job scripts, ready to be submitted

    verbose: a binary verbose flag
    """
    if(not app_path("qsub")):
        raise ApplicationNotFoundError,"qsub not found. Can't submit jobs."
    
    for file in filenames:        
        command = 'qsub %s' % file
        result = Popen(command, shell=True, universal_newlines=True,\
                           stdout=PIPE, stderr=STDOUT).stdout.read()
        if verbose:
            print result

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
