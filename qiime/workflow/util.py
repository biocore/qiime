#!/usr/bin/env python
# File created on 30 Dec 2009.
from __future__ import division
import sys
import re
from subprocess import Popen, PIPE, STDOUT
from os import makedirs, listdir
from glob import glob
from os.path import split, splitext, join, dirname, abspath
from datetime import datetime
from numpy import array
from cogent.util.misc import safe_md5
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import parse_mapping_file, parse_qiime_parameters
from qiime.util import (compute_seqs_per_library_stats,
                        get_qiime_scripts_dir,
                        create_dir, guess_even_sampling_depth,
                        get_interesting_mapping_fields,qiime_system_call,
                        get_qiime_library_version)
from biom.parse import parse_biom_table
from cogent.core.moltype import IUPAC_DNA_ambiguities
import os

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

"""
This file contains the QIIME workflow functions which string together 
independent scripts. For usage examples see the related files in the 
scripts directory:
 - 
"""

## Start utilities used by the workflow functions
def generate_log_fp(output_dir,
                    basefile_name='log',
                    suffix='txt',
                    timestamp_pattern='%Y%m%d%H%M%S'):
    timestamp = datetime.now().strftime(timestamp_pattern)
    filename = '%s_%s.%s' % (basefile_name,timestamp,suffix)
    return join(output_dir,filename)

class WorkflowError(Exception):
    pass

class WorkflowLogger(object):
    
    def __init__(self,log_fp=None,params=None,qiime_config=None,open_mode='w'):
        if log_fp:
            self._f = open(log_fp,open_mode)
        else:
            self._f = None
        start_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('Logging started at %s\n' % start_time)
        self.write('QIIME version: %s\n\n' % get_qiime_library_version())
        self.writeQiimeConfig(qiime_config)
        self.writeParams(params)
    
    def write(self,s):
        if self._f:
            self._f.write(s)
            # Flush here so users can see what step they're
            # on after each write, since some steps can take
            # a long time, and a relatively small amount of 
            # data is being written to the log files.
            self._f.flush()
        else:
            pass
    
    def writeQiimeConfig(self,qiime_config):
        if qiime_config == None:
            self.write('No qiime config provided.\n')
        else:
            self.write('qiime_config values:\n')
            for k,v in qiime_config.items():
                if v:
                    self.write('%s\t%s\n' % (k,v))
            self.write('\n')
            
    def writeParams(self,params):
        if params == None:
            self.write('No params provided.\n')
        else:
            self.write('parameter file values:\n')
            for k,v in params.items():
                for inner_k,inner_v in v.items():
                    val = inner_v or 'True'
                    self.write('%s:%s\t%s\n' % (k,inner_k,val))
            self.write('\n')
    
    def close(self):
        end_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('\nLogging stopped at %s\n' % end_time)
        if self._f:
            self._f.close()
        else:
            pass

def print_commands(commands,
                   status_update_callback,
                   logger,
                   close_logger_on_success=True):
    """Print list of commands to run """
    logger.write("Printing commands only.\n\n")
    for c in commands:
        for e in c:
            status_update_callback('#%s' % e[0])
            print '%s' % e[1]
            logger.write('# %s command\n%s\n\n' % e)
            
def call_commands_serially(commands,
                           status_update_callback,
                           logger,
                           close_logger_on_success=True):
    """Run list of commands, one after another """
    logger.write("Executing commands.\n\n")
    for c in commands:
        for e in c:
            status_update_callback('%s\n%s' % e)
            logger.write('# %s command \n%s\n\n' % e)
            stdout, stderr, return_value = qiime_system_call(e[1])
            if return_value != 0:
                msg = "\n\n*** ERROR RAISED DURING STEP: %s\n" % e[0] +\
                 "Command run was:\n %s\n" % e[1] +\
                 "Command returned exit status: %d\n" % return_value +\
                 "Stdout:\n%s\nStderr\n%s\n" % (stdout,stderr)
                logger.write(msg)
                logger.close()
                raise WorkflowError, msg
            # in the no error case, we write commands' output to the log
            # and also echo to this proc's stdout/stderr
            else:
                # write stdout and stderr to log file
                logger.write("Stdout:\n%s\nStderr:\n%s\n" % (stdout,stderr))
                # write stdout to stdout
                if stdout:
                    print stdout
                # write stderr to stderr
                if stderr:
                    sys.stderr.write(stderr)
    if close_logger_on_success: logger.close()

def print_to_stdout(s):
    print s
    
def no_status_updates(s):
    pass

def get_params_str(params):
    result = []
    for param_id, param_value in params.items():
        result.append('--%s' % (param_id))
        if param_value != None:
            result.append(param_value)
    return ' '.join(result)

def validate_and_set_jobs_to_start(params,
                                   jobs_to_start,
                                   default_jobs_to_start,
                                   parallel,
                                   option_parser):
    if (jobs_to_start != int(default_jobs_to_start)) and \
       not parallel:
        option_parser.error("Passing -O requires that -a is also passed.")
    params['parallel']['jobs_to_start'] = str(jobs_to_start)

def log_input_md5s(logger,fps):
    logger.write("Input file md5 sums:\n")
    for fp in fps:
        if fp != None:
            logger.write("%s: %s\n" % (fp, safe_md5(open(fp)).hexdigest()))
    logger.write("\n")
    

## End utilities used by the workflow functions

## Begin task-specific workflow functions

    
    
## Start reference otu picking workflow


## Begin task-specific workflow functions


def run_ampliconnoise(mapping_fp,
    output_dir, command_handler, params, qiime_config,
    logger=None, status_update_callback=print_to_stdout,
    chimera_alpha=-3.8228,chimera_beta=0.6200, sff_txt_fp=None, numnodes=2,
    suppress_perseus=True, output_filepath=None, platform='flx',
    seqnoise_resolution=None, truncate_len=None):
    """ Run the ampliconnoise pipeline
    
        The steps performed by this function are:
1. Split input sff.txt file into one file per sample

2. Run scripts required for PyroNoise

3. Run scripts required for SeqNoise

4. Run scripts requred for Perseus (chimera removal)

5. Merge output files into one file similar to the output of split_libraries.py

    output_filepath should be absolute
    seqnoise_resolution should be string
    environment variable PYRO_LOOKUP_FILE must be set correctly. Thus be
    careful passing command handlers that don't spawn child processes, as they
    may not inherit the correct environment variable setting
    """
    map_data,headers,comments = parse_mapping_file(open(mapping_fp,'U'))
    create_dir(output_dir)

    if seqnoise_resolution == None:
        if platform=='flx': seqnoise_resolution = '30.0'
        elif platform=='titanium': seqnoise_resolution = '25.0'
        else: raise RuntimeError('seqnoise_resolution not set, and no'+\
            ' default for platform '+platform)

    if truncate_len == None:
        if platform=='flx': truncate_len = '220'
        elif platform=='titanium': truncate_len = '400'
        else: raise RuntimeError('truncate_len not set, and no'+\
            ' default for platform '+platform)

    sample_names = [] # these are filenames minus extension, and are sample IDs
    primer_seqs = [] # same order as sample_names
    bc_seqs = [] # same order as sample_names
    for i in range(len(map_data)):
        sample_names.append(map_data[i][headers.index('SampleID')])
        bc_seqs.append(map_data[i][headers.index('BarcodeSequence')])
        # don't know why don't just take off the primer now. 
        # but that's done later
        # primer += (map_data[i][headers.index('LinkerPrimerSequence')])
        # for char, bases in IUPAC_DNA_ambiguities.items():
        #     primer = primer.replace(char,'['+''.join(bases)+']')

        primer = (map_data[i][headers.index('LinkerPrimerSequence')])
        for char, bases in IUPAC_DNA_ambiguities.items():
            primer = primer.replace(char,'['+''.join(bases)+']')
        primer_seqs.append(primer)

    if len(set(primer_seqs)) != 1:
        raise RuntimeError(
            'Error: only one primer per mapping file supported.')
    one_primer = primer_seqs[0]

    commands = []
    python_exe_fp = qiime_config['python_exe_fp']
    script_dir = get_qiime_scripts_dir()

    if logger == None:
        logger = WorkflowLogger(generate_log_fp(output_dir),
                                params=params,
                                qiime_config=qiime_config)
        close_logger_on_success = True
    else:
        close_logger_on_success = False
    log_input_md5s(logger,[mapping_fp,sff_txt_fp])

    # execute commands in output_dir
    called_dir = os.getcwd()
    os.chdir(output_dir)
    fh = open(os.path.join(output_dir,'map.csv'),'w')
    for i in range(len(sample_names)):
        fh.write(sample_names[i]+','+bc_seqs[i]+'\n')
    fh.close()

    # these are the fasta results, e.g. PC.636_Good.fa
    # later we merge them and copy to output file
    post_pyro_tail = '_'+truncate_len
    if suppress_perseus == True:
        fasta_result_names = [sample_name + post_pyro_tail+'_seqnoise_cd.fa'
          for sample_name in sample_names]
    else:
        fasta_result_names = [sample_name + '_Good.fa' \
          for sample_name in sample_names]

    cmd = 'cd '+output_dir # see also os.chdir above
    commands.append([('change to output dir', cmd)])
    
    cmd = 'echo $PYRO_LOOKUP_FILE > pyro_lookup_filepath.txt'
    commands.append([('confirm pyro lookup filepath environment variable',
        cmd)])


    cmd = 'SplitKeys.pl '+one_primer+' map.csv < '+\
        os.path.join(called_dir,sff_txt_fp)+\
        ' > splitkeys_log.txt 2> unassigned.fna'
    commands.append([('split sff.txt via barcodes (keys)', cmd)])

    for i, sample_name in enumerate(sample_names):

        # Build the summarize taxonomy command
        if platform == 'flx':
            cmd = 'Clean360.pl '+one_primer+' '+sample_name+' < '+\
                sample_name+'.raw'
            commands.append([('clean flows '+sample_name, cmd)])

            # these run through the whole sff file once per sample, I think
            # cmd = "FlowsFA.pl " + primer_seqs[i] + ' '+sample_name +' < '+\
            #     os.path.join(called_dir,sff_txt_fp)
            # commands.append([('extract flows '+sample_name, cmd)])
        elif platform == 'titanium':
            cmd = 'CleanMinMax.pl '+one_primer+' '+sample_name+' < '+\
                sample_name+'.raw'
            commands.append([('clean flows '+sample_name, cmd)])

            # cmd = "FlowsMinMax.pl " + primer_seqs[i] + ' '+sample_name +' < '+\
            #     os.path.join(called_dir,sff_txt_fp)
            # commands.append([('extract flows '+sample_name, cmd)])
        else:
            raise RuntimeError("platform " + platform + " not supported")

        cmd = "mpirun -np "+str(numnodes)+" PyroDist -in "+\
          sample_name+".dat -out "+sample_name+ " > "+sample_name+".pdout"
        commands.append([('pyrodist '+sample_name, cmd)])

        cmd = "FCluster -in "+sample_name+".fdist -out "+sample_name+\
          " > "+sample_name+".fcout"
        commands.append([('fcluster pyrodist '+sample_name, cmd)])

# e.g.:
# mpirun -np 2 PyroNoise -din PC.354.dat -out PC.354_pyronoise -lin
# PC.354.list -s 60.0 -c 0.01 > PC.354_pyronoise.pnout
        cmd = "mpirun -np "+str(numnodes)+" PyroNoise -din "+\
            sample_name+".dat -out "+\
            sample_name+"_pyronoise "+"-lin "+\
            sample_name+".list -s 60.0 -c 0.01 > "+\
            sample_name+"_pyronoise.pnout"
        commands.append([('pyronoise '+sample_name, cmd)])

        cmd = 'Parse.pl '+bc_seqs[i]+one_primer+' '+truncate_len+' < '+\
            sample_name+'_pyronoise_cd.fa'+' > '+ sample_name+'_'+\
            truncate_len+'.fa'
        commands.append([('truncate '+sample_name, cmd)])

        # now start with post_pyro_tail
        cmd = "mpirun -np "+str(numnodes)+" SeqDist -in "+\
            sample_name+post_pyro_tail+\
            ".fa > "+sample_name+post_pyro_tail+".seqdist"
        commands.append([('seqdist '+sample_name, cmd)])

        cmd = "FCluster -in "+sample_name+post_pyro_tail+".seqdist -out "+\
            sample_name+post_pyro_tail+"fcl > "+\
            sample_name+post_pyro_tail+".fcout"
        commands.append([('fcluster seqdist '+sample_name, cmd)])

# e.g.:
# mpirun -np 2 SeqNoise -in PC.354_pyronoise_cd.fa -din
# PC.354_pyronoise_cd.seqdist -out PC.354_pyronoise_cd_seqnoise -lin
# PC.354_pyronoise_cdfcl.list -min PC.354_pyronoise.mapping -s 30.0 -c 0.08 >
# PC.354_pyronoise_cd.snout

        cmd = "mpirun -np "+str(numnodes)+" SeqNoise -in "+\
            sample_name+post_pyro_tail+\
            ".fa -din "+sample_name+post_pyro_tail+".seqdist -out "+\
            sample_name+post_pyro_tail+\
            "_seqnoise -lin "+sample_name+post_pyro_tail+'fcl.list -min '+\
            sample_name+'_pyronoise'+\
            '.mapping -s '+seqnoise_resolution+' -c 0.08 > '+\
            sample_name+post_pyro_tail+'.snout'
        commands.append([('seqnoise '+sample_name, cmd)])


        if suppress_perseus == False: 

            cmd = 'Perseus -sin '+sample_name+post_pyro_tail+\
                '_seqnoise_cd.fa > ' +\
                sample_name+'.per'
            commands.append([('Perseus '+sample_name, cmd)])

            cmd = 'Class.pl '+sample_name+'.per '+\
                str(chimera_alpha) + ' '+ str(chimera_beta)+\
                ' > '+sample_name+'.class'
            commands.append([('Class.pl '+sample_name, cmd)])

            cmd = 'FilterGoodClass.pl '+sample_name+post_pyro_tail+\
                '_seqnoise_cd.fa '+\
                sample_name+'.class 0.5 > '+sample_name+'_Chi.fa 2> '+\
                sample_name+'_Good.fa'
            commands.append([('FilterGoodClass '+sample_name, cmd)])

        cmd = '%s %s/unweight_fasta.py -i %s -o %s -l %s' %\
         (python_exe_fp, script_dir, fasta_result_names[i], 
         sample_name+'_unw.fna', sample_name)
        commands.append([('unweight fasta '+sample_name, cmd)])

    cmd = 'cat ' +\
      ' '.join([sample_name+'_unw.fna' for sample_name in sample_names]) +\
      ' > ' + output_filepath # this should be an abs filepath
    commands.append([('cat into one fasta file', cmd)])

    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback,
                    logger=logger,
                    close_logger_on_success=close_logger_on_success)

