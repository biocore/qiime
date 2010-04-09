#!/usr/bin/env python
#parallel_blast.py: make and run parallel blast given file of seqs and db

__author__ = "Rob Knight, Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from os.path import split, splitext
from qiime.parallel.util import get_rename_command

def get_commands(infile_paths,db_path,blast_executable_path,\
    blastmat_path,e_value,word_size,num_hits,output_dir,working_dir,\
    command_prefix=None,command_suffix=None):
    
    command_prefix = command_prefix or\
     '/bin/bash; export BLASTMAT=%s;' % blastmat_path
    command_suffix = command_suffix or\
     '; exit'
    
    commands = []
    result_filepaths = []
    
    for i, infile_path in enumerate(infile_paths):
        
        infile_basename = splitext(split(infile_path)[1])[0]
        working_outfile_path = '%s/%s_blast_out.txt' %\
          (working_dir,infile_basename)
        outfile_path = '%s/%s_blast_out.txt' % (output_dir,infile_basename)
        
        rename_command = '; mv %s %s' % (working_outfile_path, outfile_path)
        
        result_filepaths.append(outfile_path)
        
        command = \
         "%s %s -p blastn -m 9 -e %s -W %s -b %s -i %s -d %s > %s %s %s" % \
         (command_prefix,\
          blast_executable_path,\
          e_value,\
          word_size,\
          num_hits, 
          infile_path,\
          db_path,\
          working_outfile_path,\
          rename_command,\
          command_suffix)
        commands.append(command)
    
    return commands, result_filepaths

