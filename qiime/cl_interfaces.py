#!/usr/bin/env python
# File created on 31 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from os import makedirs
from qiime.util import (make_option,
                        load_qiime_config,
                        parse_command_line_parameters,
                        get_options_lookup)
from qiime.parse import parse_qiime_parameters
from qiime.workflow.util import (print_commands,
                            call_commands_serially,
                            print_to_stdout,
                            no_status_updates,
                            log_input_md5s,
                            WorkflowLogger,
                            generate_log_fp)
from qiime.workflow.upstream import run_qiime_data_preparation
from qiime.command.util import (WorkflowCommand,
                                  QiimeCommand)

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

class PickOtusThroughOtuTable(WorkflowCommand):
    """
    """
    _brief_description = """A workflow script for picking OTUs through building OTU tables"""
    _script_description = """This script takes a sequence file and performs all processing steps through building the OTU table."""
    _script_usage = [("""Simple example""","""The following command will start an analysis on seqs.fna (-i), which is a post-split_libraries fasta file. The sequence identifiers in this file should be of the form <sample_id>_<unique_seq_id>. The following steps, corresponding to the preliminary data preparation, are applied: Pick de novo OTUs at 97%; pick a representative sequence for each OTU (the OTU centroid sequence); align the representative set with PyNAST; assign taxonomy with RDP classifier; filter the alignment prior to tree building - remove positions which are all gaps, and specified as 0 in the lanemask; build a phylogenetic tree with FastTree; build an OTU table. All output files will be written to the directory specified by -o, and subdirectories as appropriate. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""","""%prog -i $PWD/seqs.fna -o $PWD/otus/""")]
    _script_usage_output_to_remove = ['$PWD/otus/']
    _output_description = """This script will produce an OTU mapping file (pick_otus.py), a representative set of sequences (FASTA file from pick_rep_set.py), a sequence alignment file (FASTA file from align_seqs.py), taxonomy assignment file (from assign_taxonomy.py), a filtered sequence alignment (from filter_alignment.py), a phylogenetic tree (Newick file from make_phylogeny.py) and a biom-formatted OTU table (from make_otu_table.py)."""
    _required_options = [
        make_option('-i','--input_fp',type='existing_filepath',
            help='the input fasta file [REQUIRED]'),
        make_option('-o','--output_dir',type='new_dirpath',
            help='the output directory [REQUIRED]'),
    ]
    _optional_options = [\
        make_option('-p','--parameter_fp',type='existing_filepath',
            help='path to the parameter file, which specifies changes'+\
                ' to the default behavior. '+\
                'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .'+\
                ' [if omitted, default values will be used]'),
        make_option('-f','--force',action='store_true',\
                dest='force',help='Force overwrite of existing output directory'+\
                ' (note: existing files in output_dir will not be removed)'+\
                ' [default: %default]'),\
        make_option('-w','--print_only',action='store_true',\
                dest='print_only',help='Print the commands but don\'t call them -- '+\
                'useful for debugging [default: %default]',default=False),\
        make_option('-a','--parallel',action='store_true',\
                dest='parallel',default=False,\
                help='Run in parallel where available [default: %default]'),
        options_lookup['jobs_to_start_workflow']
    ]
    _version = __version__
    
    _input_file_parameter_ids = ['input_fp','parameter_fp']

    def run_command(self, 
                    options,
                    arguments):
    
        verbose = options['verbose']
    
        input_fp = options['input_fp']
        output_dir = options['output_dir']
        verbose = options['verbose']
        print_only = options['print_only']
    
        parallel = options['parallel']
        # No longer checking that jobs_to_start > 2, but
        # commenting as we may change our minds about this.
        #if parallel: raise_error_on_parallel_unavailable()
    
        if options['parameter_fp']:
            try:
                parameter_f = open(options['parameter_fp'])
            except IOError:
                raise QiimeCommandError,\
                 "Can't open parameters file (%s). Does it exist? Do you have read access?"\
                 % options['parameter_fp']
            params = parse_qiime_parameters(parameter_f)
        else:
            params = parse_qiime_parameters([]) 
            # empty list returns empty defaultdict for now
            
        params['parallel']['jobs_to_start'] = self._validate_jobs_to_start(
                                                            options['jobs_to_start'],
                                                            qiime_config['jobs_to_start'],
                                                            parallel)
    
        try:
            makedirs(output_dir)
        except OSError:
            if options['force']:
                pass
            else:
                # Since the analysis can take quite a while, I put this check
                # in to help users avoid overwriting previous output.
                print "Output directory already exists. Please choose "+\
                 "a different directory, or force overwrite with -f."
                exit(1)
        
        if print_only:
            command_handler = print_commands
        else:
            command_handler = call_commands_serially
    
        if verbose:
            status_update_callback = print_to_stdout
        else:
            status_update_callback = no_status_updates
    
        run_qiime_data_preparation(
         input_fp, 
         output_dir,
         command_handler=command_handler,
         params=params,
         qiime_config=qiime_config,
         parallel=parallel,\
         status_update_callback=status_update_callback)
