#!/usr/bin/env python
# File created on 02 Nov 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from os import makedirs
from qiime.util import (parse_command_line_parameters, 
                        make_option, 
                        get_options_lookup,
                        load_qiime_config,)
from qiime.parse import parse_qiime_parameters
from qiime.workflow import (validate_and_set_jobs_to_start, call_commands_serially,
                            print_commands, no_status_updates, print_to_stdout)
from qiime.pick_subsampled_reference_otus_through_otu_table import (
                        pick_subsampled_open_referenence_otus,
                        iterative_pick_subsampled_open_referenence_otus)

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()


script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""

script_info['script_usage'] = []

script_info['script_usage'].append(("","Run the subsampled open-reference OTU picking workflow on seqs1.fna using refseqs.fna as the reference collection. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/","%prog -i $PWD/seqs1.fna -r $PWD/refseqs.fna -o $PWD/ucrss/ -s 0.1 -p $PWD/ucrss_params.txt"))

script_info['script_usage'].append(("","Run the subsampled open-reference OTU picking workflow in iterative mode on seqs1.fna and seqs2.fna using refseqs.fna as the initial reference collection. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/","%prog -i $PWD/seqs1.fna,$PWD/seqs2.fna -r $PWD/refseqs.fna -o $PWD/ucrss_iter/ -s 0.1 -p $PWD/ucrss_params.txt"))

script_info['script_usage_output_to_remove'] = ['$PWD/ucrss/','$PWD/ucrss_iter/']

script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_fps',
             help='the input sequences filepath or comma-separated list of filepaths',
             type='existing_filepaths'),
 make_option('-r','--reference_fp',
             type='existing_filepath',
             help='the reference sequences'),
 make_option('-o','--output_dir',
             type='new_dirpath',
             help='the output directory'),
]
script_info['optional_options'] = [
 make_option('-p','--parameter_fp',
    help='path to the parameter file, which specifies changes'+\
        ' to the default behavior. '+\
        'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .'+\
        ' [if omitted, default values will be used]'),
 make_option('-n','--new_ref_set_id',default='New',
         help='Unique identifier for OTUs that get created in this ref set '+\
         '(this is useful to support combining of reference sets) [default:%default]'),
 make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),\
 ## print to shell script doesn't work for this workflow as there is a mix of 
 ## command calls and api calls. i need to refactor the workflow scripts...
 # make_option('-w','--print_only',action='store_true',\
 #        dest='print_only',help='Print the commands but don\'t call them -- '+\
 #        'useful for debugging [default: %default]',default=False),\
 make_option('-a','--parallel',action='store_true',\
        dest='parallel',default=False,\
        help='Run in parallel where available [default: %default]'),
 options_lookup['jobs_to_start_workflow'],
 make_option('-s','--percent_subsample',type='float',default='0.001',
        help='Percent of failure sequences to include in the subsample to '+\
        'cluster de novo (larger numbers should give more comprehensive ,'+\
        'results but will be slower) [default:%default]'),
 make_option('--prefilter_percent_id',type='float',default='0.60',
        help='Sequences are pre-clustered at this percent id against the reference'+\
              ' and any reads which fail to hit are discarded (a quality filter); '+\
              'pass 0.0 to disable [default:%default]'),
 make_option('--step1_otu_map_fp',type='existing_filepath',
             help='reference OTU picking OTU map '+\
             ' (to avoid rebuilding if one has already been built)'),
 make_option('--step1_failures_fasta_fp',type='existing_filepath',
             help='reference OTU picking failures fasta filepath '+\
             ' (to avoid rebuilding if one has already been built)'),
 make_option('--suppress_step4',action='store_true',default=False,
             help='suppress the final de novo OTU picking step '+\
             ' (may be necessary for extremely large data sets) [default: %default]'),
 make_option('--min_otu_size',type='int',default=2,
             help='the minimum otu size (in number of sequences) to retain the otu '
             '[default: %default]')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    verbose = opts.verbose
    
    input_fps = opts.input_fps
    refseqs_fp = opts.reference_fp
    output_dir = opts.output_dir
    verbose = opts.verbose
    print_only = False
    percent_subsample = opts.percent_subsample
    new_ref_set_id = opts.new_ref_set_id
    prefilter_percent_id = opts.prefilter_percent_id
    if prefilter_percent_id == 0.0:
        prefilter_percent_id = None
    
    parallel = opts.parallel
    # No longer checking that jobs_to_start > 2, but
    # commenting as we may change our minds about this.
    #if parallel: raise_error_on_parallel_unavailable()
    
    if opts.parameter_fp:
        try:
            parameter_f = open(opts.parameter_fp)
        except IOError:
            raise IOError,\
             "Can't open parameters file (%s). Does it exist? Do you have read access?"\
             % opts.parameter_fp
        params = parse_qiime_parameters(parameter_f)
    else:
        params = parse_qiime_parameters([]) 
        # empty list returns empty defaultdict for now
    
    jobs_to_start = opts.jobs_to_start
    default_jobs_to_start = qiime_config['jobs_to_start']
    validate_and_set_jobs_to_start(params,
                                   jobs_to_start,
                                   default_jobs_to_start,
                                   parallel,
                                   option_parser)
    
    try:
        makedirs(output_dir)
    except OSError:
        if opts.force:
            pass
        else:
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

    if len(input_fps) == 1:
        pick_subsampled_open_referenence_otus(input_fp=input_fps[0], 
                                  refseqs_fp=refseqs_fp,
                                  output_dir=output_dir,
                                  percent_subsample=percent_subsample,
                                  new_ref_set_id=new_ref_set_id,
                                  command_handler=command_handler,
                                  params=params,
                                  min_otu_size=opts.min_otu_size,
                                  qiime_config=qiime_config,
                                  prefilter_percent_id=prefilter_percent_id,
                                  step1_otu_map_fp=opts.step1_otu_map_fp,
                                  step1_failures_fasta_fp=opts.step1_failures_fasta_fp,
                                  parallel=parallel,
                                  suppress_step4=opts.suppress_step4,
                                  logger=None,
                                  status_update_callback=status_update_callback)
    else:    
        iterative_pick_subsampled_open_referenence_otus(input_fps=input_fps,
                              refseqs_fp=refseqs_fp,
                              output_dir=output_dir,
                              percent_subsample=percent_subsample,
                              new_ref_set_id=new_ref_set_id,
                              command_handler=command_handler,
                              params=params,
                              min_otu_size=opts.min_otu_size,
                              qiime_config=qiime_config,
                              prefilter_percent_id=prefilter_percent_id,
                              step1_otu_map_fp=opts.step1_otu_map_fp,
                              step1_failures_fasta_fp=opts.step1_failures_fasta_fp,
                              parallel=parallel,
                              suppress_step4=opts.suppress_step4,
                              logger=None,
                              status_update_callback=status_update_callback)



if __name__ == "__main__":
    main()