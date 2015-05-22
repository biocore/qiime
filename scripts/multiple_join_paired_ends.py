#!/usr/bin/env python
from __future__ import division
 
__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
 
from os.path import abspath, join
from os import walk
 
from qiime.util import (parse_command_line_parameters,
                        make_option,
                        create_dir,
                        load_qiime_config)
from qiime.workflow.util import (print_commands,
                                 call_commands_serially,
                                 generate_log_fp,
                                 WorkflowLogger,
                                 no_status_updates,
                                 get_params_str)
from qiime.parse import parse_qiime_parameters
from qiime.workflow.preprocess import get_pairs, create_commands_jpe
 
script_info={}
script_info['brief_description']="""Run join_paired_ends.py on multiple files."""
script_info['script_description']= \
"""This script runs join_paired_ends.py on data that are already demultiplexed
(split up according to sample, with one sample per pair of files). The script
supports the following types of input:

- a directory containing many files, where each file is named on a per-sample
  basis
- a directory containing many directories, where each directory is named on a
  per-sample basis
 
The script assumes that the leading/trailing characters before/after the read
number indicator (see --read1_indicator) are matched between forward and
reverse reads. For example:

- S0_L001_R1_001.fastq.gz and S0_L001_R2_001.fastq.gz would be matched up reads
- S0_L002_R1_00X.fastq.gz and S0_L002_R2_00X.fastq.gz would be matched up reads

If an optional --barcode_indicator file is used, it is searched for in the same
manner that the paired files are searched for, so if the default "_I1_" is
used, S0_L001_R1_001.fastq.gz and S0_L001_R2_001.fastq.gz would be matched up
with S0_L001_I1_001.fastq.gz as the barcode indicator file.

The output directory used for each call to join_paired_ends.py uses the base
name of the input read 1 fastq file (a single directory would be problematic
since the output names for join_paired_ends.py can be the same for different
calls). Use the parameter --include_input_dir_path to also include the input
directory name in the output directory path, which may be preferable in the
case of an input folder of folders, and --remove_filepath_in_name can be used
in this case to prevent the input read 1 fastq file base name from being used
as part of the output directory name.

"""

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("Example 1:",
     "Process an input folder of paired-up files (by filename, with the "
     "default _R1_ and _R2_ containing the forward and reverse reads "
     "filenames, respectively). An optional parameters file is passed with "
     "-p. This file can specify an optional parameter for "
     "join_paired_ends.py, such as: join_paired_ends:pe_join_method SeqPrep",
     "%prog -i input_files -o output_folder -p qiime_parameters.txt"))

script_info['script_usage'].append(
    ("Example 2:",
     "Process an input folder of folders (with the filenames having "
     "_forward_ and _reverse_ containing the forward and reverse read "
     "filenames, respectively). The individual folder names are included in "
     "the output folder names, but not the filenames. A matching barcode "
     "fastq file (indicated by _barcode_) is also included.",
     "%prog -i input_folders -o output_folder --read1_indicator '_forward_' "
     "--read2_indicator '_reverse_' --include_input_dir_path "
     "--remove_filepath_in_name -b --barcode_indicator '_barcode_'"))

script_info['script_usage'].append(
    ("Example 3:",
     "To see what commands would be executed by the script without actually "
     "running them, use the following command:",
     "%prog -i input_files -o output_folder -w"))

script_info['output_description']= ("The output of running "
                                    "join_paired_ends.py on many input files. "
                                    "See script description for more details.")

script_info['required_options']= [
    make_option('-i', '--input_dir', type='existing_dirpath',
        help='Input directory of directories, or directory of paired fastq '
        'files.'),
    make_option('-o', '--output_dir', type='new_dirpath',
        help='Base output directory to write output folders')
]
script_info['optional_options']= [
    make_option('-p', '--parameter_fp', type='existing_filepath',
        help='path to the parameter file, which specifies changes'
        ' to the default behavior of join_paired_ends.py. '
        'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters'
        ' [default: join_paired_ends.py defaults will be used]'),
    make_option('--read1_indicator', default='_R1_',
        help='Substring to search for to indicate read 1 [default: %default]'),
    make_option('--read2_indicator', default='_R2_',
        help='Substring to search for to indicate read 2 [default: %default]'),
    make_option('-b', '--match_barcodes', action='store_true',
        help='Enable searching for matching barcodes '
        '[default: %default]', default=False),
    make_option('--barcode_indicator', default='_I1_',
        help='Substring to search for to indicate barcode reads '
        '[default: %default]'),
    make_option('--leading_text', default='',
        help='Leading text to add to each join_paired_ends.py command '
        '[default: no leading text added]'),
    make_option('--trailing_text', default='',
        help='Trailing text to add to each join_paired_ends.py command '
        '[default: no trailing text added]'),
    make_option('--include_input_dir_path', action='store_true', default=False,
        help='Include the input directory name in the output directory '
        'path. Useful in cases where the file names are repeated in '
        'input folders [default: %default]'),
    make_option('--remove_filepath_in_name', action='store_true', default=False,
        help='Disable inclusion of the input filename in the output '
        'directory names. Must use --include_input_dir_path if this option '
        'is enabled [default: %default]'),
    make_option('-w', '--print_only', action='store_true',
        help='Print the commands but don\'t call them -- '
        'useful for debugging [default: %default]', default=False)
]
        
script_info['version'] = __version__
 
def main():
    option_parser, opts, args = \
        parse_command_line_parameters(suppress_verbose=True, **script_info)
        
    input_dir = opts.input_dir
    parameter_fp = opts.parameter_fp
    read1_indicator = opts.read1_indicator
    read2_indicator = opts.read2_indicator
    match_barcodes = opts.match_barcodes
    barcode_indicator = opts.barcode_indicator
    leading_text = opts.leading_text
    trailing_text = opts.trailing_text
    include_input_dir_path = opts.include_input_dir_path
    output_dir = abspath(opts.output_dir)
    remove_filepath_in_name = opts.remove_filepath_in_name
    print_only = opts.print_only
    
    if remove_filepath_in_name and not include_input_dir_path:
        option_parser.error("If --remove_filepath_in_name is enabled, "
            "--include_input_dir_path must also be enabled.")

    if opts.parameter_fp:
        with open(opts.parameter_fp, 'U') as parameter_f:
            params_dict = parse_qiime_parameters(parameter_f)
        params_str = get_params_str(params_dict['join_paired_ends'])
    else:
        params_dict = {}
        params_str = ""

    create_dir(output_dir)
    
    all_files = []
    extensions = ['.fastq.gz', '.fastq', '.fq.gz', '.fq']
    
    for root, dir, fps in walk(input_dir):
        for fp in fps:
            for extension in extensions:
                if fp.endswith(extension):
                    all_files += [abspath(join(root, fp))]
        
    pairs, bc_pairs = get_pairs(all_files, read1_indicator, read2_indicator,
        match_barcodes, barcode_indicator)

    commands = create_commands_jpe(pairs, output_dir,
        params_str, leading_text, trailing_text, include_input_dir_path,
        remove_filepath_in_name, match_barcodes, bc_pairs)
        
    qiime_config = load_qiime_config()
    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially
    logger = WorkflowLogger(generate_log_fp(output_dir),
                            params=params_dict,
                            qiime_config=qiime_config)
    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback=no_status_updates,
                    logger=logger,
                    close_logger_on_success=True)
 
if __name__ == "__main__":
    main()
