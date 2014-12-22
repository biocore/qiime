#!/usr/bin/env python
from __future__ import division
 
__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
 
from os.path import abspath, join
from os import walk
 
from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
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
from qiime.workflow.preprocess import (get_pairs,
                                       create_commands_eb)
 
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""Creates and executes commands for running extract_barcodes.py on multiple files"""
script_info['script_description']="""The purpose of this script is to automate the creation of extract_barcodes.py commands for multiple files. This script is intended to help in cases of having data already demultiplexed (split up according to sample), either as a folder of files (with different prefixes before the read number) or a folder of folders.
 
There are assumptions made by this script (necessary to pair up matched files) that the leading/trailing characters before the read number indicator (see --read1_indicator) are matched between forward and reversed reads. E.g.:
S0_L001_R1_001.fastq.gz and S0_L001_R2_001.fastq.gz would be matched up reads and
S0_L002_R1_00X.fastq.gz and S0_L002_R2_00X.fastq.gz would be matched up reads.

The output directory used for each call to extract_barcodes.py uses the base name of the input read 1 fastq file (a single directory would be problematic since the output names for extract_barcodes.py can be the same for different calls). Use the parameter --include_input_dir_path to also include the input directory name in the output directory path, which may be preferable in the case of an input folder of folders, and the --remove_filepath_in_name can be used in this case.

"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example of an input folder of files, with default options used for extract_barcodes:""",
                                    """The input folder is named 'input_paired_files', the output directory is 'output_folder'.""",
                                    """%prog -i input_stitched_files/ -o output_folder/"""))
script_info['script_usage'].append(("""Example of an input folder of folders (with the filenames having _forward_ and _reverse_ contained the forward and reverse read filenames, respectively), using the extract_barcodes option for paired reads:""",
                                    """The input folder is named 'input_folders', the output directory is 'output_folder'. Additionally, the individual folder names are included in the output folder names, but not the filenames. Note-it is important to pass the --paired_data parameter if paired data are to be used in the extract_barcodes.py commands. Additionally, the paired fastq file type for extract_barcodes is specified with a qiime_parameters.txt file (with this value specified: extract_barcodes:input_type barcode_paired_end)""",
                                    """%prog -i input_folders/ -o output_folder/ -p qiime_parameters.txt --paired_data --read1_indicator '_forward_' --read2_indicator '_reverse_' --include_input_dir_path --remove_filepath_in_name"""))

script_info['output_description']="""The commands will be executed by the script by default. The commands can be printed with the -w option, so one could write these directly to an output file by calling multiple_extract_barcodes.py X > Y where X are the input parameters, and Y is an output text file. A log file will also be generated."""
script_info['required_options']= [\
    make_option('-i', '--input_dir',type='existing_dirpath',
        help='Input directory of directories, or directory of paired fastq '+\
        'files.'),
    make_option('-o', '--base_output_dir',type='new_dirpath',
        help='Base output directory to write output folders')
]
script_info['optional_options']= [\
    make_option('-p', '--parameter_fp', type='existing_filepath',
        help='path to the parameter file, which specifies changes' +
        ' to the default behavior. ' +
        'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .' +
        ' [if omitted, default values will be used]'),
    make_option('--paired_data', action='store_true',
        help='Turn this option on if paired data are to be used. The type of '
        'paired data for extract_barcodes.py should be specified with the '
        '--parameter_fp parameter. Forward and reverse reads will be '
        'searched for via the --read1 and --read2_indicator parameters. '
        ' [default: %default]', default=False),
    make_option('--read1_indicator', type='str', default='_R1_',
        help='Substring to search for to indicate read 1. [default: %default]'),
    make_option('--read2_indicator', type='str', default='_R2_',
        help='Substring to search for to indicate read 2. [default: %default]'),
    make_option('--leading_text', type='str', default='',
        help='Added leading text to command [default: %default]'),
    make_option('--trailing_text', type='str', default='',
        help='Added trailing text to command [default: %default]'),
    make_option('--include_input_dir_path', action='store_true', default=False,
        help='Include the input directory name in the output directory '
        'path. Useful in cases where the file names are repeated in '
        'input folders [default: %default]'),
    make_option('--remove_filepath_in_name', action='store_true', default=False,
        help='Disable inclusion of the input filename in the output '
        'directory names. Must use --include_input_dir_path if this option '
        'is enabled. [default: %default]'),
    make_option('-w', '--print_only', action='store_true',
        dest='print_only', help='Print the commands but don\'t call them -- ' +
        'useful for debugging [default: %default]', default=False),
 
]
        
script_info['version'] = __version__
 
def main():
    option_parser, opts, args =\
        parse_command_line_parameters(suppress_verbose=True, **script_info)
        
    input_dir = opts.input_dir
    paired_data = opts.paired_data
    parameter_fp = opts.parameter_fp
    read1_indicator = opts.read1_indicator
    read2_indicator = opts.read2_indicator
    leading_text = opts.leading_text
    trailing_text = opts.trailing_text
    include_input_dir_path = opts.include_input_dir_path
    base_output_dir = abspath(opts.base_output_dir)
    remove_filepath_in_name = opts.remove_filepath_in_name
    print_only = opts.print_only
    
    if remove_filepath_in_name and not include_input_dir_path:
        option_parser.error("If --remove_filepath_in_name enabled, "
            "--include_input_dir_path must be enabled.")
            
    if opts.parameter_fp:
        try:
            parameter_f = open(opts.parameter_fp, 'U')
        except IOError:
            raise IOError("Can't open parameters file (%s). Does it exist? Do you have read access?"
                % opts.parameter_fp)
        params = get_params_str(parse_qiime_parameters(parameter_f)['extract_barcodes'])
        parameter_f.close()
    else:
        params = ""
    
    create_dir(base_output_dir)
                
    all_files = []
    extensions = ['.fastq.gz', '.fastq', '.fq.gz', '.fq']
    
    for root, dir, fps in walk(input_dir):
        for fp in fps:
            for extension in extensions:
                if fp.endswith(extension):
                    all_files += [abspath(join(root, fp))]
        
    if paired_data:
        all_files, bc_pairs =\
            get_pairs(all_files, read1_indicator, read2_indicator)
    
    commands = create_commands_eb(all_files, paired_data, base_output_dir,
        params, leading_text, trailing_text, include_input_dir_path,
        remove_filepath_in_name)
        
    qiime_config = load_qiime_config()
    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially
    logger = WorkflowLogger(generate_log_fp(base_output_dir),
                                params={},
                                qiime_config=qiime_config)
    close_logger_on_success = True
    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback = no_status_updates,
                    logger=logger,
                    close_logger_on_success = True)
 
if __name__ == "__main__":
    main()