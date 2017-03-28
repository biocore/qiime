#!/usr/bin/env python
from __future__ import division
 
__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
 
from glob import glob
from os.path import abspath
 
from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option,
                        create_dir,
                        load_qiime_config)
from qiime.workflow.util import (print_commands,
                                 call_commands_serially,
                                 generate_log_fp,
                                 WorkflowLogger,
                                 no_status_updates)
from qiime.generate_join_paired_ends_commands import (get_pairs,
                                                      create_commands)
 
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""Creates and executes commands for running join_paired_ends.py on multiple files"""
script_info['script_description']="""The purpose of this script is to automate the creation of join_paired_ends.py commands for multiple files. This script is intended to help in cases of having data already demultiplexed (split up according to SampleID), either as a folder of paired files (with different prefixes before the read number) or a folder of folders, each containing a paired read 1 and 2 file.
 
There are assumptions made by this script (necessary to pair up matched files) that the leading/trailing characters before the read number indicator (see --read1_indicator) are matched between forward and reversed reads. E.g.:
S0_L001_R1_001.fastq.gz and S0_L001_R2_001.fastq.gz would be matched up reads and
S0_L002_R1_00X.fastq.gz and S0_L002_R2_00X.fastq.gz would be matched up reads.

The output directory used for each call to join_paired_ends.py uses the base name of the input read 1 fastq file (a single directory would be problematic since the output names for join_paired_ends.py can be the same for different calls). Use the parameter --include_input_dir_path to also include the input directory name in the output directory path, which may be preferable in the case of an input folder of folders, and the --remove_filepath_in_name can be used in this case.

"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example of an input folder of paired-up files (by filename, with the default _R1_ and _R2_ contained the forward and reverse reads filenames, respectively):""",
                                    """The input folder is named 'input_paired_files', the output directory is 'output_folder', and an optional parameter is added of '--pe_join_method SeqPrep'""",
                                    """%prog -i input_paired_files/ -o output_folder/ --type_of_input folder_of_files --optional_params '--pe_join_method SeqPrep'"""))
script_info['script_usage'].append(("""Example of an input folder of folders (with the filenames having _forward_ and _reverse_ contained the forward and reverse read filenames, respectively):""",
                                    """The input folder is named 'input_folders', the output directory is 'output_folder', and an optional parameter is added of '--pe_join_method SeqPrep'. Additionally, the individual folder names are included in the output folder names, but not the filenames.""",
                                    """%prog -i input_folders/ -o output_folder/ --type_of_input folder_of_folders --optional_params '--pe_join_method SeqPrep' --read1_indicator '_forward_' --read2_indicator '_reverse_' --include_input_dir_path --remove_filepath_in_name"""))

script_info['output_description']="""The commands will be executed by the script by default. The commands can be printed with the -w option, so one could write these directly to an output file by calling generate_join_paired_ends_commands.py X > Y where X are the input parameters, and Y is an output text file. A log file will also be generated."""
script_info['required_options']= [\
    make_option('-i', '--input_dir',type='existing_dirpath',
        help='Input directory of directories, or directory of paired fastq '+\
        'files.'),
    make_option('-o', '--base_output_dir',type='new_dirpath',
        help='Base output directory to write output folders')
]
script_info['optional_options']= [\
    make_option('--type_of_input', type='choice', default='folder_of_files',
        help='Selects whether a folder of files (folder_of_files) or a folder '
        'of folders (folder_of_folders) will '
        'will be searched for by this script. To avoid errors or problems, '
        'only paired read files should be used with this script. '
        ' [default: %default]', choices=['folder_of_files',
        'folder_of_folders']),
    make_option('--optional_params', type='str', default='',
        help='Added optional parameter string to use in join_paired_ends.py, '
        'such as "-m SeqPrep" [default: %default]'),
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
    type_of_input = opts.type_of_input
    optional_params = opts.optional_params
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
    if type_of_input == "folder_of_files" and remove_filepath_in_name:
        option_parser.error("For input of type folder_of_files, the option "
            "--remove_filepath_in_name can not be enabled, as this would "
            "preclude creation of unique output folders.")
    
    create_dir(base_output_dir)
                
    if type_of_input == "folder_of_files":
        all_files = glob(abspath(input_dir) + "/*.fastq*")
    else:
        all_files = []
        all_files_folders = glob(abspath(input_dir) + "/*")
        for curr_folder in all_files_folders:
            all_files += (glob(abspath(curr_folder) + "/*.fastq*"))
        
    pairs = get_pairs(all_files, read1_indicator, read2_indicator)
    
    commands = create_commands(pairs, base_output_dir,
        optional_params, leading_text, trailing_text, include_input_dir_path,
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
                    close_logger_on_success=close_logger_on_success)
 
if __name__ == "__main__":
    main()