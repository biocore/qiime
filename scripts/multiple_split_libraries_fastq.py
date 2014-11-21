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
from qiime.workflow.preprocess import (get_matching_files,
                                       create_commands_slf)
 
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""Creates and executes a command for running split_libraries_fastq.py on multiple input files"""
script_info['script_description']="""The purpose of this script is to automate the creation of a split_libraries_fastq.py command for multiple files. This script is intended to help in cases of having data already demultiplexed (split up according to sample), either as a folder of files (with different prefixes before the read number) or a folder of folders.
 
There are assumptions made by this script (necessary to pair reads and barcodes/mapping files) that the leading characters before the read indicator (see --read_indicator) are matched between the read and barcode/mapping files. E.g.:
sample1_L001_R1_001.fastq.gz, sample1_L001_I1_001.fastq.gz, sample1_L001_mapping_001.txt would be matched up reads, barcode, and mapping files if "R1" is the read indicator, "I1" is the barcode indicator, and "mapping" is the mapping file indicator.

Because the file extension of the mapping file will be different than those of the input fastq files (the script looks for .txt files), the trailing text will be split on the first ".f" that is encountered to match up text following the mapping file indicator and the reads/barcodes files.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example of an input folder of folders, with options specified to pair up reads, barcodes, and mapping files, with a qiime_parameters.txt file included to use the parameter split_libraries_fastq:barcodetype 12""",
                                    """""",
                                    """%prog -i input_folders/ -o output_folder/ --demultiplexing_method 'mapping_barcode_files' --read_indicator 'reads' --barcode_indicator 'barcode' --mapping_indicator 'mapping' -p qiime_parameters.txt"""))
script_info['script_usage'].append(("""Example of an input folder of files, with the option specified to generate SampleID names based upon the filenames (default behavior is to split on the first underscore and name the SampleID by this string).""",
                                    """""",
                                    """%prog -i input_stitched_files/ -o output_folder/ --demultiplexing_method 'sampleid_by_file'"""))
script_info['script_usage'].append(("""Example of an input folder of folders, with the option specified to generate SampleID names based upon the folders containing the fastq files""",
                                    """In this case, the fastq filenames themselves are not included, only the directory names are used""",
                                    """%prog -i input_folders_no_barcodes/ -o output_folder/ --demultiplexing_method 'sampleid_by_file' --include_input_dir_path --remove_filepath_in_name"""))

script_info['output_description']="""The commands will be executed by the script by default. The commands can be printed with the -w option, so one could write these directly to an output file by calling multiple_split_libraries.py X > Y where X are the input parameters, and Y is an output text file. A log file will also be generated."""
script_info['required_options']= [\
    make_option('-i', '--input_dir',type='existing_dirpath',
        help='Input directory of directories or fastq files.'),
    make_option('-o', '--output_dir',type='new_dirpath',
        help='Output directory to write split_libraries_fastq results')
]
script_info['optional_options']= [\
    make_option('-m', '--demultiplexing_method', type='choice',
        choices = ['sampleid_by_file', 'mapping_barcode_files'],
        default = "sampleid_by_file", help=('Method for demultiplexing. Can '
        'either be "sampleid_by_file" or "mapping_barcode_files". With the '
        'sampleid_by_file option, each fastq file (and/or the directory name) '
        'will be used to generate the --sample_ids value passed to '
        'split_libraries_fastq. The mapping_barcode_files option will search '
        'for barcodes and mapping files that match the input read files. '
        ' [default: %default]')),
    make_option('-p', '--parameter_fp', type='existing_filepath',
        help='path to the parameter file, which specifies changes' +
        ' to the default behavior. ' +
        'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .' +
        ' [if omitted, default values will be used]'),
    make_option('--read_indicator', type='str', default='_R1_',
        help='Substring to search for to indicate read files. '
        ' [default: %default]'),
    make_option('--barcode_indicator', type='str', default='_I1_',
        help='Substring to search for to indicate barcodes files. '
        ' [default: %default]'),
    make_option('--mapping_indicator', type='str', default='_mapping_',
        help='Substring to search for to indicate mapping files. '
        ' [default: %default]'),
    make_option('--sampleid_indicator', type='str', default='_',
        help='String in fastq filename before this character will be used '
        'in output SampleID name. [default: %default]'),
    make_option('--include_input_dir_path', action='store_true', default=False,
        help='Include the input directory name in the output SampleID '
        'name. Useful in cases where the file names are repeated in '
        'input folders [default: %default]'),
    make_option('--remove_filepath_in_name', action='store_true', default=False,
        help='Disable inclusion of the input filename in the output '
        'SampleID names. Must use --include_input_dir_path if this option '
        'is enabled. [default: %default]'),
    make_option('--leading_text', type='str', default='',
        help='Added leading text to command [default: %default]'),
    make_option('--trailing_text', type='str', default='',
        help='Added trailing text to command [default: %default]'),
    make_option('-w', '--print_only', action='store_true',
        dest='print_only', help='Print the commands but don\'t call them -- ' +
        'useful for debugging [default: %default]', default=False),
 
]
        
script_info['version'] = __version__
 
def main():
    option_parser, opts, args =\
        parse_command_line_parameters(suppress_verbose=True, **script_info)
        
    input_dir = opts.input_dir
    demultiplexing_method = opts.demultiplexing_method
    parameter_fp = opts.parameter_fp
    read_indicator = opts.read_indicator
    barcode_indicator = opts.barcode_indicator
    mapping_indicator = opts.mapping_indicator
    sampleid_indicator = opts.sampleid_indicator
    leading_text = opts.leading_text
    trailing_text = opts.trailing_text
    include_input_dir_path = opts.include_input_dir_path
    output_dir = abspath(opts.output_dir)
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
        params = get_params_str(parse_qiime_parameters(parameter_f)['split_libraries_fastq'])
        parameter_f.close()
    else:
        params = ""
    
    create_dir(output_dir)
                
    all_fastq = []
    all_mapping = []
    
    extensions = ['.fastq.gz', '.fastq', '.fq.gz', '.fq']
    
    for root, dir, fps in walk(input_dir):
        for fp in fps:
            for extension in extensions:
                if fp.endswith(extension):
                    all_fastq += [abspath(join(root, fp))]
                    
    if demultiplexing_method == 'mapping_barcode_files':
        for root, dir, fps in walk(input_dir):
            for fp in fps:
                if fp.endswith(".txt"):
                    all_mapping += [abspath(join(root, fp))]

        all_files = get_matching_files(all_fastq, all_mapping,
            read_indicator, barcode_indicator, mapping_indicator)
    else:
        all_files = all_fastq
    
    commands = create_commands_slf(all_files, demultiplexing_method, output_dir,
        params, leading_text, trailing_text, include_input_dir_path,
        remove_filepath_in_name, sampleid_indicator)
        
    qiime_config = load_qiime_config()
    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially
    logger = WorkflowLogger(generate_log_fp(output_dir),
                                params={},
                                qiime_config=qiime_config)
    close_logger_on_success = True
    # Call the command handler on the list of commands
    command_handler(commands,
                    status_update_callback = no_status_updates,
                    logger = logger,
                    close_logger_on_success = True)
 
if __name__ == "__main__":
    main()