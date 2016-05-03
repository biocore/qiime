#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from os.path import abspath, join
from os import walk
from fnmatch import filter

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
from qiime.workflow.preprocess import (get_matching_files,
                                       create_commands_slf)

script_info={}
script_info['brief_description']="""Run split_libraries_fastq.py on multiple files."""
script_info['script_description']="""
This script runs split_libraries_fastq.py on data that are already demultiplexed
(split up according to sample, with one sample per file). The script supports
the following types of input:

- a directory containing many files, where each file is named on a per-sample
  basis (with different prefixes before the read number)
- a directory containing many directories, where each directory is named on a
  per-sample basis

This script assumes that the leading characters before the read indicator
(see --read_indicator) are matched between the read, barcode, and mapping files.
For example, sample1_L001_R1_001.fastq.gz, sample1_L001_I1_001.fastq.gz,
sample1_L001_mapping_001.txt would be matched up if "R1" is the read indicator,
"I1" is the barcode indicator, and "mapping" is the mapping file indicator.

Depending on the inputed arguments, this script can create any of the following
files in the output directory:

    * seqs.txt: this file will contain a list (one element per line) of found
      sequences filepaths.
    * barcodes.txt: this file will contain a list (one element per line) of
      found barcodes filepaths.
    * map.txt: this file wil contain a list (one element per line) of mapping
      file filepaths.
    * sample_ids.txt: this file will contain a list (one element per line)
      sample identifiers.
"""
script_info['script_usage'] = []

script_info['script_usage'].append(
    ("Example 1:",
     "Process an input folder of folders, with options specified to pair up "
     "reads, barcodes, and mapping files. A qiime_parameters.txt file is "
     "included to use the parameter split_libraries_fastq:barcode_type 12",
     "%prog -i input_folders -o output_folder --demultiplexing_method "
     "mapping_barcode_files --read_indicator reads --barcode_indicator barcode "
     "--mapping_indicator mapping -p qiime_parameters.txt"))

script_info['script_usage'].append(
    ("Example 2:",
     "Process an input folder of files, with the option specified to generate "
     "sample ids using the filenames containing the text _R1_ (behavior is to use "
     "all text before the first underscore as the sample id)",
     "%prog -i input_files -o output_folder --read_indicator '*_R1_*'"
     " --demultiplexing_method sampleid_by_file"))

script_info['script_usage'].append(
    ("Example 3:",
     "Process an input folder of folders, with an option specified to "
     "use the folder names as the sample ids. In this case, the fastq "
     "filenames themselves are not included, only the folder names are used. "
     "The target reads in this case have the text _reads_ in the filenames. ",
     "%prog -i input_folders_no_barcodes -o output_folder "
     "--demultiplexing_method sampleid_by_file --include_input_dir_path "
     "--remove_filepath_in_name --read_indicator '*_reads_*'"))

script_info['script_usage'].append(
 ("Example 4:",
  "To see what commands would be executed by the script without actually "
  "running them, use the following command (target reads include _R1_ in the filename):",
  "%prog -i input_files -o output_folder -w --read_indicator '*_R1_*'"))

script_info['output_description']= (
    "The output of running split_libraries_fastq.py on many input files. "
    "See script description for more details.")

script_info['required_options']= [
    make_option('-i', '--input_dir',type='existing_dirpath',
        help='Input directory of directories or fastq files.'),
    make_option('-o', '--output_dir',type='new_dirpath',
        help='Output directory to write split_libraries_fastq.py results')
]
script_info['optional_options']= [
    make_option('-m', '--demultiplexing_method', type='choice',
        choices = ['sampleid_by_file', 'mapping_barcode_files'],
        default = "sampleid_by_file", help=('Method for demultiplexing. Can '
        'either be "sampleid_by_file" or "mapping_barcode_files". With the '
        'sampleid_by_file option, each fastq file (and/or directory name) '
        'will be used to generate the --sample_ids value passed to '
        'split_libraries_fastq.py. The mapping_barcode_files option will '
        'search for barcodes and mapping files that match the input read '
        'files [default: %default]')),
    make_option('-p', '--parameter_fp', type='existing_filepath',
        help='path to the parameter file, which specifies changes'
        ' to the default behavior of split_libraries_fastq.py. '
        'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters'
        ' [default: split_libraries_fastq.py defaults will be used]'),
    make_option('--read_indicator', default='*_R1_*',
        help='Substring to search for to indicate read files, when '
        '--demultiplexing_method is sampleid_by_file, wildcards can be used, e.g. \'*\' '
        'for all files. If multiple fastq files are present as in the case after joining '
        'reads, then one can select for joined reads with a value such as '
        '\'*fastqjoin.join*\''
        ' [default: %default]'),
    make_option('--barcode_indicator', default='_I1_',
        help='Substring to search for to indicate barcode files'
        ' [default: %default]'),
    make_option('--mapping_indicator', default='_mapping_',
        help='Substring to search for to indicate mapping files'
        ' [default: %default]'),
    make_option('--mapping_extensions', default='txt,tsv',
        help='Comma-separated list of file extensions used to identify mapping '
        'files. Only applies when --demultiplexing_method is '
        '"mapping_barcode_files" [default: %default]'),
    make_option('--sampleid_indicator', default='_',
        help='Text in fastq filename before this value will be used '
        'as output sample ids [default: %default]'),
    make_option('--include_input_dir_path', action='store_true', default=False,
        help='Include the input directory name in the output sample id '
        'name. Useful in cases where the file names are repeated in '
        'input folders [default: %default]'),
    make_option('--remove_filepath_in_name', action='store_true', default=False,
        help='Disable inclusion of the input filename in the output '
        'sample id names. Must use --include_input_dir_path if this option '
        'is enabled [default: %default]'),
    make_option('--leading_text', default='',
        help='Leading text to add to each split_libraries_fastq.py command '
        '[default: no leading text added]'),
    make_option('--trailing_text', default='',
        help='Trailing text to add to each split_libraries_fastq.py command '
        '[default: no trailing text added]'),
    make_option('-w', '--print_only', action='store_true',
        help='Print the commands but don\'t call them -- ' +
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
    mapping_extensions = opts.mapping_extensions.split(',')
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
        with open(opts.parameter_fp, 'U') as parameter_f:
            params_dict = parse_qiime_parameters(parameter_f)
        params_str = get_params_str(params_dict['split_libraries_fastq'])
    else:
        params_dict = {}
        params_str = ""

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
                for mapping_extension in mapping_extensions:
                    if fp.endswith(mapping_extension):
                        all_mapping += [abspath(join(root, fp))]

        all_files = get_matching_files(all_fastq, all_mapping,
            read_indicator, barcode_indicator, mapping_indicator)
    else:
        # Filter down files to only the target files, raise error if nothing found
        all_files = filter(all_fastq, read_indicator)
        if not all_files:
            raise ValueError,("No reads detected-please check the values indicated with "
                "the --read_indicator parameter. Set as '*' to include all files, or use "
                "a value such as '*fastqjoin.join*' to detect only the reads that are "
                "joined after join_paired_ends.py.")

    commands = create_commands_slf(all_files, demultiplexing_method, output_dir,
        params_str, leading_text, trailing_text, include_input_dir_path,
        remove_filepath_in_name, sampleid_indicator)

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
