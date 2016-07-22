#!/usr/bin/env python
# File created on 20 Jun 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Jose A. Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from os import makedirs
from qiime.util import load_qiime_config
from qiime.parse import parse_qiime_parameters
from qiime.workflow.util import (print_commands,
                                 call_commands_serially,
                                 print_to_stdout,
                                 no_status_updates,
                                 validate_and_set_jobs_to_start)
from qiime.workflow.ampliconnoise import run_ampliconnoise
import os


qiime_config = load_qiime_config()

# summarize_taxa_through_plots.py
options_lookup = get_options_lookup()
script_info = {}
script_info['brief_description'] = "Run AmpliconNoise"
script_info['script_description'] = """
The steps performed by this script are:

1. Split input sff.txt file into one file per sample

2. Run scripts required for PyroNoise

3. Run scripts required for SeqNoise

4. Run scripts requred for Perseus (chimera removal)

5. Merge output files into one file similar to the output of split_libraries.py

This script produces a denoised fasta sequence file such as:
>PC.355_41
CATGCTGCCTC...
...
>PC.636_23
CATGCTGCCTC...
...

Additionally, the intermediate results of the ampliconnoise pipeline are
written to an output directory.

Ampliconnoise must be installed and correctly configured, and parallelized
steps will be called with mpirun, not qiime's start_parallel_jobs_torque.py script.

"""
script_info['script_usage'] = [(
    "",
    "Run ampliconnoise, write output to anoise_out.fna, compatible with output of split_libraries.py",
    "%prog -i Fasting_Example.sff.txt -m Fasting_Map.txt -o anoise_out.fna")]
script_info[
    'output_description'] = "a fasta file of sequences, with labels as:'>sample1_0' , '>sample1_1' ..."
script_info['required_options'] = [
    options_lookup['mapping_fp'],
    make_option(
        '-i',
        '--sff_filepath',
        type='existing_filepath',
        help='sff.txt filepath'),
    make_option(
        '-o',
        '--output_filepath',
        type='new_filepath',
        help='the output file'),


]
script_info['optional_options'] = [
    make_option(
        '-n',
        '--np',
        type='int',
        default=2,
        help='number of processes to use for mpi steps. Default: %default'),
    make_option(
        '--chimera_alpha',
        type='float',
        default=-
        3.8228,
        help='alpha value to Class.pl used for chimera removal  Default: %default'),
    make_option(
        '--chimera_beta',
        type='float',
        default=0.6200,
        help='beta value to Class.pl used for chimera removal  Default: %default'),
    make_option(
        '--seqnoise_resolution',
        type='string',
        default=None,
        help='-s parameter passed to seqnoise. Default is 25.0 for titanium, 30.0 for flx'),
    make_option(
        '-d',
        '--output_dir',
        type='new_dirpath',
        default=None,
        help='directory for ampliconnoise intermediate results. Default is output_filepath_dir'),
    make_option('-p', '--parameter_fp', type='existing_filepath',
                help='path to the parameter file, which specifies changes' +
                ' to the default behavior. ' +
                'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters.' +
                ' [if omitted, default values will be used]'),
    make_option('-f', '--force', action='store_true', default=False,
                dest='force', help='Force overwrite of existing output directory' +
                ' (note: existing files in output_dir will not be removed)' +
                ' [default: %default]'),
    make_option('-w', '--print_only', action='store_true', default=False,
                dest='print_only', help='Print the commands but don\'t call them -- ' +
                'useful for debugging [default: %default]'),
    make_option(
        '--suppress_perseus',
        action='store_true',
        default=False,
        help='omit perseus from ampliconnoise workflow'),
    make_option(
        '--platform',
        type='choice',
        choices=['titanium',
                 'flx'],
        default='flx',
        help="sequencing technology, options are 'titanium','flx'. [default: %default]"),
    make_option(
        '--truncate_len',
        type='int',
        default=None,
        help="Specify a truncation length for ampliconnoise.  Note that is this is not specified, the truncate length is chosen by the --platform option (220 for FLX, 400 for Titanium) [default: %default]")


]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    if opts.output_dir is None:
        opts.output_dir = opts.output_filepath + '_dir'

    if opts.parameter_fp:
        try:
            parameter_f = open(opts.parameter_fp, 'U')
        except IOError:
            raise IOError("Can't open parameters file (%s). Does it exist? Do you have read access?"
                          % opts.parameter_fp)
        params = parse_qiime_parameters(parameter_f)
        parameter_f.close()
    else:
        params = parse_qiime_parameters([])
        # empty list returns empty defaultdict for now

    try:
        makedirs(opts.output_dir)
    except OSError:
        if opts.force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            option_parser.error("Output directory already exists. Please choose"
                                " a different directory, or force overwrite with -f.")

    if opts.print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially

    if opts.verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates

    # set env variable
    if opts.platform == 'flx':
        existing_pyro_fp = os.environ['PYRO_LOOKUP_FILE']
        new_pyro_fp = os.path.join(
            os.path.split(existing_pyro_fp)[0],
            'LookUp_E123.dat')
        os.environ['PYRO_LOOKUP_FILE'] = new_pyro_fp
    elif opts.platform == 'titanium':
        existing_pyro_fp = os.environ['PYRO_LOOKUP_FILE']
        new_pyro_fp = os.path.join(
            os.path.split(existing_pyro_fp)[0],
            'LookUp_Titanium.dat')
        os.environ['PYRO_LOOKUP_FILE'] = new_pyro_fp
    else:
        raise RuntimeError(
            'could not find PYRO_LOOKUP_FILE for platform ' +
            platform)

    if opts.truncate_len:
        try:
            truncate_len_int_check = int(opts.truncate_len)
            truncate_len = str(truncate_len_int_check)
        except ValueError:
            raise ValueError(("If specified, truncate_len must be int type."))
    else:
        truncate_len = None

    run_ampliconnoise(
        mapping_fp=opts.mapping_fp,
        output_dir=os.path.abspath(opts.output_dir),
        command_handler=command_handler,
        params=params,
        qiime_config=qiime_config,
        status_update_callback=status_update_callback,
        chimera_alpha=opts.chimera_alpha,
        chimera_beta=opts.chimera_beta,
        sff_txt_fp=opts.sff_filepath,
        numnodes=opts.np,
        suppress_perseus=opts.suppress_perseus,
        output_filepath=os.path.abspath(opts.output_filepath),
        platform=opts.platform,
        seqnoise_resolution=opts.seqnoise_resolution,
        truncate_len=truncate_len
    )

if __name__ == "__main__":
    main()
