#!/usr/bin/env python
# File created on 12 Jan 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from shutil import copyfile
from os import makedirs
from os.path import basename, join

from qiime_default_reference import (get_reference_sequences,
                                      get_reference_tree)

from qiime.util import (load_qiime_config, parse_command_line_parameters,
    get_options_lookup, make_option)
from qiime.parse import parse_qiime_parameters
from qiime.workflow.upstream import run_pick_closed_reference_otus
from qiime.workflow.util import (print_commands, call_commands_serially,
    print_to_stdout, no_status_updates, validate_and_set_jobs_to_start)

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

if get_reference_sequences() == qiime_config['pick_otus_reference_seqs_fp']:
    reference_fp_help = (
             "The reference sequences [default: %default]. " +
             "NOTE: If you do not pass -r to this script, you will be using "
             "QIIME's default reference sequences. In this case, QIIME will "
             "copy the corresponding reference tree to the output directory. "
             "This is the tree that should be used to perform phylogenetic "
             "diversity analyses (e.g., with core_diversity_analyses.py).")
else:
    reference_fp_help = "The reference sequences [default: %default]."

script_info = {}
script_info[
    'brief_description'] = "Closed-reference OTU picking/Shotgun UniFrac workflow."
script_info['script_description'] = """
This script picks OTUs using a closed reference and constructs an OTU table.
Taxonomy is assigned using a pre-defined taxonomy map of reference sequence OTU
to taxonomy. If full-length genomes are provided as the reference sequences,
this script applies the Shotgun UniFrac method.

**Note:** If most or all of your sequences are failing to hit the reference,
your sequences may be in the reverse orientation with respect to your reference
database. To address this, you should add the following line to your parameters
file (creating one, if necessary) and pass this file as -p:

pick_otus:enable_rev_strand_match True

Be aware that this doubles the amount of memory used.

"""

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("",
     "Pick OTUs, assign taxonomy, and create an OTU table against a reference set of OTUs. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).",
     "%prog -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus_w_tax/ -t $PWD/taxa.txt"))

script_info['script_usage'].append(
    ("",
     "Pick OTUs and create an OTU table against a reference set of OTUs without adding taxonomy assignments. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).",
     "%prog -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus/"))

script_info['script_usage'].append(
    ("",
     "Pick OTUs, assign taxonomy, and create an OTU table against a reference set of OTUs using usearch_ref. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).",
     "%prog -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus_usearch/ -p $PWD/usearch_params.txt -t $PWD/taxa.txt"))

script_info['script_usage'].append(
    ("",
     "Pick OTUs using usearch_ref, assign taxonomy, and create an OTU table "
     "against a reference set of OTUs using usearch_ref. ALWAYS SPECIFY ABSOLUTE "
     "FILE PATHS (absolute path represented here as $PWD, but will generally look "
     "something like /home/ubuntu/my_analysis/).",
     "%prog -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus_usearch_ref/ "
     "-p $PWD/usearch5.2_params.txt -t $PWD/taxa.txt"))

script_info['script_usage'].append(
    ("",
     "Pick OTUs, assign taxonomy, and create an OTU table against a "
     "reference set of OTUs using sortmerna. ALWAYS SPECIFY ABSOLUTE "
     "FILE PATHS (absolute path represented here as $PWD, but will "
     "generally look something like /home/ubuntu/my_analysis/). ",
     "%prog -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus_sortmerna/ "
     "-p $PWD/sortmerna_params.txt -t $PWD/taxa.txt"))

script_info[
    'script_usage_output_to_remove'] = [
    '$PWD/otus/',
    '$PWD/otus_w_tax/',
    '$PWD/otus_usearch/',
    '$PWD/otus_sortmerna/',
    '$PWD/otus_usearch_ref/']

script_info['output_description'] = ""

script_info['required_options'] = [
    make_option(
        '-i',
        '--input_fp',
        type='existing_filepath',
        help='the input sequences'),
    make_option(
        '-o',
        '--output_dir',
        type='new_dirpath',
        help='the output directory'),
]
script_info['optional_options'] = [
    make_option('-r', '--reference_fp', type='existing_filepath',
                help=reference_fp_help,
                default=qiime_config['pick_otus_reference_seqs_fp']),
    make_option('-p', '--parameter_fp', type='existing_filepath',
                help='path to the parameter file, which specifies changes' +
                ' to the default behavior. ' +
                'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .' +
                ' [if omitted, default values will be used]'),
    make_option('-t', '--taxonomy_fp', type='existing_filepath',
        help='the taxonomy map [default: %default]',
        default=qiime_config['assign_taxonomy_id_to_taxonomy_fp']),
    make_option('-s', '--assign_taxonomy', action='store_true',
                default=False,
                help='Assign taxonomy to each sequence using '
                'assign_taxonomy.py (this will override --taxonomy_fp, if provided) '
                '[default: %default]'),
    make_option('-f', '--force', action='store_true',
                dest='force', help='Force overwrite of existing output directory' +
                ' (note: existing files in output_dir will not be removed)' +
                ' [default: %default]'),
    make_option('-w', '--print_only', action='store_true',
                dest='print_only', help='Print the commands but don\'t call them -- ' +
                'useful for debugging [default: %default]', default=False),
    make_option('-a', '--parallel', action='store_true',
                dest='parallel', default=False,
                help='Run in parallel where available [default: %default]'),
    options_lookup['jobs_to_start_workflow'],
    make_option('--suppress_taxonomy_assignment', action='store_true',
                default=False, help='skip the taxonomy assignment step, resulting in '
                'an OTU table without taxonomy (this will override --taxonomy_fp '
                'and --assign_taxonomy, if provided) [default: %default]'),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose

    input_fp = opts.input_fp
    reference_fp = opts.reference_fp
    taxonomy_fp = opts.taxonomy_fp
    output_dir = opts.output_dir
    verbose = opts.verbose
    print_only = opts.print_only
    assign_taxonomy = opts.assign_taxonomy

    if opts.suppress_taxonomy_assignment:
        assign_taxonomy = False
        taxonomy_fp = None

    parallel = opts.parallel
    # No longer checking that jobs_to_start > 2, but
    # commenting as we may change our minds about this.
    # if parallel: raise_error_on_parallel_unavailable()

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

    jobs_to_start = opts.jobs_to_start
    default_jobs_to_start = qiime_config['jobs_to_start']
    validate_and_set_jobs_to_start(params,
                                   jobs_to_start,
                                   default_jobs_to_start,
                                   parallel,
                                   option_parser)

    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially
        try:
            makedirs(output_dir)
        except OSError:
            if opts.force:
                pass
            else:
                option_parser.error("Output directory already exists. Please choose"
                                    " a different directory, or force overwrite with -f.")

    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates

    run_pick_closed_reference_otus(
        input_fp,
        reference_fp,
        output_dir,
        taxonomy_fp,
        assign_taxonomy=assign_taxonomy,
        command_handler=command_handler,
        params=params,
        qiime_config=qiime_config,
        parallel=parallel,
        status_update_callback=status_update_callback)

    if get_reference_sequences() == reference_fp:
        reference_tree_fp = get_reference_tree()
        fn = basename(reference_tree_fp)
        copyfile(reference_tree_fp, join(output_dir, fn))


if __name__ == "__main__":
    main()
