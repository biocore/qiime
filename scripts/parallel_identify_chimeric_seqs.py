#!/usr/bin/env python
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jens Reeder", "Jai Ram Rideout",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.util import make_option
from os.path import split, splitext, join
from qiime.util import (make_option, parse_command_line_parameters,
                        load_qiime_config, get_options_lookup)

from qiime.parallel.identify_chimeric_seqs import ParallelChimericSequenceIdentifier

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Parallel chimera detection"""
script_info[
    'script_description'] = """This script works like the identify_chimeric_seqs.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage'] = []
# copied from identify_chimeric_seqs.py
script_info['script_usage'].append(("""blast_fragments example""", """For each sequence provided as input, the blast_fragments method splits the input sequence into n roughly-equal-sized, non-overlapping fragments, and assigns taxonomy to each fragment against a reference database. The BlastTaxonAssigner (implemented in assign_taxonomy.py) is used for this. The taxonomies of the fragments are compared with one another (at a default depth of 4), and if contradictory assignments are returned the sequence is identified as chimeric. For example, if an input sequence was split into 3 fragments, and the following taxon assignments were returned:

==========  ==========================================================
fragment1:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
fragment2:  Archaea;Euryarchaeota;Halobacteriales;uncultured
fragment3:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
==========  ==========================================================

The sequence would be considered chimeric at a depth of 3 (Methanobacteriales vs. Halobacteriales), but non-chimeric at a depth of 2 (all Euryarchaeota).

blast_fragments begins with the assumption that a sequence is non-chimeric, and looks for evidence to the contrary. This is important when, for example, no taxonomy assignment can be made because no blast result is returned. If a sequence is split into three fragments, and only one returns a blast hit, that sequence would be considered non-chimeric. This is because there is no evidence (i.e., contradictory blast assignments) for the sequence being chimeric. This script can be run by the following command, where the resulting data is written to $PWD/blast_fragments_chimeric_seqs.txt and using default parameters (i.e., number of fragments ("-n 3"), taxonomy depth ("-d 4") and maximum E-value ("-e 1e-30")). ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""", """%prog -i $PWD/inseqs.fasta -t $PWD/id_to_tax.txt -r $PWD/refseqs.fasta -o $PWD/blast_fragments_chimeric_seqs.txt -m blast_fragments"""))

script_info['script_usage'].append(("""ChimeraSlayer Example:""",
                                    """Identify chimeric sequences using the ChimeraSlayer algorithm against a user provided reference database. The input sequences need to be provided in aligned (Py)Nast format and the reference database needs to be provided as aligned FASTA (-a). Note that the reference database needs to be the same that was used to build the alignment of the input sequences! ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""",
                                    """%prog -i $PWD/inseqs_aligned.fasta -o $PWD/chimera_slayer_chimeric_seqs.txt"""))

script_info[
    'output_description'] = """The result of parallel_identify_chimeric_seqs.py is a text file that identifies which sequences are chimeric."""

script_info['required_options'] = [
    options_lookup['fasta_as_primary_input'],
]

chimera_detection_method_choices = ['blast_fragments', 'ChimeraSlayer']

script_info['optional_options'] = [
    make_option('-a', '--aligned_reference_seqs_fp',
                type='existing_filepath',
                default=qiime_config['pynast_template_alignment_fp'],
                help='Path to (Py)Nast aligned reference sequences. '
                'REQUIRED when method ChimeraSlayer [default: %default]'),

    make_option('-t', '--id_to_taxonomy_fp',
                type='existing_filepath',
                help='Path to tab-delimited file mapping sequences to assigned '
                'taxonomy. Each assigned taxonomy is provided as a comma-separated '
                'list. [default: %default; REQUIRED when method is blast_fragments]'),

    make_option('-r', '--reference_seqs_fp',
                type='existing_filepath',
                help='Path to reference sequences (used to build a blast db when method blast_fragments). '
                '[default: %default; REQUIRED when method blast_fragments' +
                ' if no blast_db is provided;]'),

    make_option('-b', '--blast_db', type='blast_db',
                help='Database to blast against. Must provide either --blast_db or '
                '--reference_seqs_fp when method is blast_fragments [default: %default]'),

    make_option('-m', '--chimera_detection_method',
                type='choice', help='Chimera detection method. Choices: ' +
                " or ".join(chimera_detection_method_choices) +
                '. [default:%default]',
                choices=chimera_detection_method_choices, default='ChimeraSlayer'),

    make_option('-n', '--num_fragments',
                type='int', help='Number of fragments to split sequences into' +
                ' (i.e., number of expected breakpoints + 1) [default: %default]',
                default=3),

    make_option('-d', '--taxonomy_depth',
                type='int', help='Number of taxonomic divisions to consider' +
                ' when comparing taxonomy assignments [default: %default]',
                default=4),

    make_option('-e', '--max_e_value',
                type='float', help='Max e-value to assign taxonomy' +
                ' [default: %default]', default=1e-30),

    make_option('--min_div_ratio',
                type='float', help='min divergence ratio ' +
                '(passed to ChimeraSlayer). If set to None uses ' +
                'ChimeraSlayer default value. ' +
                ' [default: %default]', default=None),

    make_option('-o', '--output_fp',
                type='new_filepath',
                help='Path to store output [default: derived from input_seqs_fp]'),

    # Define parallel-script-specific parameters
    options_lookup['jobs_to_start'],
    options_lookup['retain_temp_files'],
    options_lookup['suppress_submit_jobs'],
    options_lookup['poll_directly'],
    options_lookup['cluster_jobs_fp'],
    options_lookup['suppress_polling'],
    options_lookup['job_prefix'],
    options_lookup['seconds_to_sleep']
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create dict of command-line options.
    params = eval(str(opts))

    # Additional option checks (copied from scripts/identify_chimeric_seqs.py).
    if opts.chimera_detection_method == 'blast_fragments':
        if not (opts.blast_db or opts.reference_seqs_fp):
            option_parser.error('Must provide either --blast_db or' +
                                ' --reference_seqs_fp and --id_to_taxonomy_fp when' +
                                ' method is blast_fragments.')
        if not opts.id_to_taxonomy_fp:
            option_parser.error('Must provide --id_to_taxonomy_fp when method' +
                                ' is blast_fragments.')

        if opts.num_fragments < 2:
            option_parser.error('Invalid number of fragments (-n %d) Must be >= 2.'
                                % opts.num_fragments)
    elif opts.chimera_detection_method == 'ChimeraSlayer':
        if not opts.aligned_reference_seqs_fp:
            option_parser.error("Must provide --aligned_reference_seqs_fp " +
                                "when using method ChimeraSlayer")

    # Set the output_fp if not set.
    output_fp = opts.output_fp
    if not output_fp:
        input_basename = splitext(split(opts.input_fasta_fp)[1])[0]
        output_fp = '%s_chimeric.txt' % input_basename
        params['output_fp'] = output_fp

    # Find the output dir path based on the output file path.
    output_dir, _ = split(output_fp)
    if output_dir == "":
        output_dir = "./"

    parallel_runner = ParallelChimericSequenceIdentifier(
        cluster_jobs_fp=opts.cluster_jobs_fp,
        jobs_to_start=opts.jobs_to_start,
        retain_temp_files=opts.retain_temp_files,
        suppress_polling=opts.suppress_polling,
        seconds_to_sleep=opts.seconds_to_sleep)

    parallel_runner(opts.input_fasta_fp,
                    output_dir,
                    params,
                    job_prefix=opts.job_prefix,
                    poll_directly=opts.poll_directly,
                    suppress_submit_jobs=opts.suppress_submit_jobs)


if __name__ == "__main__":
    main()
