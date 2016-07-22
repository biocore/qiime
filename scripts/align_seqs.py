#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso, Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from os import makedirs
from os.path import exists, splitext, split, isdir
from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option, get_pynast_version, load_qiime_config,
                        create_dir)
from qiime.align_seqs import alignment_module_names, alignment_method_constructors,\
    pairwise_alignment_methods, CogentAligner, compute_min_alignment_length

options_lookup = get_options_lookup()

qiime_config = load_qiime_config()

blast_db_default_help =\
    qiime_config['pynast_template_alignment_blastdb'] or\
    'created on-the-fly from template_alignment'

alignment_method_choices = \
    alignment_method_constructors.keys() + alignment_module_names.keys()
pairwise_alignment_method_choices = pairwise_alignment_methods.keys()

script_info = {}
script_info[
    'brief_description'] = """Align sequences using a variety of alignment methods"""
script_info['script_description'] = """
This script aligns the sequences in a FASTA file to each other or to a template sequence alignment, depending on the method chosen. Currently, there are three methods which can be used by the user:

1. PyNAST (Caporaso et al., 2009) - The default alignment method is PyNAST, a python implementation of the NAST alignment algorithm.  The NAST algorithm aligns each provided sequence (the "candidate" sequence) to the best-matching sequence in a pre-aligned database of sequences (the "template" sequence).  Candidate sequences are not permitted to introduce new gap characters into the template database, so the algorithm introduces local mis-alignments to preserve the existing template sequence.

2. MUSCLE (Edgar, 2004) - MUSCLE is an alignment method which stands for MUltiple Sequence Comparison by Log-Expectation.

3. INFERNAL (Nawrocki, Kolbe, & Eddy, 2009) - Infernal ("INFERence of RNA ALignment") is for an alignment method for using RNA structure and sequence similarities.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Alignment with PyNAST""", """The default alignment method is PyNAST, a python implementation of the NAST alignment algorithm. The NAST algorithm aligns each provided sequence (the \"candidate\" sequence) to the best-matching sequence in a pre-aligned database of sequences (the \"template\" sequence). Candidate sequences are not permitted to introduce new gap characters into the template database, so the algorithm introduces local mis-alignments to preserve the existing template sequence. The quality thresholds are the minimum requirements for matching between a candidate sequence and a template sequence. The set of matching template sequences will be searched for a match that meets these requirements, with preference given to the sequence length. By default, the minimum sequence length is 150 and the minimum percent id is 75%. The minimum sequence length is much too long for typical pyrosequencing reads, but was chosen for compatibility with the original NAST tool.

The following command can be used for aligning sequences using the PyNAST method, where we supply the program with a FASTA file of unaligned sequences (i.e. resulting FASTA file from pick_rep_set.py, a FASTA file of pre-aligned sequences (this is the template file, which is typically the Greengenes core set - available from http://greengenes.lbl.gov/), and the results will be written to the directory \"pynast_aligned/\":""", """%prog -i $PWD/unaligned.fna -t $PWD/core_set_aligned.fasta.imputed -o $PWD/pynast_aligned_defaults/"""))
script_info['script_usage'].append(
    ("""""",
     """Alternatively, one could change the minimum sequence length ("-e") requirement and minimum sequence identity (\"-p\"), using the following command:""",
     """%prog -i $PWD/unaligned.fna -t core_set_aligned.fasta.imputed -o $PWD/pynast_aligned/ -e 500 -p 95.0"""))

script_info['script_usage'].append(
    ("""Alignment with MUSCLE""",
     """One could also use the MUSCLE algorithm. The following command can be used to align sequences (i.e. the resulting FASTA file from pick_rep_set.py), where the output is written to the directory \"muscle_alignment/\":""",
     """%prog -i $PWD/unaligned.fna -m muscle -o $PWD/muscle_alignment/"""))

script_info['script_usage'].append(("""Alignment with Infernal""", """An alternative alignment method is to use Infernal. Infernal is similar to the PyNAST method, in that you supply a template alignment, although Infernal has several distinct differences. Infernal takes a multiple sequence alignment with a corresponding secondary structure annotation. This input file must be in Stockholm alignment format. There is a fairly good description of the Stockholm format rules at: http://en.wikipedia.org/wiki/Stockholm_format. Infernal will use the sequence and secondary structural information to align the candidate sequences to the full reference alignment. Similar to PyNAST, Infernal will not allow for gaps to be inserted into the reference alignment. Using Infernal is slower than other methods, and therefore is best used with sequences that do not align well using PyNAST.

The following command can be used for aligning sequences using the Infernal method, where we supply the program with a FASTA file of unaligned sequences, a STOCKHOLM file of pre-aligned sequences and secondary structure (this is the template file - an example file can be obtained from: ftp://ftp.microbio.me/qiime/seed.16s.reference_model.sto.zip), and the results will be written to the directory \"infernal_aligned/\":""", """%prog -m infernal -i $PWD/unaligned.fna -t $PWD/seed.16s.reference_model.sto -o $PWD/infernal_aligned/"""))
script_info['output_description'] = """All aligners will output a fasta file containing the alignment and log file in the directory specified by --output_dir (default <alignment_method>_aligned). PyNAST additionally outputs a failures file, containing the sequences which failed to align. So the result of %prog will be up to three files, where the prefix of each file depends on the user supplied FASTA file:

1. \"..._aligned.fasta\" - This is a FASTA file containing all aligned sequences.

2. \"..._failures.fasta\" - This is a FASTA file containing all sequences which did not meet all the criteria specified. (PyNAST only)

3. \"..._log.txt\" - This is a log file containing information pertaining to the results obtained from a particular method (e.g. BLAST percent identity, etc.)."""

script_info['required_options'] = [
    options_lookup['fasta_as_primary_input']
]


script_info['optional_options'] = []

# Check if PyNAST is installed - if get_pynast_version() returns None,
# PyNAST is not installed.
pynast_installed = get_pynast_version() is not None

if pynast_installed:
    script_info['optional_options'].append(
        make_option('-m', '--alignment_method',
                    type='choice', help='Method for aligning' +
                    ' sequences. Valid choices are: ' +
                    ', '.join(alignment_method_choices) +
                    ' [default: %default]',
                    choices=alignment_method_choices,
                    default='pynast'))
    script_info['optional_options'].append(
        make_option('-a', '--pairwise_alignment_method',
                    type='choice', help='method for performing pairwise ' +
                    'alignment in PyNAST. Valid choices are ' +
                    ', '.join(pairwise_alignment_method_choices) +
                    ' [default: %default]',
                    choices=pairwise_alignment_method_choices,
                    default='uclust'))
    script_info['optional_options'].append(
        make_option('-t', '--template_fp',
                    type='existing_filepath',
                    help='Filepath for template alignment [default: %default]',
                    default=qiime_config['pynast_template_alignment_fp']))
    script_info['optional_options'].append(
        make_option('-e', '--min_length',
                    type='int', help='Minimum sequence ' +
                    'length to include in alignment [default: 75% of the' +
                    ' median input sequence length]',
                    default=-1))
    script_info['optional_options'].append(
        make_option('-p', '--min_percent_id',
                    type='float', help='Minimum percent ' +
                    'sequence identity to closest blast hit to include sequence in' +
                    ' alignment [default: %default]', default=0.75))
    script_info['optional_options'].append(
        make_option('-d', '--blast_db', type='blast_db',
                    dest='blast_db', help='Database to blast against when -m pynast ' +
                    '[default: %s]' % blast_db_default_help,
                    default=qiime_config['pynast_template_alignment_blastdb']))
else:
    alignment_method_choices.remove('pynast')
    script_info['optional_options'].append(
        make_option('-m', '--alignment_method',
                    type='choice', help='Method for aligning' +
                    ' sequences. Valid choices are: ' +
                    ', '.join(alignment_method_choices) +
                    ' [default: %default]',
                    choices=alignment_method_choices,
                    default='muscle'))

script_info['optional_options'] += [

    make_option('--muscle_max_memory', type='int',
                help='Maximum memory allocation for the muscle alignment method ' +
                '(MB) [default: 80% of available memory, as detected by MUSCLE]'),

    make_option('-o', '--output_dir', type='new_dirpath',
                help='Path to store ' +
                'result file [default: <ALIGNMENT_METHOD>_aligned]'),
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.min_percent_id <= 1.0:
        opts.min_percent_id *= 100

    if not (1.0 <= opts.min_percent_id <= 100.0):
        option_parser.error('Minimum percent sequence identity must be' +
                            ' between 1.0 and 100.0: %2.2f' % opts.min_percent_id)

    if not opts.template_fp and opts.alignment_method == 'pynast':
        option_parser.error(
            'PyNAST requires a template alignment to be passed via -t')

    input_seqs_filepath = opts.input_fasta_fp
    alignment_method = opts.alignment_method
    output_dir = opts.output_dir or alignment_method + '_aligned'
    create_dir(output_dir, fail_on_exist=False)

    # compute the minimum alignment length if a negative value was
    # provided (the default)
    min_length = opts.min_length
    if min_length < 0:
        min_length = compute_min_alignment_length(
            open(input_seqs_filepath, 'U'))

    fpath, ext = splitext(input_seqs_filepath)
    input_dir, fname = split(fpath)

    result_path = output_dir + '/' + fname + "_aligned.fasta"
    log_path = output_dir + '/' + fname + "_log.txt"
    failure_path = output_dir + '/' + fname + "_failures.fasta"

    if alignment_method in alignment_method_constructors:
        # try/except was causing problems here, so replacing with
        # an explicit check
        # define the aligner params
        aligner_constructor =\
            alignment_method_constructors[alignment_method]
        aligner_type = alignment_method
        params = {'min_len': min_length,
                  'min_pct': opts.min_percent_id,
                  'template_filepath': opts.template_fp,
                  'blast_db': opts.blast_db,
                  'pairwise_alignment_method': opts.pairwise_alignment_method}
        # build the aligner object
        aligner = aligner_constructor(params)
        # apply the aligner
        aligner(input_seqs_filepath, result_path=result_path,
                log_path=log_path, failure_path=failure_path)
    else:
        # define the aligner params
        aligner = CogentAligner({
            'Module': alignment_module_names.get(alignment_method, 'Unknown'),
            'Method': alignment_method
        })
        if alignment_method == "muscle":
            if opts.muscle_max_memory is not None:
                aligner.Params["-maxmb"] = str(opts.muscle_max_memory)
        # build the aligner
        aligner_type = 'Cogent'
        # apply the aligner
        aligner(result_path, seq_path=input_seqs_filepath,
                log_path=log_path)

if __name__ == "__main__":
    main()
