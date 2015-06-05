#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso, Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Jens Reeder",
               "Jose Antonio Navas Molina"]
__credits__ = [
    "Greg Caporaso",
    "Daniel McDonald",
    "Jens Reeder",
    "William Walters"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import split, splitext, join, abspath
from multiprocessing import cpu_count

from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option, create_dir)
from qiime.identify_chimeric_seqs import (blast_fragments_identify_chimeras,
                                          chimeraSlayer_identify_chimeras, usearch61_chimera_check)

options_lookup = get_options_lookup()

# identify_chimeric_seqs.py
script_info = {}
script_info[
    'brief_description'] = """Identify chimeric sequences in input FASTA file"""
script_info['script_description'] = """A FASTA file of sequences, can be screened to remove chimeras (sequences generated due to the PCR amplification of multiple templates or parent sequences). QIIME currently includes a taxonomy-assignment-based approach, blast_fragments, for identifying sequences as chimeric and the ChimeraSlayer algorithm.

1. Blast_fragments approach:

The reference sequences (-r) and id-to-taxonomy map (-t) provided are the same format as those provided to assign_taxonomy.py. The reference sequences are in fasta format, and the id-to-taxonomy map contains tab-separated lines where the first field is a sequence identifier, and the second field is the taxonomy separated by semi-colons (e.g., Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium). The reference collection should be derived from a chimera-checked database (such as the full greengenes database), and filtered to contain only sequences at, for example, a maximum of 97% sequence identity.

2. ChimeraSlayer:

ChimeraSlayer uses BLAST to identify potential chimera parents and computes the optimal branching alignment of the query against two parents.
We suggest to use the pynast aligned representative sequences as input.

3. usearch61:

usearch61 performs both de novo (abundance based) chimera and reference based detection.  Unlike the other two chimera checking software, unclustered sequences should be used as input rather than a representative sequence set, as these sequences need to be clustered to get abundance data.  The results can be taken as the union or intersection of all input sequences not flagged as chimeras.  For details, see: http://drive5.com/usearch/usearch_docs.html
"""

script_info['script_usage'] = []
script_info['script_usage'].append(("""blast_fragments example""", """For each sequence provided as input, the blast_fragments method splits the input sequence into n roughly-equal-sized, non-overlapping fragments, and assigns taxonomy to each fragment against a reference database. The BlastTaxonAssigner (implemented in assign_taxonomy.py) is used for this. The taxonomies of the fragments are compared with one another (at a default depth of 4), and if contradictory assignments are returned the sequence is identified as chimeric. For example, if an input sequence was split into 3 fragments, and the following taxon assignments were returned:

==========  ==========================================================
fragment1:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
fragment2:  Archaea;Euryarchaeota;Halobacteriales;uncultured
fragment3:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
==========  ==========================================================

The sequence would be considered chimeric at a depth of 3 (Methanobacteriales vs. Halobacteriales), but non-chimeric at a depth of 2 (all Euryarchaeota).

blast_fragments begins with the assumption that a sequence is non-chimeric, and looks for evidence to the contrary. This is important when, for example, no taxonomy assignment can be made because no blast result is returned. If a sequence is split into three fragments, and only one returns a blast hit, that sequence would be considered non-chimeric. This is because there is no evidence (i.e., contradictory blast assignments) for the sequence being chimeric. This script can be run by the following command, where the resulting data is written to the directory "identify_chimeras/" and using default parameters (e.g. chimera detection method ("-m blast_fragments"), number of fragments ("-n 3"), taxonomy depth ("-d 4") and maximum E-value ("-e 1e-30")):""", """%prog -i repr_set_seqs.fasta -t taxonomy_assignment.txt -r ref_seq_set.fna -m blast_fragments -o chimeric_seqs_blast.txt"""))

script_info[
    'script_usage'].append(("""ChimeraSlayer Example:""", """Identify chimeric sequences using the ChimeraSlayer algorithm against a user provided reference data base. The input sequences need to be provided in aligned (Py)Nast format. The reference data base needs to be provided as aligned FASTA (-a). Note that the reference database needs to be the same that was used to build the alignment of the input sequences!""",
                            """%prog -m ChimeraSlayer -i repr_set_seqs_aligned.fasta -a ref_seq_set_aligned.fasta -o chimeric_seqs_cs.txt"""))

script_info[
    'script_usage'].append(("""usearch61 Example:""", """Identify chimeric sequences using the usearch61 algorithm against a user provided reference data base.  The input sequences should be the demultiplexed (not clustered rep set!) sequences, such as those output from split_libraries.py. The input sequences need to be provided as unaligned fasta in the same orientation as the query sequences.""",
                            """%prog -m usearch61 -i seqs.fna -r ref_sequences.fasta -o usearch61_chimera_checking/"""))

script_info[
    'output_description'] = """The result of identify_chimeric_seqs.py is a text file that identifies which sequences are chimeric."""
script_info['required_options'] = [options_lookup['fasta_as_primary_input']]

chimera_detection_method_choices = ['blast_fragments', 'ChimeraSlayer',
                                    'usearch61']

script_info['optional_options'] = [
    make_option('-t', '--id_to_taxonomy_fp', type='existing_filepath',
                help='Path to tab-delimited file mapping sequences to assigned '
                'taxonomy. Each assigned taxonomy is provided as a comma-separated '
                'list. [default: %default; REQUIRED when method is blast_fragments]'),

    make_option('-r', '--reference_seqs_fp', type='existing_filepath',
                help='Path to reference sequences (used to build a blast db when '
                'method blast_fragments or reference database for usearch61). '
                '[default: %default; REQUIRED when method blast_fragments' +
                ' if no blast_db is provided, suppress requirement for usearch61 '
                'with --suppress_usearch61_ref;]'),

    make_option('-a', '--aligned_reference_seqs_fp', type='existing_filepath',
                help='Path to (Py)Nast aligned reference sequences. '
                'REQUIRED when method ChimeraSlayer [default: %default]'),

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

    make_option('-R', '--min_div_ratio',
                type='float', help='min divergence ratio ' +
                '(passed to ChimeraSlayer). If set to None uses ' +
                'ChimeraSlayer default value. ' +
                ' [default: %default]', default=None),

    make_option('-k', '--keep_intermediates',
                action='store_true', help='Keep intermediate files, ' +
                'useful for debugging ' +
                ' [default: %default]', default=False),

    make_option('--suppress_usearch61_intermediates', action='store_true',
                help='Use to suppress retention of usearch intermediate files/logs.'
                '[default: %default]', default=False),

    make_option('--suppress_usearch61_ref', action='store_true',
                help='Use to suppress reference based chimera detection with usearch61 '
                '[default: %default]', default=False),

    make_option('--suppress_usearch61_denovo', action='store_true',
                help='Use to suppress de novo based chimera detection with usearch61 '
                '[default: %default]', default=False),

    make_option('--split_by_sampleid', action='store_true',
                help='Enable to split sequences by initial SampleID, requires that fasta '
                'be in demultiplexed format, e.g., >Sample.1_0, >Sample.2_1, >Sample.1_2, '
                'with the initial string before first underscore matching SampleIDs. If '
                'not in this format, could cause unexpected errors. [default: %default]',
                default=False),

    make_option('--non_chimeras_retention', default='union',
                help=("usearch61 only - selects "
                      "subsets of sequences detected as non-chimeras to retain after "
                      "de novo and reference based chimera detection.  Options are "
                      "intersection or union.  union will retain sequences that are "
                      "flagged as non-chimeric from either filter, while intersection "
                      "will retain only those sequences that are flagged as non-"
                      "chimeras from both detection methods. [default: %default]"),
                type='string'),

    make_option('--usearch61_minh', default=0.28, help=("Minimum score (h). "
                                                        "Increasing this value tends to reduce the number of false "
                                                        "positives and decrease sensitivity."
                                                        "[default: %default]"), type='float'),

    make_option('--usearch61_xn', default=8.0, help=("Weight of 'no' vote. "
                                                     "Increasing this value tends to the number of false positives "
                                                     "(and also sensitivity). Must be > 1."
                                                     "[default: %default]"), type='float'),

    make_option('--usearch61_dn', default=1.4, help=("Pseudo-count prior for "
                                                     "'no' votes. (n). Increasing this value tends to the number of "
                                                     "false positives (and also sensitivity). Must be > 0."
                                                     "[default: %default]"), type='float'),

    make_option('--usearch61_mindiffs', default=3, help=("Minimum number of "
                                                         "diffs in a segment. Increasing this value tends to reduce the "
                                                         "number of false positives while reducing sensitivity to very "
                                                         "low-divergence chimeras. Must be > 0."
                                                         "[default: %default]"), type='int'),

    make_option('--usearch61_mindiv', default=0.8, help=("Minimum divergence, "
                                                         "i.e. 100% - identity between the query and closest reference "
                                                         "database sequence. Expressed as a percentage, so the default "
                                                         "is 0.8, which allows chimeras that are up to 99.2% similar to "
                                                         "a reference sequence. This value is chosen to improve "
                                                         "sensitivity to very low-divergence chimeras.  Must be > 0."
                                                         "[default: %default]"), type='float'),

    make_option('--usearch61_abundance_skew', default=2.0, help=("Abundance "
                                                                 "skew setting for de novo chimera detection with usearch61. Must "
                                                                 "be > 0."
                                                                 " [default: %default]"), type='float'),

    make_option('--percent_id_usearch61', default=0.97,
                help=("Percent identity threshold for clustering "
                      "with usearch61, expressed as a fraction between 0 and "
                      "1. [default: %default]"), type='float'),

    make_option('--minlen', default=64, help=("Minimum length of sequence "
                "allowed for usearch61 [default: %default]"), type='int'),

    make_option('--word_length', default=8,
                help="word length value for usearch61. "
                "[default: %default]", type='int'),

    make_option('--max_accepts', default=1,
                help="max_accepts value to usearch61. "
                "[default: %default]", type='int'),

    make_option('--max_rejects', default=8,
                help="max_rejects value for usearch61.  "
                "[default: %default]", type='int'),

    make_option('-o', '--output_fp', type='new_filepath',
                help='Path to store output, output filepath in the case of '
                'blast_fragments and ChimeraSlayer, or directory in case of usearch61 '
                ' [default: derived from input_seqs_fp]'),

    make_option('--threads', default='one_per_cpu', help=(
                "Specify number of threads per core to be used for  "
                "usearch61 commands that utilize multithreading. By default, "
                "will calculate the number of cores to utilize so a single "
                "thread will be used per CPU. Specify a fractional number, e.g."
                " 1.0 for 1 thread per core, or 0.5 for a single thread on "
                "a two core CPU. Only applies to usearch61. "
                "[default: %default]"))
]
script_info['version'] = __version__


def main():
    """Run chimera checker with given options>"""

    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # additional option checks
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
            option_parser.error("Must provide --aligned_reference_seqs_fp "
                                "when using method ChimeraSlayer")
    elif opts.chimera_detection_method == 'usearch61':
        if opts.suppress_usearch61_ref and opts.suppress_usearch61_denovo:
            option_parser.error("Supressing both de novo and reference "
                                "chimera detection not allowed.")
        if not opts.reference_seqs_fp and not opts.suppress_usearch61_ref:
            option_parser.error("--reference_seqs_fp required for reference "
                                "based chimera detection, suppress reference based chimera "
                                "detection with --suppress_usearch61_ref")
        if opts.reference_seqs_fp:
            try:
                temp_f = open(opts.reference_seqs_fp, "U")
                temp_f.close()
            except IOError:
                raise IOError("Unable to open --reference_seqs_fp, please "
                              "check filepath and permissions.")
        if opts.non_chimeras_retention not in ['intersection', 'union']:
            option_parser.error("--non_chimeras_retention must be either "
                                "'union' or 'intersection'")
        if opts.usearch61_xn <= 1:
            option_parser.error("--usearch61_xn must be > 1")
        if opts.usearch61_dn <= 0:
            option_parser.error("--usearch61_dn must be > 0")
        if opts.usearch61_mindiffs <= 0:
            option_parser.error("--usearch61_mindiffs must be > 0")
        if opts.usearch61_mindiv <= 0:
            option_parser.error("--usearch61_mindiv must be > 0")
        if opts.usearch61_abundance_skew <= 0:
            option_parser.error("--usearch61_abundance_skew must be > 0")

    verbose = opts.verbose  # not used yet ...
    input_seqs_fp = opts.input_fasta_fp
    id_to_taxonomy_fp = opts.id_to_taxonomy_fp
    reference_seqs_fp = opts.reference_seqs_fp
    chimera_detection_method = opts.chimera_detection_method
    num_fragments = opts.num_fragments
    output_fp = opts.output_fp
    taxonomy_depth = opts.taxonomy_depth
    max_e_value = opts.max_e_value
    blast_db = opts.blast_db
    keep_intermediates = opts.keep_intermediates
    threads = opts.threads

    # calculate threads as 1 per CPU, or use float of input value
    if threads == 'one_per_cpu':
        threads = float(1 / cpu_count())
    else:
        # Make sure input is a float
        try:
            threads = float(threads)
        except ValueError:
            option_parser.error("--threads must be a float value if "
                                "default 'one_per_cpu' value overridden.")

    if not output_fp:
        if chimera_detection_method == "usearch61":
            output_dir = "usearch61_chimeras/"
            create_dir(output_dir, fail_on_exist=False)
        else:
            input_basename = splitext(split(input_seqs_fp)[1])[0]
            output_fp = '%s_chimeric.txt' % input_basename
    elif chimera_detection_method == "usearch61":
        output_dir = output_fp
        create_dir(output_dir, fail_on_exist=False)

    if chimera_detection_method == 'blast_fragments':
        blast_fragments_identify_chimeras(input_seqs_fp,
                                          id_to_taxonomy_fp,
                                          reference_seqs_fp, blast_db=blast_db,
                                          num_fragments=opts.num_fragments,
                                          max_e_value=max_e_value,
                                          output_fp=output_fp,
                                          taxonomy_depth=taxonomy_depth)
    elif chimera_detection_method == 'ChimeraSlayer':
        chimeraSlayer_identify_chimeras(input_seqs_fp,
                                        output_fp=output_fp,
                                        db_FASTA_fp=opts.reference_seqs_fp,
                                        db_NAST_fp=opts.aligned_reference_seqs_fp,
                                        min_div_ratio=opts.min_div_ratio,
                                        keep_intermediates=keep_intermediates)
    elif chimera_detection_method == 'usearch61':
        usearch61_chimera_check(input_seqs_fp,
                                output_dir=output_dir,
                                reference_seqs_fp=reference_seqs_fp,
                                suppress_usearch61_intermediates=opts.suppress_usearch61_intermediates,
                                suppress_usearch61_ref=opts.suppress_usearch61_ref,
                                suppress_usearch61_denovo=opts.suppress_usearch61_denovo,
                                split_by_sampleid=opts.split_by_sampleid,
                                non_chimeras_retention=opts.non_chimeras_retention,
                                usearch61_minh=opts.usearch61_minh,
                                usearch61_xn=opts.usearch61_xn,
                                usearch61_dn=opts.usearch61_dn,
                                usearch61_mindiffs=opts.usearch61_mindiffs,
                                usearch61_mindiv=opts.usearch61_mindiv,
                                usearch61_abundance_skew=opts.usearch61_abundance_skew,
                                percent_id_usearch61=opts.percent_id_usearch61,
                                minlen=opts.minlen,
                                word_length=opts.word_length,
                                max_accepts=opts.max_accepts,
                                max_rejects=opts.max_rejects,
                                verbose=opts.verbose,
                                threads=threads)


if __name__ == "__main__":
    main()
