#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger",
               "Jens Reeder", "William Walters", "Jose Carlos Clemente Litran",
               "Jai Ram Rideout", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import splitext, split, exists, abspath, isfile
from os import makedirs
from multiprocessing import cpu_count

from skbio.util import remove_files
from qiime.util import (parse_command_line_parameters, create_dir,
                         make_option, load_qiime_config)
from qiime.sort import sort_fasta_by_abundance
from qiime.pick_otus  import otu_picking_method_constructors,\
    otu_picking_method_choices, MothurOtuPicker

qiime_config = load_qiime_config()

script_info = {}
script_info['brief_description'] = """OTU picking"""
script_info['script_description'] = """The OTU picking step assigns similar sequences to operational taxonomic units, or OTUs, by clustering sequences based on a user-defined similarity threshold. Sequences which are similar at or above the threshold level are taken to represent the presence of a taxonomic unit (e.g., a genus, when the similarity threshold is set at 0.94) in the sequence collection.

Currently, the following clustering methods have been implemented in QIIME:

1.  cd-hit (Li & Godzik, 2006; Li, Jaroszewski, & Godzik, 2001), which applies a \"longest-sequence-first list removal algorithm\" to cluster sequences.

2.  blast (Altschul, Gish, Miller, Myers, & Lipman, 1990), which compares and clusters each sequence against a reference database of sequences.

3.  Mothur (Schloss et al., 2009), which requires an input file of aligned sequences.  The input file of aligned sequences may be generated from an input file like the one described below by running align_seqs.py.  For the Mothur method, the clustering algorithm may be specified as nearest-neighbor, furthest-neighbor, or average-neighbor.  The default algorithm is furthest-neighbor.

4.  prefix/suffix [Qiime team, unpublished], which will collapse sequences which are identical in their first and/or last bases (i.e., their prefix and/or suffix). The prefix and suffix lengths are provided by the user and default to 50 each.

5.  Trie [Qiime team, unpublished], which collapsing identical sequences and sequences which are subsequences of other sequences.

6.  uclust (Edgar, RC 2010), creates \"seeds\" of sequences which generate clusters based on percent identity.

7.  uclust_ref (Edgar, RC 2010), as uclust, but takes a reference database to use as seeds.  New clusters can be toggled on or off.

8.  usearch (Edgar, RC 2010, version v5.2.236), creates \"seeds\" of sequences which generate clusters based on percent identity, filters low abundance clusters, performs de novo and reference based chimera detection.

9.  usearch_ref (Edgar, RC 2010, version v5.2.236), as usearch, but takes a reference database to use as seeds.  New clusters can be toggled on or off.

Quality filtering pipeline with usearch 5.X is described as usearch_qf "usearch quality filter", described here: http://qiime.org/tutorials/usearch_quality_filter.html

8.  usearch61 (Edgar, RC 2010, version v6.1.544), creates \"seeds\" of sequences which generate clusters based on percent identity.

9.  usearch61_ref (Edgar, RC 2010, version v6.1.544), as usearch61, but takes a reference database to use as seeds.  New clusters can be toggled on or off.

10. sumaclust (Mercier, C. et al., 2014, version 1.0), creates \"seeds\" of sequences which generate clusters based on similarity threshold.

11. sortmerna_v2 (Kopylova, E. et al., 2012), takes a reference database to use as seeds.

12. swarm (Mahe, F. et al., 2014), creates \"seeds\" of sequences which generate clusters based on a resolution threshold.


Chimera checking with usearch 6.X is implemented in identify_chimeric_seqs.py.  Chimera checking should be done first with usearch 6.X, and the filtered resulting fasta file can then be clustered.


The primary inputs for pick_otus.py are:

1. A FASTA file containing sequences to be clustered

2. An OTU threshold (default is 0.97, roughly corresponding to species-level OTUs);

3. The method to be applied for clustering sequences into OTUs.

pick_otus.py takes a standard fasta file as input.

"""

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""Example (uclust method, default):""",
     """Using the seqs.fna file generated from split_libraries.py and outputting the results to the directory \"picked_otus_default/\", while using default parameters (0.97 sequence similarity, no reverse strand matching):""",
     """%prog -i seqs.fna -o picked_otus_default"""))

script_info['script_usage'].append(
    ("""""",
     """To change the percent identity to a lower value, such as 90%, and also enable reverse strand matching, the command would be the following:""",
     """%prog -i seqs.fna -o picked_otus_90_percent_rev/ -s 0.90 -z"""))

script_info['script_usage'].append(
    ("""Uclust Reference-based OTU picking example""",
     """uclust_ref can be passed via -m to pick OTUs against a reference set where sequences within the similarity threshold to a reference sequence will cluster to an OTU defined by that reference sequence, and sequences outside of the similarity threshold to a reference sequence will form new clusters. OTU identifiers will be set to reference sequence identifiers when sequences cluster to reference sequences, and 'qiime_otu_<integer>' for new OTUs. Creation of new clusters can be suppressed by passing -C, in which case sequences outside of the similarity threshold to any reference sequence will be listed as failures in the log file, and not included in any OTU.""",
     """%prog -i seqs.fna -r refseqs.fasta -m uclust_ref --denovo_otu_id_prefix qiime_otu_"""))

script_info['script_usage'].append(
    ("""Example (cdhit method):""",
     """Using the seqs.fna file generated from split_libraries.py and outputting the results to the directory \"cdhit_picked_otus/\", while using default parameters (0.97 sequence similarity, no prefix filtering):""",
     """%prog -i seqs.fna -m cdhit -o cdhit_picked_otus/"""))

script_info['script_usage'].append(
    ("""""",
     """Currently the cd-hit OTU picker allows for users to perform a pre-filtering step, so that highly similar sequences are clustered prior to OTU picking. This works by collapsing sequences which begin with an identical n-base prefix, where n is specified by the -n parameter. A commonly used value here is 100 (e.g., -n 100). So, if using this filter with -n 100, all sequences which are identical in their first 100 bases will be clustered together, and only one representative sequence from each cluster will be passed to cd-hit. This is used to greatly decrease the run-time of cd-hit-based OTU picking when working with very large sequence collections, as shown by the following command:""",
     """%prog -i seqs.fna -m cdhit -o cdhit_picked_otus_filter/ -n 100"""))

script_info['script_usage'].append(("""""", """Alternatively, if the user would like to collapse identical sequences, or those which are subsequences of other sequences prior to OTU picking, they can use the trie prefiltering (\"-t\") option as shown by the following command.

Note: It is highly recommended to use one of the prefiltering methods when analyzing large datasets (>100,000 seqs) to reduce run-time.""", """%prog -i seqs.fna -m cdhit -o cdhit_picked_otus_trie_prefilter/ -t"""))

script_info['script_usage'].append(("""BLAST OTU-Picking Example:""", """OTUs can be picked against a reference database using the BLAST OTU picker. This is useful, for example, when different regions of the SSU RNA have sequenced and a sequence similarity based approach like cd-hit therefore wouldn't work. When using the BLAST OTU picking method, the user must supply either a reference set of sequences or a reference database to compare against. The OTU identifiers resulting from this step will be the sequence identifiers in the reference database. This allows for use of a pre-existing tree in downstream analyses, which again is useful in cases where different regions of the 16s gene have been sequenced.

The following command can be used to blast against a reference sequence set, using the default E-value and sequence similarity (0.97) parameters:""", """%prog -i seqs.fna -o blast_picked_otus/ -m blast -r refseqs.fasta"""))

script_info['script_usage'].append(
    ("""""",
     """If you already have a pre-built BLAST database, you can pass the database prefix as shown by the following command:""",
     """%prog -i seqs.fna -o blast_picked_otus_prebuilt_db/ -m blast -b refseqs.fasta"""))

script_info['script_usage'].append(
    ("""""",
     """If the user would like to change the sequence similarity (\"-s\") and/or the E-value (\"-e\") for the blast method, they can use the following command:""",
     """%prog -i seqs.fna -o blast_picked_otus_90_percent/ -m blast -r refseqs.fasta -s 0.90 -e 1e-30"""))

script_info['script_usage'].append(
    ("""Prefix-suffix OTU Picking Example:""",
     """OTUs can be picked by collapsing sequences which begin and/or end with identical bases (i.e., identical prefixes or suffixes).  This OTU picker is currently likely to be of limited use on its own, but will be very useful in collapsing very similar sequences in a chained OTU picking strategy that is currently in development. For example, the user will be able to pick OTUs with this method, followed by representative set picking, and then re-pick OTUs on their representative set. This will allow for highly similar sequences to be collapsed, followed by running a slower OTU picker. This ability to chain OTU pickers is not yet supported in QIIME. The following command illustrates how to pick OTUs by collapsing sequences which are identical in their first 50 and last 25 bases:""",
     """%prog -i seqs.fna -o prefix_suffix_picked_otus/ -m prefix_suffix -p 50 -u 25"""))

script_info['script_usage'].append(("""Mothur OTU Picking Example:""", """The Mothur program (http://www.mothur.org/) provides three clustering algorithms for OTU formation: furthest-neighbor (complete linkage), average-neighbor (group average), and nearest-neighbor (single linkage). Details on the algorithms may be found on the Mothur website and publications (Schloss et al., 2009). However, the running times of Mothur's clustering algorithms scale with the number of sequences squared, so the program may not be feasible for large data sets.

The following command may be used to create OTUs based on a furthest-neighbor algorithm (the default setting) using aligned sequences as input:""", """%prog -i seqs.aligned.fna -o mothur_picked_otus/ -m mothur"""))

script_info['script_usage'].append(
    ("""""",
     """If you prefer to use a nearest-neighbor algorithm instead, you may specify this with the '-c' flag:""",
     """%prog -i seqs.aligned.fna -o mothur_picked_otus_nn/ -m mothur -c nearest"""))

script_info['script_usage'].append(
    ("""""",
     """The sequence similarity parameter may also be specified. For example, the following command may be used to create OTUs at the level of 90% similarity:""",
     """%prog -i seqs.aligned.fna -o mothur_picked_otus_90_percent/ -m mothur -s 0.90"""))

script_info['script_usage'].append(
    ("""usearch """,
     """Usearch (http://www.drive5.com/usearch/) provides clustering, chimera checking, and quality filtering. The following command specifies a minimum cluster size of 2 to be used during cluster size filtering:""",
     """%prog -i seqs.fna -m usearch --word_length 64 --db_filepath refseqs.fasta -o usearch_qf_results/ --minsize 2"""))

script_info['script_usage'].append(
    ("""usearch example where reference-based chimera detection is disabled, and minimum cluster size filter is reduced from default (4) to 2:""",
     """""",
     """%prog -i seqs.fna -m usearch --word_length 64 --suppress_reference_chimera_detection --minsize 2 -o usearch_qf_results_no_ref_chim_detection/"""))

script_info['script_usage'].append(
    ("Use de novo OTU-picker Swarm:",
     "Using the seqs.fna file generated from split_libraries.py and "
     "outputting the results to the directory \"$PWD/picked_otus_swarm/\", "
     "while using default parameters (resolution = 1) ",
     "%prog -i $PWD/seqs.fna -m swarm -o $PWD/picked_otus_swarm"))

script_info['script_usage_output_to_remove'] = ['$PWD/picked_otus_swarm/']

script_info['output_description'] = """The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an (arbitrary) cluster identifier, and the remaining fields correspond to sequence identifiers assigned to that cluster. Sequence identifiers correspond to those provided in the input FASTA file.  Usearch (i.e. usearch quality filter) can additionally have log files for each intermediate call to usearch.

Example lines from the resulting .txt file:

=   ====    ====    ====
0   seq1    seq5
1   seq2
2   seq3
3   seq4    seq6    seq7
=   ====    ====    ====

This result implies that four clusters were created based on 7 input sequences. The first cluster (cluster id 0) contains two sequences, sequence ids seq1 and seq5; the second cluster (cluster id 1) contains one sequence, sequence id seq2; the third cluster (cluster id 2) contains one sequence, sequence id seq3, and the final cluster (cluster id 3) contains three sequences, sequence ids seq4, seq6, and seq7.

The resulting .log file contains a list of parameters passed to the pick_otus.py script along with the output location of the resulting .txt file."""

script_info['required_options'] = [
    make_option('-i', '--input_seqs_filepath', type='existing_filepath',
                help='Path to input sequences file'),
]

script_info['optional_options'] = [
    make_option('-m', '--otu_picking_method', type='choice',
                choices=otu_picking_method_choices, default="uclust",
                help='Method for picking OTUs.  Valid choices are: ' +
                      ', '.join(otu_picking_method_choices) +
                      '. The mothur method requires an input file '
                      'of aligned sequences.  usearch will enable the usearch quality '
                      'filtering pipeline. [default: %default]'),

    make_option('-c', '--clustering_algorithm', type='choice',
                choices=MothurOtuPicker.ClusteringAlgorithms, default='furthest',
                help='Clustering algorithm for mothur otu picking method.  Valid '
                      'choices are: ' +
                      ', '.join(MothurOtuPicker.ClusteringAlgorithms) +
                      '. [default: %default]'),

    make_option('-M', '--max_cdhit_memory', type='int', default=400,
                help='Maximum available memory to cd-hit-est (via the program\'s -M '
                      'option) for cdhit OTU picking method (units of Mbyte) '
                      '[default: %default]'),

    make_option('-o', '--output_dir', type='new_dirpath',
                help='Path to store result file '
                      '[default: ./<OTU_METHOD>_picked_otus/]'),

    make_option('-r', '--refseqs_fp', type='existing_filepath',
                help='Path to reference sequences to search against when using -m '
                      'blast, -m sortmerna, -m uclust_ref, -m usearch_ref, or -m '
                      'usearch61_ref [default: %default]',
                      default=qiime_config['pick_otus_reference_seqs_fp']),

    make_option('-b', '--blast_db', type='blast_db',
                help='Pre-existing database to blast against when using -m blast '
                      '[default: %default]'),

    make_option('-e', '--max_e_value_blast', type='float', default=1e-10,
                help='Max E-value when clustering with BLAST '
                     '[default: %default]'),

    # SortMeRNA specific parameters
    make_option('--sortmerna_db', type='string',
                help='Pre-existing database to search against when using '
                     '-m sortmerna [default: %default]'),

    make_option('--sortmerna_e_value', type='float', default=1,
                help='Maximum E-value when clustering [default = %default]'),

    make_option('--sortmerna_coverage', type='float', default=0.97,
                help='Mininum percent query coverage (of an alignment) '
                     'to consider a hit, expressed as a fraction between 0 '
                     'and 1 [default: %default]'),

    make_option('--sortmerna_tabular', default=False, action='store_true',
                help='Output alignments in the Blast tabular format '
                     'with two additional columns including the CIGAR '
                     'string and the percent query coverage '
                     '[default: %default]'),

    make_option('--sortmerna_best_N_alignments', type='int', default=1,
                help='Must be set together with --sortmerna_tabular. '
                     'This option specifies how many alignments per read '
                     'will be written [default: %default]'),

    make_option('--sortmerna_max_pos', type='int', default=10000,
                help='The maximum number of positions per seed to store '
                     ' in the indexed database [default: %default]'),
    # end SortMeRNA specific parameters

    make_option('--min_aligned_percent',
                help='Minimum percent of query sequence that can be aligned '
                     'to consider a hit, expressed as a fraction between 0 '
                     'and 1 (BLAST OTU picker only) '
                     '[default: %default]',
                default=0.50, type='float'),

    make_option('-s', '--similarity', type='float', default=0.97,
                help=('Sequence similarity threshold (for blast, cdhit, uclust, '
                      'uclust_ref, usearch, usearch_ref, usearch61, usearch61_ref, '
                      'sumaclust, and sortmerna) [default: %default]')),

    make_option('--sumaclust_exact', action='store_true', default=False,
                help='A sequence is assigned to the best matching seed '
                     'rather than the first matching seed passing the '
                     'similarity threshold [default: %default]'),

    make_option('--sumaclust_l', action='store_true', default=True,
                help='Reference sequence length if the shortest '
                     '[default: %default]'),

    make_option('--denovo_otu_id_prefix', default="denovo", type='string',
                help='OTU identifier prefix (string) for the de novo '
                     'OTU pickers (sumaclust, swarm and uclust) '
                     '[default: %default, OTU ids are ascending'
                     'integers]'),

    # Swarm specific parameters
    make_option('--swarm_resolution', default=1, type='int',
                help='Maximum number of differences allowed between '
                     'two amplicons, meaning that two amplicons will '
                     'be grouped if they have integer (or less) '
                     'differences (see Swarm manual at '
                     'https://github.com/torognes/swarm for more details). '
                     '[default: %default]'),
    # end Swarm specific parameters

    make_option('-q', '--trie_reverse_seqs', action='store_true',
                default=False,
                help='Reverse seqs before picking OTUs with the Trie OTU picker for '
                      'suffix (rather than prefix) collapsing [default: %default]'),

    make_option('-n', '--prefix_prefilter_length', type='int', default=None,
                help='Prefilter data so seqs with identical first '
                      'prefix_prefilter_length are automatically grouped into a single '
                      'OTU.  This is useful for large sequence collections where OTU '
                      'picking doesn\'t scale well [default: %default; 100 is a good '
                      'value]'),

    make_option('-t', '--trie_prefilter', action='store_true',
                default=False,
                help='prefilter data so seqs which are identical prefixes of a longer '
                      'seq are automatically grouped into a single OTU; useful for '
                      'large sequence collections where OTU picking doesn\'t scale '
                      'well [default: %default]'),

    make_option('-p', '--prefix_length', type='int', default=50,
                help='Prefix length when using the prefix_suffix otu picker; '
                      'WARNING: CURRENTLY DIFFERENT FROM prefix_prefilter_length '
                      '(-n)! [default: %default]'),

    make_option('-u', '--suffix_length', type='int', default=50,
                help='Suffix length when using the prefix_suffix otu picker '
                      '[default: %default]'),

    make_option('-z', '--enable_rev_strand_match', action='store_true',
                default=False,
                help='Enable reverse strand matching for uclust, uclust_ref, '
                      'usearch, usearch_ref, usearch61, or usearch61_ref otu picking, '
                      'will double the amount of memory used. [default: %default]'),

    make_option('-D', '--suppress_presort_by_abundance_uclust',
                action='store_true',
                default=False,
                help='Suppress presorting of sequences by abundance when picking'
                      ' OTUs with uclust or uclust_ref [default: %default]'),

    make_option('-A', '--optimal_uclust', action='store_true',
                default=False,
                help='Pass the --optimal flag to uclust for uclust otu'
                      ' picking. [default: %default]'),

    make_option('-E', '--exact_uclust', action='store_true',
                default=False,
                help='Pass the --exact flag to uclust for uclust otu'
                      ' picking. [default: %default]'),

    make_option('-B', '--user_sort', action='store_true',
                default=False,
                help='Pass the --user_sort flag to uclust for uclust otu'
                      ' picking. [default: %default]'),

    make_option('-C', '--suppress_new_clusters', action='store_true',
                default=False,
                help="Suppress creation of new clusters using seqs that don't"
                " match reference when using -m uclust_ref, -m usearch61_ref, or "
                "-m usearch_ref [default: %default]"),

    make_option('--max_accepts', default='default',
                help="max_accepts value to uclust, uclust_ref, usearch61, and "
                "usearch61_ref.  By default, will use value suggested by "
                "method (uclust: 1, usearch61: 1) [default: %default]"),

    make_option('--max_rejects', default='default',
                help="max_rejects value for uclust, uclust_ref, usearch61, and "
                "usearch61_ref.  With default settings, will use value "
                "recommended by clustering method used "
                "(uclust: 8, usearch61: 8 for usearch_fast_cluster option,"
                " 32 for reference and smallmem options) "
                "[default: %default]"),

    make_option('--stepwords', type='int', default=8,
                help="stepwords value to uclust and "
                "uclust_ref [default: %default]"),

    make_option('--word_length', default='default',
                help="word length value for uclust, uclust_ref, and "
                "usearch, usearch_ref, usearch61, and usearch61_ref. "
                "With default setting, will use the setting recommended by "
                "the method (uclust: 8, usearch: 64, usearch61: 8).  int "
                "value can be supplied to override this setting. "
                "[default: %default]"),

    make_option('--suppress_uclust_stable_sort', default=False,
                action='store_true', help="Don't pass --stable-sort to "
                                           "uclust [default: %default]"),

    make_option('--suppress_prefilter_exact_match',
                default=False, action='store_true',
                help="Don't collapse exact matches before calling "
                  "sortmerna, sumaclust or uclust [default: %default]"),

    make_option('-d', '--save_uc_files', default=True, action='store_false',
                help="Enable preservation of intermediate uclust (.uc) files "
                      "that are used to generate clusters via uclust.  Also enables "
                      "preservation of all intermediate files created by usearch "
                      " and usearch61. [default: %default]"),

    make_option('-j', '--percent_id_err', default=0.97,
                help="Percent identity threshold for cluster error detection "
                      "with usearch, expressed as a fraction between 0 and "
                      "1. [default: %default]", type='float'),

    make_option('-g', '--minsize', default=4, help="Minimum cluster size "
                                                    "for size filtering with usearch. [default: %default]",
                type='int'),

    make_option('-a', '--abundance_skew', default=2.0, help="Abundance skew "
                                                             "setting for de novo chimera detection with usearch. "
                                                             "[default: %default]", type='float'),

    make_option('-f', '--db_filepath', type='existing_filepath', default=None,
                help="Reference database of fasta sequences for reference "
                      "based chimera detection with usearch. [default: %default]"),

    make_option('--perc_id_blast', default=0.97, help="Percent ID for "
                                                       "mapping OTUs created by usearch back to original sequence"
                                                       " IDs [default: %default]", type='float'),

    make_option('--de_novo_chimera_detection', help="Deprecated:  de novo chimera detection performed by default, "
        "pass --suppress_de_novo_chimera_detection to disable."
        " [default: %default]"),

    make_option('-k', '--suppress_de_novo_chimera_detection', default=False,
                help="Suppress de novo chimera detection in usearch. "
                      "[default: %default]", action='store_true'),

    make_option('--reference_chimera_detection',
                help="Deprecated:  Reference based chimera detection performed "
                      "by default, pass --supress_reference_chimera_detection to "
                      "disable [default: %default]"),

    make_option('-x', '--suppress_reference_chimera_detection', default=False,
                help="Suppress reference based chimera detection in usearch. "
                      "[default: %default]", action='store_true'),

    make_option('--cluster_size_filtering', help="Deprecated, "
                                                  "cluster size filtering enabled by default, pass "
                                                  "--suppress_cluster_size_filtering to disable."
                                                  "  [default: %default]"),

    make_option('-l', '--suppress_cluster_size_filtering', default=False,
                help="Suppress cluster size filtering in usearch.  "
                      "[default: %default]", action='store_true'),

    make_option('--remove_usearch_logs', default=False, help="Disable "
                                                              "creation of logs when usearch is called.  Up to nine logs are "
                                                              "created, depending on filtering steps enabled.  "
                                                              "[default: %default]", action='store_true'),

    make_option('--derep_fullseq', default=False, help="Dereplication "
                "of full sequences, instead of subsequences. Faster than "
                                                        "the default --derep_subseqs in usearch. "
                                                        "[default: %default]", action='store_true'),

    make_option('-F', '--non_chimeras_retention', default='union',
                help="Selects "
                      "subsets of sequences detected as non-chimeras to retain after "
                      "de novo and reference based chimera detection.  Options are "
                      "intersection or union.  union will retain sequences that are "
                      "flagged as non-chimeric from either filter, while intersection "
                      "will retain only those sequences that are flagged as non-"
                      "chimeras from both detection methods. [default: %default]",
                type='string'),

    make_option('--minlen', default=64, help="Minimum length of sequence "
                "allowed for usearch, usearch_ref, usearch61, and "
                "usearch61_ref. [default: %default]", type='int'),

    make_option('--usearch_fast_cluster', default=False, help="Use fast "
                "clustering option for usearch or usearch61_ref with new "
                "clusters.  --enable_rev_strand_match can not be enabled "
                "with this option, and the only valid option for "
                "usearch61_sort_method is 'length'.  This option uses more "
                "memory than the default option for de novo clustering."
                " [default: %default]", action='store_true'),

    make_option('--usearch61_sort_method', default='abundance', help=
                "Sorting method for usearch61 and usearch61_ref.  Valid "
                "options are abundance, length, or None.  If the "
                "--usearch_fast_cluster option is enabled, the only sorting "
                "method allowed in length. [default: %default]", type='str'),

    make_option('--sizeorder', default=False, help=
                "Enable size based preference in clustering with usearch61. "
                "Requires that --usearch61_sort_method be abundance. "
                "[default: %default]", action='store_true'),

    make_option('--threads', default=1, help=
                "Specify number of threads (1 thread per core) to be used for usearch61, "
                "sortmerna, sumaclust and swarm commands that utilize multithreading. "
                "[default: %default]")
]

script_info['version'] = __version__


def main():
    # Parse the command line parameters
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create local copies of the options to avoid repetitive lookups
    prefix_prefilter_length = opts.prefix_prefilter_length
    otu_picking_method = opts.otu_picking_method
    prefix_length = opts.prefix_length
    suffix_length = opts.suffix_length
    trie_prefilter = opts.trie_prefilter
    trie_reverse_seqs = opts.trie_reverse_seqs
    refseqs_fp = opts.refseqs_fp
    blast_db = opts.blast_db
    similarity = opts.similarity
    enable_rev_strand_match = opts.enable_rev_strand_match
    suppress_presort_by_abundance_uclust = \
        opts.suppress_presort_by_abundance_uclust
    optimal_uclust = opts.optimal_uclust
    user_sort = opts.user_sort
    max_accepts = opts.max_accepts
    max_rejects = opts.max_rejects
    stepwords = opts.stepwords
    word_length = opts.word_length
    min_aligned_percent = opts.min_aligned_percent
    uclust_stable_sort = not opts.suppress_uclust_stable_sort
    save_uc_files = opts.save_uc_files
    prefilter_identical_sequences =\
        not opts.suppress_prefilter_exact_match
    derep_fullseq = opts.derep_fullseq
    chimeras_retention = opts.non_chimeras_retention
    verbose = opts.verbose
    threads = opts.threads
    denovo_otu_id_prefix = opts.denovo_otu_id_prefix

    # sortmerna specific parameters
    sortmerna_db = opts.sortmerna_db
    sortmerna_e_value = opts.sortmerna_e_value
    sortmerna_coverage = opts.sortmerna_coverage
    sortmerna_tabular = opts.sortmerna_tabular
    sortmerna_best_N_alignments = opts.sortmerna_best_N_alignments
    sortmerna_max_pos = opts.sortmerna_max_pos

    # swarm specific parameters
    swarm_resolution = opts.swarm_resolution

    # usearch specific parameters
    percent_id_err = opts.percent_id_err
    minsize = opts.minsize
    abundance_skew = opts.abundance_skew
    db_filepath = opts.db_filepath
    perc_id_blast = opts.perc_id_blast
    de_novo_chimera_detection = not opts.suppress_de_novo_chimera_detection
    reference_chimera_detection = not opts.suppress_reference_chimera_detection
    cluster_size_filtering = not opts.suppress_cluster_size_filtering
    remove_usearch_logs = opts.remove_usearch_logs
    minlen = opts.minlen

    # usearch61 specific parameters
    # also uses: enable_rev_strand_match, refseqs_fp, suppress_new_clusters,
    # save_uc_files, remove_usearch_logs, minlen, maxaccepts, maxrejects,
    # word_len
    usearch_fast_cluster = opts.usearch_fast_cluster
    usearch61_sort_method = opts.usearch61_sort_method
    sizeorder = opts.sizeorder

    # sumaclust specific parameters
    sumaclust_exact = opts.sumaclust_exact
    sumaclust_l = opts.sumaclust_l

    # Set default values according to clustering method
    if word_length != "default":
        try:
            word_length = int(word_length)
        except ValueError:
            raise ValueError("--word_length must either be 'default' "
                             "or an int value")
    if word_length == "default":
        if otu_picking_method in ["usearch", "usearch_ref"]:
            word_length = 64
        else:
            # default setting for usearch61, uclust, uclust_ref
            word_length = 8

    if max_accepts != "default":
        try:
            max_accepts = int(max_accepts)
        except ValueError:
            option_parser.error("--max_accepts must either be 'default' "
                                "or an int value")
    if max_accepts == "default":
        # default setting for usearch61, uclust, uclust_ref
        max_accepts = 1

    if max_rejects != "default":
        try:
            max_rejects = int(max_rejects)
        except ValueError:
            option_parser.error("--max_rejects must be either 'default' "
                                "or an int value")
    if max_rejects == "default":
        if otu_picking_method in ["uclust", "uclust_ref"]:
            max_rejects = 8
        # usearch61 settings, depends upon fast clustering option
        else:
            if usearch_fast_cluster:
                max_rejects = 8
            else:
                max_rejects = 32

    # Check for logical/compatible inputs
    if user_sort and not suppress_presort_by_abundance_uclust:
        option_parser.error(
            "Cannot pass -B/--user_sort without -D/--suppress_presort_by_abundance_uclust, as your input would be resorted by abundance. To presort your own sequences before passing to uclust, pass -DB.")

    if abundance_skew <= 1:
        option_parser.error('abundance skew must be > 1')

    # Check for logical inputs
    if otu_picking_method in ['usearch', 'usearch_ref'] and \
            reference_chimera_detection and not db_filepath:
        option_parser.error('No reference filepath specified with '
                            '--db_filepath option. Disable reference based chimera detection '
                            'with --suppress_reference_chimera_detection or specify a reference '
                            'fasta file with --db_filepath.')

    if chimeras_retention not in ['intersection', 'union']:
        option_parser.error('--chimeras_retention must be either union or '
                            'intersection.')

    if usearch61_sort_method not in ['length', 'abundance', 'None']:
        option_parser.error("--usearch61_sort_method must be one of the "
                            "following: length, abundance, None")

    if otu_picking_method in ['usearch61', 'usearch61_ref']:
        if usearch_fast_cluster:
            if enable_rev_strand_match:
                option_parser.error("--enable_rev_strand_match can not be "
                                    "enabled when using --usearch_fast_cluster.")
            if usearch61_sort_method != "length":
                option_parser.error("--usearch61_sort_method must be 'length' "
                                    "when --usearch_fast_cluster used.")

    if otu_picking_method in ['usearch61']:
        if opts.suppress_new_clusters:
            option_parser.error("--suppress_new_clusters cannot be enabled when "
                                "using usearch61 as the OTU picking method as this is strictly "
                                "de novo.  Use --otu_picking_method usearch61_ref and a reference "
                                "database for closed reference OTU picking.")

    if sizeorder:
        if usearch61_sort_method != 'abundance':
            option_parser.error("To use --sizeorder, usearch61_sort_method must "
                                "be abundance.")

    # Test that db_filepath can be opened to avoid wasted time
    if db_filepath:
        try:
            tmp_db_filepath = open(db_filepath, "U")
            tmp_db_filepath.close()
            db_filepath = abspath(db_filepath)
        except IOError:
            raise IOError('Unable to open %s, please check path/permissions' %
                          db_filepath)

    # Input validation to throw a useful error message on common mistakes
    if not (0.0 <= similarity <= 1.0):
        option_parser.error("%r is an invalid sequence similarity threshold. "
                            "--similarity must be between 0.0 and 1.0 "
                            "(inclusive)." % similarity)

    if (otu_picking_method == 'cdhit' and
            similarity < 0.80):
        option_parser.error('cdhit requires similarity >= 0.80.')

    if (otu_picking_method == 'blast' and
            refseqs_fp is None and
            blast_db is None):
        option_parser.error('blast requires refseqs_fp or blast_db')

    if otu_picking_method in ['uclust_ref', 'usearch_ref', 'usearch61_ref']:
        if (refseqs_fp is None):
            option_parser.error('uclust_ref, usearch_ref, usearch61_ref ' +
                                ' requires refseqs_fp')
        elif not exists(refseqs_fp):
            option_parser.error('refseqs_fp %s does not exist' % refseqs_fp)
        else:
            refseqs_fp = abspath(refseqs_fp)

    # number of threads to use
    if threads == 1:
        # Make sure input is an integer
        try:
            threads = int(threads)
        except ValueError:
            option_parser.error("--threads must be a integer value.")


    if otu_picking_method == 'sortmerna':

        # check sortmerna_e_value is a float and positive
        try:
            sortmerna_e_value = float(sortmerna_e_value)
        except ValueError:
            option_parser.error("--sortmerna_e_value must be a float.")
        if sortmerna_e_value < 0:
            option_parser.error("--sortmerna_e_value must be positive.")

        # check sortmerna_coverage is a float and in [0,1]
        try:
            sortmerna_coverage = float(sortmerna_coverage)
        except ValueError:
            option_parser.error('--sortmerna_coverage must be a float.')
        if sortmerna_e_value < 0:
            option_parser.error('--sortmerna_coverage must be positive.')

        # check sortmerna_best_N_alignments is an integer
        try:
            sortmerna_best_N_alignments = int(sortmerna_best_N_alignments)
        except ValueError:
            option_parser.error('--sortmerna_best_N_alignments must '
                                'be an integer value.')
        if sortmerna_best_N_alignments < 0:
            option_parser.error('--sortmerna_best_N_alignments must '
                                'be a positive value.')

        # sortmerna_tabular must be set if sortmerna_best_N_alignments > 1;
        # sortmerna_best_N_alignments = 1 will always be passed to sortmerna,
        # with or without sortmerna_tabular, as at least 1 best match is
        # required to build an OTU map
        elif sortmerna_best_N_alignments > 1 and \
             sortmerna_tabular is False:
             option_parser.error('--sortmerna_tabular must be set together '
                                 'with --sortmerna_best_N_alignments.')

        # check FASTA reference file or the indexed database (with the FASTA reference file) were provided
        if refseqs_fp is None:
            option_parser.error('sortmerna always requires refseqs_fp (with or without sortmerna_db)')
        elif sortmerna_db:
            if isfile(sortmerna_db + '.stats') is False:
                option_parser.error('%s does not exist, make sure you have indexed '
                                    'the database using indexdb_rna' % (sortmerna_db + '.stats'))

    if otu_picking_method == 'swarm':
        # check resolution is a positive integer
        if swarm_resolution < 1:
            option_parser.error('--swarm_resolution=INT must '
                                'be a positive integer value.')

    # End input validation

    # use the otu_picking_method value to get the otu picker constructor
    otu_picker_constructor =\
        otu_picking_method_constructors[otu_picking_method]

    # split the input filepath into components used for generating
    # the output file name
    input_seqs_filepath = abspath(opts.input_seqs_filepath)
    input_seqs_dir, input_seqs_filename = split(input_seqs_filepath)
    input_seqs_basename, ext = splitext(input_seqs_filename)

    # create the output directory name (if not provided) and
    # create it if it doesn't already exist
    output_dir = opts.output_dir or otu_picking_method + '_picked_otus'
    create_dir(output_dir, fail_on_exist=False)

    # Create the output and log file names
    result_path = '%s/%s_otus.txt' % (output_dir, input_seqs_basename)
    log_path = '%s/%s_otus.log' % (output_dir, input_seqs_basename)
    failure_path = '%s/%s_failures.txt' % (output_dir, input_seqs_basename)

    # Perform OTU picking -- parameters and calls are made
    # on a per-method basis

    # cd-hit
    if otu_picking_method == 'cdhit':
        params = {'Similarity': similarity,
                  '-M': opts.max_cdhit_memory}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path,
                   prefix_prefilter_length=prefix_prefilter_length,
                   trie_prefilter=trie_prefilter)

    # uclust (de novo)
    elif otu_picking_method == 'uclust':
        params = {'Similarity': similarity,
                  'enable_rev_strand_matching': opts.enable_rev_strand_match,
                  'optimal': opts.optimal_uclust,
                  'exact': opts.exact_uclust,
                  # suppress_sort=True when seqs are or will be pre-sorted
                  'suppress_sort': user_sort,
                  'presort_by_abundance':
                  not suppress_presort_by_abundance_uclust,
                  'max_accepts': max_accepts,
                  'max_rejects': max_rejects,
                  'stepwords': stepwords,
                  'word_length': word_length,
                  'new_cluster_identifier': opts.denovo_otu_id_prefix,
                  'stable_sort': uclust_stable_sort,
                  'save_uc_files': save_uc_files,
                  'output_dir': output_dir,
                  'prefilter_identical_sequences': prefilter_identical_sequences}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path, HALT_EXEC=False)

    # usearch (usearch_qf)
    elif otu_picking_method == 'usearch':
        params = {'percent_id': similarity,
                  'maxrejects': max_rejects,
                  'w': word_length,
                  'save_intermediate_files': save_uc_files,
                  'output_dir': output_dir,
                  'percent_id_err': percent_id_err,
                  'minsize': minsize,
                  'abundance_skew': abundance_skew,
                  'db_filepath': db_filepath,
                  'perc_id_blast': perc_id_blast,
                  'de_novo_chimera_detection': de_novo_chimera_detection,
                  'reference_chimera_detection': reference_chimera_detection,
                  'cluster_size_filtering': cluster_size_filtering,
                  'remove_usearch_logs': remove_usearch_logs,
                  'derep_fullseq': derep_fullseq,
                  'chimeras_retention': chimeras_retention,
                  'verbose': verbose,
                  'minlen': minlen,
                  'rev': enable_rev_strand_match}

        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath, result_path=result_path,
                   log_path=log_path, HALT_EXEC=False)

    # usearch (usearch_qf) with reference OTU picking
    elif otu_picking_method == 'usearch_ref':
        params = {'percent_id': similarity,
                  'maxrejects': max_rejects,
                  'w': word_length,
                  'save_intermediate_files': save_uc_files,
                  'output_dir': output_dir,
                  'percent_id_err': percent_id_err,
                  'minsize': minsize,
                  'abundance_skew': abundance_skew,
                  'db_filepath': db_filepath,
                  'perc_id_blast': perc_id_blast,
                  'de_novo_chimera_detection': de_novo_chimera_detection,
                  'reference_chimera_detection': reference_chimera_detection,
                  'cluster_size_filtering': cluster_size_filtering,
                  'remove_usearch_logs': remove_usearch_logs,
                  'suppress_new_clusters': opts.suppress_new_clusters,
                  'derep_fullseq': derep_fullseq,
                  'chimeras_retention': chimeras_retention,
                  'verbose': verbose,
                  'minlen': minlen,
                  'rev': enable_rev_strand_match}

        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath, result_path=result_path,
                   refseqs_fp=refseqs_fp, failure_path=failure_path,
                   log_path=log_path, HALT_EXEC=False)

    # usearch 6.1 (de novo OTU picking only)
    elif otu_picking_method == 'usearch61':
        otu_prefix = opts.denovo_otu_id_prefix or 'denovo'
        params = {
            'percent_id': similarity,
            'wordlength': word_length,
            'save_intermediate_files': save_uc_files,
            'output_dir': output_dir,
            'remove_usearch_logs': remove_usearch_logs,
            'verbose': verbose,
            'minlen': minlen,
            'rev': enable_rev_strand_match,
            'usearch_fast_cluster': usearch_fast_cluster,
            'usearch61_sort_method': usearch61_sort_method,
            'usearch61_maxrejects': max_rejects,
            'usearch61_maxaccepts': max_accepts,
            'sizeorder': sizeorder,
            'threads': threads
        }

        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath, result_path=result_path,
                   log_path=log_path, otu_prefix=otu_prefix,
                   HALT_EXEC=False)

    # usearch 6.1 reference OTU picking
    elif otu_picking_method == 'usearch61_ref':
        otu_prefix = opts.denovo_otu_id_prefix or 'denovo'
        params = {
            'percent_id': similarity,
            'wordlength': word_length,
            'save_intermediate_files': save_uc_files,
            'output_dir': output_dir,
            'remove_usearch_logs': remove_usearch_logs,
            'verbose': verbose,
            'minlen': minlen,
            'rev': enable_rev_strand_match,
            'usearch_fast_cluster': usearch_fast_cluster,
            'usearch61_sort_method': usearch61_sort_method,
            'usearch61_maxrejects': max_rejects,
            'usearch61_maxaccepts': max_accepts,
            'sizeorder': sizeorder,
            'suppress_new_clusters': opts.suppress_new_clusters,
            'threads': threads
        }

        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath, refseqs_fp, result_path=result_path,
                   log_path=log_path, failure_path=failure_path,
                   otu_prefix=otu_prefix, HALT_EXEC=False)

    # uclust (reference-based)
    elif otu_picking_method == 'uclust_ref':
        params = {'Similarity': similarity,
                  'enable_rev_strand_matching': opts.enable_rev_strand_match,
                  'optimal': opts.optimal_uclust,
                  'exact': opts.exact_uclust,
                  # suppress_sort=True when seqs are or will be pre-sorted
                  'suppress_sort': user_sort,
                  'presort_by_abundance':
                  not suppress_presort_by_abundance_uclust,
                  'suppress_new_clusters': opts.suppress_new_clusters,
                  'max_accepts': max_accepts,
                  'max_rejects': max_rejects,
                  'stepwords': stepwords,
                  'word_length': word_length,
                  'new_cluster_identifier': opts.denovo_otu_id_prefix,
                  'stable_sort': uclust_stable_sort,
                  'save_uc_files': save_uc_files,
                  'output_dir': output_dir,
                  'prefilter_identical_sequences':
                  prefilter_identical_sequences,
                  'chimeras_retention': chimeras_retention}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath, refseqs_fp,
                   result_path=result_path, log_path=log_path,
                   failure_path=failure_path)

    # prefix/suffix
    elif otu_picking_method == 'prefix_suffix':
        otu_picker = otu_picker_constructor({})
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path,
                   prefix_length=prefix_length, suffix_length=suffix_length)

    # mothur
    elif otu_picking_method == 'mothur':
        params = {'Similarity': similarity,
                  'Algorithm': opts.clustering_algorithm}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path)

    # trie
    elif otu_picking_method == 'trie':
        params = {'Reverse': trie_reverse_seqs}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path)

    # blast
    elif otu_picking_method == 'blast':
        params = {'max_e_value': opts.max_e_value_blast,
                  'Similarity': similarity,
                  'min_aligned_percent': min_aligned_percent}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path,
                   blast_db=blast_db, refseqs_fp=refseqs_fp)

    # sortmerna
    elif otu_picking_method == 'sortmerna':
        params = {'max_e_value': sortmerna_e_value,
                  'similarity': similarity,
                  'coverage': sortmerna_coverage,
                  'threads': threads,
                  'blast': sortmerna_tabular,
                  'best': sortmerna_best_N_alignments,
                  'max_pos': sortmerna_max_pos,
                  'prefilter_identical_sequences':
                  prefilter_identical_sequences}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path,
                   sortmerna_db=sortmerna_db, refseqs_fp=refseqs_fp,
                   failure_path=failure_path)

    # sumaclust
    elif otu_picking_method == 'sumaclust':
        params = {'similarity': similarity,
                  'exact': sumaclust_exact,
                  'threads': threads,
                  'l': sumaclust_l,
                  'prefilter_identical_sequences':
                  prefilter_identical_sequences,
                  'denovo_otu_id_prefix': denovo_otu_id_prefix}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path)

    # swarm
    elif otu_picking_method == 'swarm':
        params = {'resolution': swarm_resolution,
                  'threads': threads,
                  'denovo_otu_id_prefix': denovo_otu_id_prefix}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path)

    # other -- shouldn't be able to get here as a KeyError would have
    # been raised earlier
    else:
        raise ValueError("Unknown OTU picking method: %s" % otu_picking_method)


if __name__ == "__main__":
    main()
