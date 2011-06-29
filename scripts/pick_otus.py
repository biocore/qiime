#!/usr/bin/env python
# File created on 09 Feb 2010

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Rob Knight","Greg Caporaso", "Kyle Bittinger",
               "Jens Reeder","William Walters"] 
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from os.path import splitext, split, exists
from os import makedirs
from qiime.util import make_option
from qiime.util import get_tmp_filename
from cogent.util.misc import remove_files
from qiime.util import (parse_command_line_parameters, create_dir)
from qiime.sort import sort_fasta_by_abundance
from qiime.pick_otus  import otu_picking_method_constructors,\
 otu_picking_method_choices, MothurOtuPicker

script_info={}
script_info['brief_description'] = """OTU picking"""
script_info['script_description'] = """The OTU picking step assigns similar sequences to operational taxonomic units, or OTUs, by clustering sequences based on a user-defined similarity threshold. Sequences which are similar at or above the threshold level are taken to represent the presence of a taxonomic unit (e.g., a genus, when the similarity threshold is set at 0.94) in the sequence collection.

Currently, the following clustering methods have been implemented in QIIME:

1. cd-hit (Li & Godzik, 2006; Li, Jaroszewski, & Godzik, 2001), which applies a \"longest-sequence-first list removal algorithm\" to cluster sequences.  

2. blast (Altschul, Gish, Miller, Myers, & Lipman, 1990), which compares and clusters each sequence against a reference database of sequences.

3. Mothur (Schloss et al., 2009), which requires an input file of aligned sequences.  The input file of aligned sequences may be generated from an input file like the one described below by running align_seqs.py.  For the Mothur method, the clustering algorithm may be specified as nearest-neighbor, furthest-neighbor, or average-neighbor.  The default algorithm is furthest-neighbor.

4. prefix/suffix [Qiime team, unpublished], which will collapse sequences which are identical in their first and/or last bases (i.e., their prefix and/or suffix). The prefix and suffix lengths are provided by the user and default to 50 each.

5. Trie [Qiime team, unpublished], which collapsing identical sequences and sequences which are subsequences of other sequences.

6. uclust (Robert Edgar, unpublished, 2009), creates \"seeds\" of sequences which generate clusters based on percent identity.

The primary inputs for pick_otus.py are:

1. A FASTA file containing sequences to be clustered

2. An OTU threshold (default is 0.97, roughly corresponding to species-level OTUs);

3. The method to be applied for clustering sequences into OTUs.

pick_otus.py takes a standard fasta file as input.

"""
script_info['script_usage'] = []

script_info['script_usage'].append(("""Example (uclust method, default):""","""Using the seqs.fna file generated from split_libraries.py and outputting the results to the directory \"picked_otus/\", while using default parameters (0.97 sequence similarity, no reverse strand matching):""","""pick_otus.py -i seqs.fna -o picked_otus/"""))


script_info['script_usage'].append(("""""","""To change the percent identity to a lower value, such as 90%, and also enable reverse strand matching, the script would be the following:""","""pick_otus.py -i seqs.fna -o picked_otus/ -s 0.90 -z"""))

script_info['script_usage'].append(("""Uclust Reference-based OTU picking example""","""uclust_ref can be passed via -m to pick OTUs against a reference set where sequences within the similarity threshold to a reference sequence will cluster to an OTU defined by that reference sequence, and sequences outside of the similarity threshold to a reference sequence will form new clusters. OTU identifiers will be set to reference sequence identifiers when sequences cluster to reference sequences, and 'qiime_otu_<integer>' for new OTUs. Creation of new clusters can be suppressed by passing -C, in which case sequences outside of the similarity threshold to any reference sequence will be listed as failures in the log file, and not included in any OTU.""","""pick_otus.py -i seqs.fna -r core_set_unaligned.fasta_11_8_07 -m uclust_ref"""))

script_info['script_usage'].append(("""Example (cdhit method):""","""Using the seqs.fna file generated from split_libraries.py and outputting the results to the directory \"picked_otus/\", while using default parameters (0.97 sequence similarity, no prefix filtering):""","""pick_otus.py -i seqs.fna -m cdhit -o picked_otus/"""))

script_info['script_usage'].append(("""""","""Currently the cd-hit OTU picker allows for users to perform a pre-filtering step, so that highly similar sequences are clustered prior to OTU picking. This works by collapsing sequences which begin with an identical n-base prefix, where n is specified by the -n parameter. A commonly used value here is 100 (e.g., -n 100). So, if using this filter with -n 100, all sequences which are identical in their first 100 bases will be clustered together, and only one representative sequence from each cluster will be passed to cd-hit. This is used to greatly increase the run-time of cd-hit-based OTU picking when working with very large sequence collections, as shown by the following command:""","""pick_otus.py -i seqs.fna -m cdhit -o picked_otus/ -n 100"""))

script_info['script_usage'].append(("""""","""Alternatively, if the user would like to collapse identical sequences, or those which are subsequences of other sequences prior to OTU picking, they can use the trie prefiltering (\"-t\") option as shown by the following command:""","""pick_otus.py -i seqs.fna -m cdhit -o picked_otus/ -t"""))

script_info['script_usage'].append(("""""","""Note: It is highly recommended to use one of the prefiltering methods when analyzing large dataset (>100,000 seqs) to reduce run-time.""",""""""))


script_info['script_usage'].append(("""BLAST OTU-Picking Example:""","""OTUs can be picked against a reference database using the BLAST OTU picker. This is useful, for example, when different regions of the SSU RNA have sequenced and a sequence similarity based approach like cd-hit therefore wouldn't work. When using the BLAST OTU picking method, the user must supply either a reference set of sequences or a reference database to compare against. The OTU identifiers resulting from this step will be the sequence identifiers in the reference database. This allows for use of a pre-existing tree in downstream analyses, which again is useful in cases where different regions of the 16s gene have been sequenced.

The following command can be used to blast against a reference sequence set, using the default E-value and sequence similarity (0.97) parameters:""","""pick_otus.py -i seqs.fna -o picked_otus/ -m blast -r ref_seq_set.fna"""))

script_info['script_usage'].append(("""""","""If you already have a pre-built BLAST database, you can pass the database prefix as shown by the following command:""","""pick_otus.py -i seqs.fna -o picked_otus/ -m blast -b ref_database"""))

script_info['script_usage'].append(("""""","""If the user would like to change the sequence similarity (\"-s\") and/or the E-value (\"-e\") for the blast method, they can use the following command:""","""pick_otus.py -i seqs.fna -o picked_otus/ -m blast -s 0.90 -e 1e-30"""))

script_info['script_usage'].append(("""Prefix-suffix OTU Picking Example:""","""OTUs can be picked by collapsing sequences which being and/or end with identical bases (i.e., identical prefixes or suffixes). This OTU picker is currently likely to be of limited use on its own, but will be very useful in collapsing very similar sequences in a chained OTU picking strategy that is currently in development. For example, user will be able to pick OTUs with this method, followed by representative set picking, and then re-pick OTUs on their representative set. This will allow for highly similar sequences to be collapsed, followed by running a slower OTU picker. This ability to chain OTU pickers is not yet supported in QIIME. The following command illustrates how to pick OTUs by collapsing sequences which are identical in their first 50 and last 25 bases:""","""pick_otus.py -i seqs.fna -o picked_otus/ -m prefix_suffix -p 50 -u 25"""))

script_info['script_usage'].append(("""Mothur OTU Picking Example:""","""The Mothur program (http://www.mothur.org/) provides three clustering algorithms for OTU formation: furthest-neighbor (complete linkage), average-neighbor (group average), and nearest-neighbor (single linkage). Details on the algorithms may be found on the Mothur website and publications (Schloss et al., 2009). However, the running times of Mothur's clustering algorithms scale with the number of sequences squared, so the program may not be feasible for large data sets.

The following command may be used to create OTU's based on a furthest-neighbor algorithm (the default setting):""","""pick_otus.py -i seqs.fna -o picked_otus/ -m mothur"""))

script_info['script_usage'].append(("""""","""If you prefer to use a nearest-neighbor algorithm instead, you may specify this with the '-c' flag:""","""pick_otus.py -i seqs.fna -o picked_otus/ -m mothur -c nearest"""))

script_info['script_usage'].append(("""""","""The sequence similarity parameter may also be specified. For example, the following command may be used to create OTU's at the level of 95% similarity:""","""pick_otus.py -i seqs.fna -o picked_otus/ -m mothur -s 0.90"""))

script_info['output_description'] = """The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an (arbitrary) cluster identifier, and the remaining fields correspond to sequence identifiers assigned to that cluster. Sequence identifiers correspond to those provided in the input FASTA file.

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
    make_option('-i', '--input_seqs_filepath',
        help='Path to input sequences file'),
    ]

script_info['optional_options'] = [
    make_option('-m', '--otu_picking_method', type='choice',
        choices=otu_picking_method_choices, default = "uclust",
        help=('Method for picking OTUs.  Valid choices are: ' +\
              ', '.join(otu_picking_method_choices) +\
              '. The mothur method requires an input file ' +\
              'of aligned sequences [default: %default]')),
    make_option('-c', '--clustering_algorithm', type='choice',
        choices=MothurOtuPicker.ClusteringAlgorithms, default='furthest',
        help=('Clustering algorithm for mothur otu picking method.  Valid ' +\
              'choices are: ' +\
              ', '.join(MothurOtuPicker.ClusteringAlgorithms) +\
              '. [default: %default]')),
    make_option('-M', '--max_cdhit_memory', type=int, default=400,
        help=('Maximum available memory to cd-hit-est (via the program\'s -M '
              'option) for cdhit OTU picking method (units of Mbyte) '
              '[default: %default]')),
    make_option('-o', '--output_dir',\
        help=('Path to store result file '
              '[default: ./<OTU_METHOD>_picked_otus/]')),
    make_option('-r', '--refseqs_fp',
        help=('Path to reference sequences to search against when using -m '
              'blast or -m uclust_ref [default: %default]')),
    make_option('-b', '--blast_db',
        help=('Pre-existing database to blast against when using -m blast '
              '[default: %default]')),
    make_option('--min_aligned_percent',
        help=('Minimum percent of query sequence that can be aligned to consider a hit '
              ' (BLAST OTU picker only) [default: %default]'),default=0.50,type='float'),
    make_option('-s', '--similarity', type='float', default=0.97,
        help=('Sequence similarity threshold (for cdhit, uclust, or uclust_ref) '
              '[default: %default]')),
    make_option('-e', '--max_e_value', type='float', default=1e-10,
        help=('Max E-value when clustering with BLAST [default: %default]')),
    make_option('-q', '--trie_reverse_seqs', action='store_true',
        default=False,
        help=('Reverse seqs before picking OTUs with the Trie OTU picker for '
              'suffix (rather than prefix) collapsing [default: %default]')),
    make_option('-n', '--prefix_prefilter_length', type=int, default=None,
        help=('Prefilter data so seqs with identical first '
              'prefix_prefilter_length are automatically grouped into a single '
              'OTU.  This is useful for large sequence collections where OTU '
              'picking doesn\'t scale well [default: %default; 100 is a good '
              'value]')),
    make_option('-t', '--trie_prefilter', action='store_true',
        default=False,
        help=('prefilter data so seqs which are identical prefixes of a longer '
              'seq are automatically grouped into a single OTU; useful for '
              'large sequence collections where OTU picking doesn\'t scale '
              'well [default: %default]')),
    make_option('-p', '--prefix_length', type=int, default=50,
        help=('Prefix length when using the prefix_suffix otu picker; '
              'WARNING: CURRENTLY DIFFERENT FROM prefix_prefilter_length (-n)! '
              '[default: %default]')),
    make_option('-u', '--suffix_length', type=int, default=50,
        help=('Suffix length when using the prefix_suffix otu picker '
              '[default: %default]')),
    make_option('-z', '--enable_rev_strand_match', action='store_true',
        default=False,
        help=('Enable reverse strand matching for uclust otu picking, '
              'will double the amount of memory used. [default: %default]')),
    make_option('-D','--suppress_presort_by_abundance_uclust', action='store_true', 
              default=False,
              help=('Suppress presorting of sequences by abundance when picking'
              ' OTUs with uclust or uclust_ref [default: %default]')),
    make_option('-A','--optimal_uclust', action='store_true', 
              default=False,
              help=('Pass the --optimal flag to uclust for uclust otu'
              ' picking. [default: %default]')),
    make_option('-E','--exact_uclust', action='store_true', 
              default=False,
              help=('Pass the --exact flag to uclust for uclust otu'
              ' picking. [default: %default]')),
    make_option('-B','--user_sort', action='store_true', 
              default=False,
              help=('Pass the --user_sort flag to uclust for uclust otu'
              ' picking. [default: %default]')),
    make_option('-C','--suppress_new_clusters',action='store_true',default=False,
              help="Suppress creation of new clusters using seqs that don't" +
              " match reference when using -m uclust_ref [default: %default]"),
    make_option('--max_accepts',type='int',default=20,
              help="max_accepts value to uclust and "
                   "uclust_ref [default: %default]"),
    make_option('--max_rejects',type='int',default=500,
              help="max_rejects value to uclust and "
                   "uclust_ref [default: %default]"),
   make_option('--stepwords',type='int',default=20,
             help="stepwords value to uclust and "
                  "uclust_ref [default: %default]"),
   make_option('--word_length',type='int',default=12,
             help="w value to uclust and "
                  "uclust_ref [default: %default]"),
    make_option('--uclust_otu_id_prefix',default=None,
              help=("OTU identifier prefix (string) for the de novo uclust" 
                    " OTU picker [default: %default, OTU ids are ascending"
                    " integers]")),
    make_option('--uclust_stable_sort',default=True,action='store_true',
              help=("Deprecated: stable sort enabled by default, pass "
                    "--uclust_suppress_stable_sort to disable [default: %default]")),
    make_option('--suppress_uclust_stable_sort',default=False,action='store_true',
        help=("Don't pass --stable-sort to uclust [default: %default]")),
    make_option('--suppress_uclust_prefilter_exact_match',
                default=False,action='store_true',
        help=("Don't collapse exact matches before calling uclust [default: %default]")),
    make_option('-d', '--save_uc_files', default=True, action='store_false',
              help=("Enable preservation of intermediate uclust (.uc) files "
              "that are used to generate clusters via uclust. "
              "[default: %default]"))
    ]

script_info['version'] = __version__

def main():
    # Parse the command line parameters
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)
    
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
    prefilter_identical_sequences = not opts.suppress_uclust_prefilter_exact_match
    

    
    # Input validation to throw a useful error message on common mistakes
    if (otu_picking_method == 'cdhit' and
        similarity < 0.80):
        option_parser.error('cdhit requires similarity >= 0.80.')
    
    if (otu_picking_method == 'blast' and
        refseqs_fp == None and
        blast_db == None):
           option_parser.error('blast requires refseqs_fp or blast_db')
    
    if (otu_picking_method == 'uclust_ref'):
        if (refseqs_fp == None):
            option_parser.error('uclust_ref requires refseqs_fp')
        elif not exists(refseqs_fp):
            option_parser.error('refseqs_fp %s does not exist' %refseqs_fp)
    # End input validation
    
    
    # use the otu_picking_method value to get the otu picker constructor
    otu_picker_constructor =\
     otu_picking_method_constructors[otu_picking_method]
    
    # split the input filepath into components used for generating
    # the output file name
    input_seqs_filepath = opts.input_seqs_filepath
    input_seqs_dir, input_seqs_filename = split(input_seqs_filepath)
    input_seqs_basename, ext = splitext(input_seqs_filename)
    
    # create the output directory name (if not provided) and 
    # create it if it doesn't already exist
    output_dir = opts.output_dir or otu_picking_method + '_picked_otus'
    create_dir(output_dir, fail_on_exist=False)
    
    # Create the output and log file names
    result_path = '%s/%s_otus.txt' % (output_dir,input_seqs_basename)
    log_path = '%s/%s_otus.log' % (output_dir,input_seqs_basename)
    failure_path = '%s/%s_failures.txt' % (output_dir,input_seqs_basename)
    
    ## Perform OTU picking -- parameters and calls are made 
    ## on a per-method basis
    
    ## cd-hit
    if otu_picking_method == 'cdhit':
        params = {'Similarity':opts.similarity,
                  '-M':opts.max_cdhit_memory}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path,log_path=log_path,
                   prefix_prefilter_length=prefix_prefilter_length,
                   trie_prefilter=trie_prefilter)
    
    ## uclust (de novo)
    elif otu_picking_method == 'uclust':
        params = {'Similarity':opts.similarity,
        'enable_rev_strand_matching':opts.enable_rev_strand_match,
        'optimal':opts.optimal_uclust,
        'exact':opts.exact_uclust,
        # suppress_sort=True when seqs are or will be pre-sorted
        'suppress_sort':user_sort,
        'presort_by_abundance': not suppress_presort_by_abundance_uclust,
        'max_accepts':max_accepts,
        'max_rejects':max_rejects,
        'stepwords':stepwords,
        'word_length':word_length,
        'new_cluster_identifier':opts.uclust_otu_id_prefix,
        'stable_sort':uclust_stable_sort,
        'save_uc_files':save_uc_files,
        'output_dir':output_dir,
        'prefilter_identical_sequences':prefilter_identical_sequences}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path,log_path=log_path,HALT_EXEC=False)
             
    ## uclust (reference-based)
    elif otu_picking_method == 'uclust_ref':
        params = {'Similarity':opts.similarity,
        'enable_rev_strand_matching':opts.enable_rev_strand_match,
        'optimal':opts.optimal_uclust,
        'exact':opts.exact_uclust,
        # suppress_sort=True when seqs are or will be pre-sorted
        'suppress_sort':user_sort,
        'presort_by_abundance': not suppress_presort_by_abundance_uclust,
        'suppress_new_clusters':opts.suppress_new_clusters,
        'max_accepts':max_accepts,
        'max_rejects':max_rejects,
        'stepwords':stepwords,
        'word_length':word_length,
        'new_cluster_identifier':opts.uclust_otu_id_prefix,
        'stable_sort':uclust_stable_sort,
        'save_uc_files':save_uc_files,
        'output_dir':output_dir,
        'prefilter_identical_sequences':prefilter_identical_sequences}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,refseqs_fp,
                   result_path=result_path,log_path=log_path,
                   failure_path=failure_path)
         
    ## prefix/suffix
    elif otu_picking_method == 'prefix_suffix':
        otu_picker = otu_picker_constructor({})
        otu_picker(input_seqs_filepath,
                   result_path=result_path,log_path=log_path,
                   prefix_length=prefix_length,suffix_length=suffix_length)
         
    ## mothur
    elif otu_picking_method == 'mothur':
        params = {'Similarity': opts.similarity,
                  'Algorithm': opts.clustering_algorithm}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path)
                   
    ## trie
    elif otu_picking_method == 'trie':
        params = {'Reverse':trie_reverse_seqs}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path)
                   
    ## blast
    elif otu_picking_method == 'blast':
        params = {'max_e_value':opts.max_e_value,
                  'Similarity': opts.similarity,
                  'min_aligned_percent':min_aligned_percent}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path,
                   blast_db=opts.blast_db,refseqs_fp=opts.refseqs_fp)
    
    ## other -- shouldn't be able to get here as a KeyError would have
    ## been raised earlier
    else:
        raise ValueError, "Unknown OTU picking method: %s" % otu_picking_method

if __name__ == "__main__":
    main()

