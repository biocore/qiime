#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger", "Antonio Gonzalez Pena", "David Soergel"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option,
                        get_rdp_jarpath,
                        load_qiime_config,
                        get_tmp_filename,
                        remove_files)
from os import system, remove, path, mkdir
from os.path import split, splitext
from qiime.assign_taxonomy import (
    BlastTaxonAssigner, MothurTaxonAssigner, RdpTaxonAssigner,
    RtaxTaxonAssigner, Tax2TreeTaxonAssigner, validate_rdp_version,
    )

assignment_method_constructors = {
    'blast': BlastTaxonAssigner,
    'mothur': MothurTaxonAssigner,
    'rdp': RdpTaxonAssigner,
    'rtax': RtaxTaxonAssigner,
    'tax2tree': Tax2TreeTaxonAssigner
}

assignment_method_choices = ['rdp','blast','rtax','mothur', 'tax2tree']

options_lookup = get_options_lookup()

qiime_config = load_qiime_config()

script_info={}
script_info['brief_description']="""Assign taxonomy to each sequence"""
script_info['script_description']="""Contains code for assigning taxonomy, using several techniques.

Given a set of sequences, %prog attempts to assign the taxonomy of each sequence. Currently there are three methods implemented: assignment with BLAST, assignment with the RDP classifier, and assignment with the RTAX classifier. The output of this step is a mapping of input sequence identifiers (1st column of output file) to taxonomy (2nd column) and quality score (3rd column). The sequence identifier of the best BLAST hit is also included if the blast method is used (4th column).

Example reference data sets and id_to_taxonomy maps can be found in the Greengenes OTUs. To get the latest build of those click the "Most recent Greengenes OTUs" link on the top right of http://blog.qiime.org. After downloading and unzipping you can use the following following files as -r and -t. As of this writing the latest build was gg_otus_4feb2011, but that portion of path to these files will change with future builds. Modify these paths accordining when calling %prog.

-r gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta
-t gg_otus_4feb2011/taxonomies/greengenes_tax_rdp_train.txt (best for retraining the RDP classifier)
-t gg_otus_4feb2011/taxonomies/greengenes_tax.txt (best for BLAST taxonomy assignment)
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Sample Assignment with BLAST:""","""
Taxonomy assignments are made by searching input sequences against a blast database of pre-assigned reference sequences. If a satisfactory match is found, the reference assignment is given to the input sequence. This method does not take the hierarchical structure of the taxonomy into account, but it is very fast and flexible. If a file of reference sequences is provided, a temporary blast database is built on-the-fly. The quality scores assigned by the BLAST taxonomy assigner are e-values.

To assign the sequences to the representative sequence set, using a reference set of sequences and a taxonomy to id assignment text file, where the results are output to default directory "blast_assigned_taxonomy", you can run the following command:""","""%prog -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt"""))
script_info['script_usage'].append(("""""","""Optionally, the user could changed the E-value ("-e"), using the following command:""","""%prog -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt -e 0.01"""))
script_info['script_usage'].append(("""Assignment with the RDP Classifier:""","""The RDP Classifier program (Wang, Garrity, Tiedje, & Cole, 2007) assigns taxonomies by matching sequence segments of length 8 to a database of previously assigned sequences. It uses a naive bayesian algorithm, which means that for each potential assignment, it attempts to calculate the probability of the observed matches, assuming that the assignment is correct and that the sequence segments are completely independent. The RDP Classifier is distributed with a pre-built database of assigned sequence, which is used by default. The quality scores provided by the RDP classifier are confidence values.

Note: If a reference set of sequences and taxonomy to id assignment file are provided, the script will use them to generate a new training dataset for the RDP Classifier on-the-fly.  Because of the RDP Classifier's implementation, all lineages in the training dataset must contain the same number of ranks.

To assign the representative sequence set, where the output directory is "rdp_assigned_taxonomy", you can run the following command:""","""%prog -i repr_set_seqs.fasta -m rdp"""))
script_info['script_usage'].append(("""""","""Alternatively, the user could change the minimum confidence score ("-c"), using the following command:""","""%prog -i repr_set_seqs.fasta -m rdp -c 0.85"""))
script_info['script_usage'].append(("""Sample Assignment with RTAX:""","""
Taxonomy assignments are made by searching input sequences against a fasta database of pre-assigned reference sequences. All matches are collected which match the query within 0.5% identity of the best match.  A taxonomy assignment is made to the lowest rank at which more than half of these hits agree.  Note that both unclustered read fasta files are required as inputs in addition to the representative sequence file.

To make taxonomic classifications of the representative sequences, using a reference set of sequences and a taxonomy to id assignment text file, where the results are output to default directory "rtax_assigned_taxonomy", you can run the following command:""","""%prog -i rtax_repr_set_seqs.fasta -m rtax --read_1_seqs_fp read_1.seqs.fna --read_2_seqs_fp read_2.seqs.fna -r rtax_ref_seq_set.fna -t rtax_id_to_taxonomy.txt"""))
script_info['script_usage'].append(("""Sample Assignment with Mothur:""", """The Mothur software provides a naive bayes classifier similar to the RDP Classifier.  A set of training sequences and id-to-taxonomy assignments must be provided.  Unlike the RDP Classifier, sequences in the training set may be assigned at any level of the taxonomy.

To make taxonomic classifications of the representative sequences, where the results are output to default directory \"mothur_assigned_taxonomy\", you can run the following command:""", "%prog -i mothur_repr_set_seqs.fasta -m mothur -r mothur_ref_seq_set.fna -t mothur_id_to_taxonomy.txt"))
script_info['output_description']="""The consensus taxonomy assignment implemented here is the most detailed lineage description shared by 90% or more of the sequences within the OTU (this level of agreement can be adjusted by the user). The full lineage information for each sequence is one of the output files of the analysis. In addition, a conflict file records cases in which a phylum-level taxonomy assignment disagreement exists within an OTU (such instances are rare and can reflect sequence misclassification within the greengenes database)."""
script_info['required_options']=[\
   options_lookup['fasta_as_primary_input']\
]

default_reference_seqs_fp = qiime_config['assign_taxonomy_reference_seqs_fp']
default_id_to_taxonomy_fp = qiime_config['assign_taxonomy_id_to_taxonomy_fp']

script_info['optional_options']=[\
 make_option('-t', '--id_to_taxonomy_fp',type="existing_filepath",
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a semicolon-separated'
         ' list. For assignment with rdp, each assigned taxonomy must be '
         'exactly 6 levels deep. [default: %s; REQUIRED when method is blast]'\
          % default_id_to_taxonomy_fp,
          default=default_id_to_taxonomy_fp),\
 make_option('-r', '--reference_seqs_fp',type="existing_filepath",
        help='Path to reference sequences.  For assignment with blast, these '
        'are used to generate a blast database. For assignment with rdp, they '
        'are used as training sequences for the classifier. '
        '[default: %s; REQUIRED if -b is not provided when method is blast]'\
        % default_reference_seqs_fp,
        default=default_reference_seqs_fp),\
 make_option('-p', '--training_data_properties_fp',type="existing_filepath",
        help='Path to ".properties" file in pre-compiled training data for the '
        'RDP Classifier.  This option is overridden by the -t and -r options. '
        '[default: %default]'),\
 make_option('--read_1_seqs_fp',type="existing_filepath",
        help='Path to fasta file containing the first read from paired-end '
        'sequencing, prior to OTU clustering (used for RTAX only). '
        '[default: %default]'),\
 make_option('--read_2_seqs_fp',type="existing_filepath",
        help='Path to fasta file containing a second read from paired-end '
        'sequencing, prior to OTU clustering (used for RTAX only). '
        '[default: %default]'),\
 make_option('--single_ok',action="store_true",
        help='When classifying paired ends, allow fallback to single-ended '
        'classification when the mate pair is lacking (used for RTAX only). '
        '[default: %default]',default=False),\
 make_option('--no_single_ok_generic',action="store_true",
        help='When classifying paired ends, do not allow fallback to single-ended '
        'classification when the mate pair is overly generic (used for RTAX only). '
        '[default: %default]',default=False),\
 make_option('--read_id_regex',type="string",
        help='Used to parse the result of OTU clustering, to get the read_1_id '
        'for each clusterID.  (used for RTAX only). '
        '[default: %default]',default="\\S+\\s+(\\S+)"),\
 make_option('--amplicon_id_regex',type="string",
        help='Used to parse the result of split_libraries, to get the ampliconID '
        'for each read_1_id.  Two groups capture read_1_id and ampliconID, '
        'respectively.  (used for RTAX only). '
        '[default: %default]',default="(\\S+)\\s+(\\S+?)\/"),\
 make_option('--header_id_regex',type="string",
        help='Used to choose the part of the header in the OTU clustering file '
        'that Rtax reports back as the ID.  The default uses the amplicon ID, '
        'not including /1 or /3, as the primary key for the query sequences. '
        '(used for RTAX only). '
        '[default: %default]',default="\\S+\\s+(\\S+?)\/"),\
 make_option('-m', '--assignment_method', type='choice',
        help='Taxon assignment method, either blast, mothur, rdp, or rtax '
        '[default:%default]',
        choices=assignment_method_choices, default="rdp"),\
 make_option('-b', '--blast_db', type='string',
        help='Database to blast against.  Must provide either --blast_db or '
        '--reference_seqs_db for assignment with blast [default: %default]'),\
 make_option('-c', '--confidence', type='float',
        help='Minimum confidence to record an assignment, only used for rdp '
        'and mothur methods [default: %default]', default=0.80),\
 make_option('--rdp_max_memory', default=1500, type='int',
        help='Maximum memory allocation, in MB, for Java virtual machine when '
        'using the rdp method.  Increase for large training sets [default: %default]'),\
 make_option('-e', '--e_value', type='float',
        help='Maximum e-value to record an assignment, only used for blast '
        'method [default: %default]',default=0.001),\
 make_option('--tree_fp', type='existing_filepath',
        help='The filepath to a prebuilt tree containing both the representative '
        'and reference sequences. Required for Tax2Tree assignment.'),\
 make_option('-o','--output_dir', type='new_dirpath',\
          help='Path to store result file '+\
          '[default: <ASSIGNMENT_METHOD>_assigned_taxonomy]')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    assignment_method = opts.assignment_method

    if assignment_method == 'blast':
        if not opts.id_to_taxonomy_fp:
            option_parser.error('Option --id_to_taxonomy_fp is required when '
                         'assigning with blast.')
        if not (opts.reference_seqs_fp or opts.blast_db):
            option_parser.error('Either a blast db (via -b) or a collection of '
                         'reference sequences (via -r) must be passed to '
                         'assign taxonomy using blast.')

    if assignment_method == 'rdp':
        try:
            validate_rdp_version()
        except RuntimeError, e:
            option_parser.error(e)

        if opts.id_to_taxonomy_fp is not None:
            if opts.reference_seqs_fp is None:
                option_parser.error(
                    'A filepath for reference sequences must be '
                    'specified (via -r) along with the id_to_taxonomy '
                    'file to train the Rdp Classifier.')
        elif opts.reference_seqs_fp is not None:
                option_parser.error(
                    'A filepath for an id to taxonomy map must be '
                    'specified (via -t) along with the reference '
                    'sequences fp to train the Rdp Classifier.')
        else:
            pass

    if assignment_method == 'rtax':
        if opts.id_to_taxonomy_fp is None or opts.reference_seqs_fp is None:
            option_parser.error('RTAX classification requires both a filepath for '
                         'reference sequences (via -r) and an id_to_taxonomy '
                         'file (via -t).')
        if opts.read_1_seqs_fp is None: # or opts.read_2_seqs_fp is None:
            option_parser.error('RTAX classification requires the FASTA files '
                         'produced by split_illumina_fastq.py for both reads, '
                         'in addition to the cluster representatives.  Pass '
                         'these via --read_1_seqs_fp and --read_2_seqs_fp.')

    if assignment_method == 'mothur':
        if None in [opts.id_to_taxonomy_fp, opts.reference_seqs_fp]:
            option_parser.error(
                'Mothur classification requires both a filepath for '
                'reference sequences (via -r) and an id_to_taxonomy '
                'file (via -t).')

    if assignment_method == 'tax2tree':
        if opts.tree_fp is None:
            option_parser.error('Tax2Tree classification requires a '
            'filepath to a prebuilt tree (via --tree_fp) containing '
            'both the representative and reference sequences. Check '
            'Tax2Tree documentation for help building a tree.')
        if opts.id_to_taxonomy_fp is None:
            option_parser.error('Tax2Tree classification requires a '
            'filepath for an id_to_taxonomy file (via -t).')

    taxon_assigner_constructor =\
     assignment_method_constructors[assignment_method]
    input_sequences_filepath = opts.input_fasta_fp

    try:
        id_to_taxonomy_fp = opts.id_to_taxonomy_fp
        params = {'id_to_taxonomy_filepath':id_to_taxonomy_fp}
    except IndexError:
        params = {}

    # Build the output filenames
    output_dir = opts.output_dir or assignment_method + '_assigned_taxonomy'
    try:
        mkdir(output_dir)
    except OSError:
        # output_dir already exists
        pass

    fpath, ext = splitext(input_sequences_filepath)
    input_dir, fname = split(fpath)
    result_path = output_dir + '/' + fname + '_tax_assignments.txt'
    log_path = output_dir + '/' + fname + '_tax_assignments.log'

    if assignment_method == 'blast':
        # one of these must have a value, otherwise we'd have
        # an optparse error
        if opts.blast_db:
            params['blast_db'] = opts.blast_db
        else:
            params['reference_seqs_filepath'] = opts.reference_seqs_fp
        params['Max E value'] = opts.e_value

    elif assignment_method == 'mothur':
        params['Confidence'] = opts.confidence
        params['id_to_taxonomy_fp'] = opts.id_to_taxonomy_fp
        params['reference_sequences_fp'] = opts.reference_seqs_fp

    elif assignment_method == 'rdp':
        params['Confidence'] = opts.confidence
        params['id_to_taxonomy_fp'] = opts.id_to_taxonomy_fp
        params['reference_sequences_fp'] = opts.reference_seqs_fp
        params['training_data_properties_fp'] = opts.training_data_properties_fp
        params['max_memory'] = "%sM" % opts.rdp_max_memory

    elif assignment_method == 'rtax':
       params['id_to_taxonomy_fp'] = opts.id_to_taxonomy_fp
       params['reference_sequences_fp'] = opts.reference_seqs_fp
       params['read_1_seqs_fp'] = opts.read_1_seqs_fp
       params['read_2_seqs_fp'] = opts.read_2_seqs_fp
       params['single_ok'] = opts.single_ok
       params['no_single_ok_generic'] = opts.no_single_ok_generic
       params['header_id_regex'] = opts.header_id_regex
       params['read_id_regex'] = opts.read_id_regex
       params['amplicon_id_regex'] = opts.amplicon_id_regex

    elif assignment_method == 'tax2tree':
        params['id_to_taxonomy_fp'] = opts.id_to_taxonomy_fp
        params['tree_fp'] = opts.tree_fp

    else:
        # should not be able to get here as an unknown classifier would
        # have raised an optparse error
        exit(1)
    temp_result_path = get_tmp_filename(prefix='assign-tax')
    taxon_assigner = taxon_assigner_constructor(params)
    taxon_assigner(input_sequences_filepath,\
     result_path=temp_result_path,log_path=log_path)
    
    ## This is an ugly hack, and needs to be pushed upstream to
    ## the taxon assigners. The output taxonomy maps that are returned by the 
    ## taxon assigners contain the full sequence headers as the first field
    ## (so including "comment" text in the fasta headers), but for consistency
    ## with the input taxonomy maps, should only contain the sequence identifier.
    ## This modifies those entries to contain only the sequence identifer, 
    ## discarding any comment information. The formatting of these result files
    ## needs to be centralized, and at that stage this processing should
    ## happen there rather than here.
    result_f = open(result_path,'w')
    for line in open(temp_result_path,'U'):
        fields = line.strip().split('\t')
        seq_id = fields[0].split()[0]
        result_f.write('%s\t%s\n' % (seq_id,'\t'.join(fields[1:])))
    result_f.close()
    remove_files([temp_result_path])
        


if __name__ == "__main__":
    main()
