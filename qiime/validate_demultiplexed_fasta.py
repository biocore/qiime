#!/usr/bin/env python
# File created Feb 1 2012

from __future__ import division

__author__ = "William Anton Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Anton Walters"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "William Anton Walters"
__email__ = "william.a.walters@gmail.com"
__status__ = "Development"

from collections import defaultdict
from os.path import split, join

from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.record import RecordError
from cogent.parse.tree import DndParser

from qiime.check_id_map import process_id_map
from qiime.split_libraries import expand_degeneracies


def get_mapping_details(mapping_fp):
    """ Returns SampleIDs, Barcodes, Primer seqs from mapping file
    
    mapping_fp: filepath to mapping file
    """
    
    mapping_f = open(mapping_fp, "U")
    
    # Only using the id_map and the errors from parsing the mapping file.
    headers, id_map, desc, run_description, errors, warnings = \
        process_id_map(mapping_f)
        
    mapping_f.close()
    
    # Errors means problems with SampleIDs or headers
    if errors:
        raise ValueError,('Error in mapping file, please validate '+\
         'mapping file with check_id_map.py')
         
    sample_ids = id_map.keys()
    
    barcode_seqs = []
    raw_linkerprimer_seqs = []
    
    for curr_id in id_map:
        barcode_seqs.append(id_map[curr_id]['BarcodeSequence'])
        raw_linkerprimer_seqs.append(id_map[curr_id]['LinkerPrimerSequence'])
    
    # remove duplicates    
    raw_linkerprimer_seqs = set(raw_linkerprimer_seqs)
    
    linker_primer_seqs = expand_degeneracies(raw_linkerprimer_seqs)
    
    return set(sample_ids), set(barcode_seqs), set(linker_primer_seqs)
    
def verify_valid_fasta_format(input_fasta_fp):
    """ Tests fasta filepath to determine if valid format
    
    input_fasta_fp:  fasta filepath
    """
    
    fasta_f = open(input_fasta_fp, "U")
    
    try:
        for label, seq in MinimalFastaParser(fasta_f):
            continue
    except RecordError:
        raise RecordError,("Input fasta file not valid fasta format.  Error "+\
         "found at %s label and %s sequence " % (label, seq))
        
    fasta_f.close()
        
        
def get_fasta_labels(input_fasta_fp):
    """ Returns the fasta labels (text before whitespace) as a list
    
    input_fasta_fp: fasta filepath
    """
    
    fasta_labels = []
    
    fasta_f = open(input_fasta_fp, "U")
    
    for label, seq in MinimalFastaParser(fasta_f):
        fasta_labels.append(label.split()[0])
        
    return fasta_labels
    
def get_dup_labels_perc(fasta_labels):
    """ Calculates percentage of sequences with duplicate labels
    
    fasta_labels: list of fasta labels
    """
    fasta_labels_count = len(fasta_labels)
    fasta_labels_derep = len(set(fasta_labels))
    
    perc_dup = "%1.3f" %\
     ((fasta_labels_count-fasta_labels_derep)/fasta_labels_count)
    return perc_dup
    
def check_labels_sampleids(fasta_labels,
                           sample_ids,
                           total_seq_count):
    """ Returns percent of valid fasta labels and that do not match SampleIDs
    
    fasta_labels: list of fasta labels
    sample_ids: set of sample IDs from mapping file
    total_seq_count: int of total sequences in fasta file
    """
    valid_id_count = 0
    matches_sampleid_count = 0
    
    for label in fasta_labels:
        curr_label = label.split('_')
        
        # Should be length 2, if not skip other processing
        if len(curr_label) != 2:
            continue
            
        valid_id_count += 1
        
        if curr_label[0] in sample_ids:
            matches_sampleid_count += 1
            
    perc_not_valid = "%1.3f" %\
     ((total_seq_count - valid_id_count)/total_seq_count)
    perc_nosampleid_match = "%1.3f" %\
     ((total_seq_count - matches_sampleid_count)/total_seq_count)
     
    return perc_not_valid, perc_nosampleid_match
    
def check_fasta_seqs(input_fasta_fp,
                     barcodes,
                     linkerprimerseqs,
                     total_seq_count,
                     valid_chars = frozenset(['A', 'T', 'C', 'G', 'N', 'a',
                      't', 'c', 'g', 'n'])):
    """ Returns perc of seqs w/ invalid chars, barcodes, or primers present
    
    input_fasta_fp:  fasta filepath
    barcodes: set of barcodes from the mapping file
    linkerprimerseqs: set of linkerprimersequences from the mapping file
    total_seq_count: total number of sequences in fasta file
    valid_chars: Currently allowed DNA chars
    """
    
    input_fasta_f = open(input_fasta_fp, "U")
    
    invalid_chars_count = 0
    barcodes_count = 0
    linkerprimers_count = 0
    
    for label,seq in MinimalFastaParser(input_fasta_f):
        
        # Only count one offending problem
        for curr_nt in seq:
            if curr_nt not in valid_chars:
                invalid_chars_count += 1
                break
        
        for curr_bc in barcodes:
            if curr_bc in seq:
                barcodes_count += 1
                break
                
        for curr_primer in linkerprimerseqs:
            if curr_primer in seq:
                linkerprimers_count += 1
                break
    
    perc_invalid_chars = "%1.3f" %\
     (invalid_chars_count/total_seq_count)
    perc_barcodes_detected = "%1.3f" %\
     (barcodes_count/total_seq_count)
    perc_primers_detected = "%1.3f" %\
     (linkerprimers_count/total_seq_count)
     
    return perc_invalid_chars, perc_barcodes_detected, perc_primers_detected
    
def run_fasta_checks(input_fasta_fp,
                     mapping_fp,
                     tree_fp=None,
                     tree_subset=False,
                     tree_exact_match=False,
                     same_seq_lens=False,
                     all_ids_found=False):
    """ Returns dictionary of records for different fasta checks
    
    input_fasta_fp: fasta filepath
    mapping_fp: mapping filepath
    tree_fp: newick tree filepath
    tree_subset: If True, will test that SampleIDs are a subset of the tree tips
    tree_exact_match: If True, will test that SampleIDs are an exact match to
     the tree
    same_seq_lens: If True, will determine if sequences are of different lens.
    all_ids_found: If True, will determine if all SampleIDs are represented
     in the sequence labels."""

    # Stores details of various checks
    fasta_report = {}
    
    # get sets of data for testing fasta labels/seqs
    sample_ids, barcodes, linkerprimerseqs = get_mapping_details(mapping_fp)
    
    fasta_labels = get_fasta_labels(input_fasta_fp)
    
    total_seq_count = len(fasta_labels)
    
    fasta_report['duplicate_labels'] = get_dup_labels_perc(fasta_labels)
     
    fasta_report['invalid_labels'], fasta_report['nosample_ids_map'] =\
     check_labels_sampleids(fasta_labels, sample_ids, total_seq_count)
     
    fasta_report['invalid_seq_chars'], fasta_report['barcodes_detected'],\
     fasta_report['linkerprimers_detected'] = check_fasta_seqs(input_fasta_fp,
     barcodes, linkerprimerseqs, total_seq_count)
     
    return fasta_report
     
def write_log_file(output_dir,
                   input_fasta_fp,
                   fasta_report):
    """ Formats data report, writes log
    
    output_dir: output directory
    input_fasta_fp: input fasta filepath, used to generate log name
    fasta_report: dictionary of percentages for different checks
    """
    
    output_fp = join(output_dir,
     split(input_fasta_fp)[1] + "_report.log")
     
    output_f = open(output_fp, "w")
    
    output_f.write("# fasta file %s validation report\n" % input_fasta_fp)
    
    output_f.write("Percent duplicate labels: %s\n" %\
     fasta_report['duplicate_labels'])
    output_f.write("Percent QIIME-incompatible fasta labels: %s\n" %\
     fasta_report['invalid_labels'])
    output_f.write("Percent of labels that fail to map to SampleIDs: %s\n" %\
     fasta_report['nosample_ids_map'])
    output_f.write("Percent of sequences with invalid characters: %s\n" %\
     fasta_report['invalid_seq_chars'])
    output_f.write("Percent of sequences with barcodes detected: %s\n" %\
     fasta_report['barcodes_detected'])
    output_f.write("Percent of sequences with primers detected: %s\n" %\
     fasta_report['linkerprimers_detected'])
     
    
    
        
    
def validate_fasta(input_fasta_fp,
                   mapping_fp,
                   output_dir,
                   tree_fp=None,
                   tree_subset=False,
                   tree_exact_match=False,
                   same_seq_lens=False,
                   all_ids_found=False):
    """ Main function for validating demultiplexed fasta file
    
    input_fasta_fp: fasta filepath
    mapping_fp: mapping filepath
    output_dir: output directory
    tree_fp: newick tree filepath
    tree_subset: If True, will test that SampleIDs are a subset of the tree tips
    tree_exact_match: If True, will test that SampleIDs are an exact match to
     the tree
    same_seq_lens: If True, will determine if sequences are of different lens.
    all_ids_found: If True, will determine if all SampleIDs are represented
     in the sequence labels.
    """

    # First test is valid fasta format, can't do other tests if file can't be
    # parsed so will bail out here if invalid
    
    verify_valid_fasta_format(input_fasta_fp)
    
    fasta_report = run_fasta_checks(input_fasta_fp, mapping_fp, tree_fp,
     tree_subset, tree_exact_match, same_seq_lens, all_ids_found)
     
    write_log_file (output_dir, input_fasta_fp, fasta_report)
    
    
        
    
    
    
    
    
