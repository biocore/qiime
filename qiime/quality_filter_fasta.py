#!/usr/bin/env python
from __future__ import division

# Handles quality filtering of demultiplexed fasta files, logs errors

__author__ = "Rob Knight and Micah Hamady"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ =  ["Rob Knight", "Micah Hamady", "Greg Caporaso",
                "Kyle Bittinger", "Jesse Stombaugh","William Walters",
                "Jens Reeder", "Jose Antonio Navas Molina",
                "Jai Ram Rideout"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Development"

from string import upper
from itertools import izip
from os.path import join
from os import rename
from gzip import GzipFile
from warnings import filterwarnings
filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from numpy import arange, histogram, mean
from cogent import DNA, LoadSeqs
from cogent.align.align import make_dna_scoring_dict, local_pairwise
from cogent.parse.fasta import MinimalFastaParser
from cogent.seqsim.sequence_generators import SequenceGenerator, IUPAC_DNA
from cogent.core.moltype import IUPAC_DNA_ambiguities

from qiime.parse import QiimeParseError, MinimalQualParser
from qiime.check_id_map import process_id_map
from qiime.format import format_histograms_two_categories
from qiime.util import qiime_open

def quality_filter_sequences(mapping_file,
                             fasta_files,
                             qual_files = set(),
                             min_seq_len = 200,
                             max_seq_len = 1000,
                             min_qual_score = 25,
                             retain_primer = False,
                             max_ambig = 6,
                             suppress_ambig_check = False,
                             max_homopolymer = 6,
                             suppress_homopolymer_check = False,
                             max_primer_mismatch = 0,
                             output_dir = ".",
                             qual_score_window = 0,
                             discard_bad_windows = False,
                             suppress_primer_check = False,
                             reverse_primers = "disable",
                             reverse_primer_mismatches = 0,
                             record_qual_scores = False,
                             truncate_ambi_bases = False,
                             suppress_sampleid_check = False,
                             enable_all_checks = False,
                             local_align_forward_primer = False):
    """ Main function for quality filtering sequences
    
    mapping_file: Filepath for metadata mapping file
    fasta_files: set of fasta filepaths
    qual_files: set of qual filepaths
    min_seq_len: Minimum sequence length
    max_seq_len: Max sequence length
    min_qual_score: Minimum overall qual score, also used for sliding window
    keep_primer: If True, forward primer is not removed from output
    max_ambig: Max ambiguous "N" bases allowed
    suppress_ambig_check: If True, disables ambiguous base check
    max_homopolymer: Max homopolymer length allowed
    suppress_homopolymer_check: If True, disabled homopolymer check
    max_primer_mismatch: Maximum allowed forward primer mismatches
    output_dir: output directory
    qual_score_window: Sliding qual window check, 0 disables, any other size
     is the size of the window
    discard_bad_windows: If True, and low qual sliding window detected, 
     discards the sequence rather than writing
    suppress_primer_check: If True, suppresses forward primer mismatch and
     removal check
    reverse_primers: Can be 'disable', 'truncate_only', and 'truncate_remove'.
     Uses ReversePrimer from mapping data to do local alignment of reverse
     primer, and will truncate at that point in the sequence. If 
     truncate_remove option is enabled, will discard reads where reverse
     primer can not be found.
    reverse_primer_mismatches: Allowed mismatches in reverse primer.
    record_qual_scores: If True, will write matching quality score sequence to
     trimmed fasta file.
    truncate_ambi_bases: If True, will ignore max_homopolymer value, and instead
     will truncate at first ambiguous "N" base encountered.
    suppress_sampleid_check: If True, will not test for unique SampleIDs in 
     the mapping file, and will use all forward and (if present) reverse 
     primers in the mapping file during tests. 
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    local_align_forward_primer: If True, will attempt to do local alignment to
     find the forward primer rather than slicing initial bases.
    """
    mapping_data = open(mapping_file, "U")
    header, mapping_data = check_map(mapping_data, suppress_primer_check,
     reverse_primers, suppress_sampleid_check) 
    ids_primers, ids_rev_primers = get_ids_primers(header, mapping_data,
     reverse_primers, suppress_sampleid_check)

    final_log_data, detailed_quality_data, seq_length_data =\
     process_seqs(fasta_files,
     qual_files, ids_primers, ids_rev_primers, min_seq_len,
     max_seq_len, min_qual_score, retain_primer, max_ambig,
     suppress_ambig_check, max_homopolymer, suppress_homopolymer_check,
     max_primer_mismatch, output_dir, qual_score_window, discard_bad_windows,
     suppress_primer_check, reverse_primers, reverse_primer_mismatches,
     record_qual_scores, truncate_ambi_bases, suppress_sampleid_check,
     enable_all_checks, local_align_forward_primer)
        
    formatted_log_data = format_log_data(final_log_data, seq_length_data,
     min_seq_len, max_seq_len, min_qual_score, retain_primer, max_ambig,
     suppress_ambig_check, max_homopolymer, suppress_homopolymer_check,
     max_primer_mismatch, output_dir, qual_score_window, discard_bad_windows,
     suppress_primer_check, reverse_primers, reverse_primer_mismatches,
     record_qual_scores, truncate_ambi_bases, suppress_sampleid_check,
     enable_all_checks)
    log_file = open(join(output_dir, "quality_filter_log.txt"), "w")
    log_file.write(formatted_log_data)
    log_file.close()
    
    if enable_all_checks:
        formatted_detailed_data = format_detailed_data(detailed_quality_data)
        detailed_quality_report = open(join(output_dir,
         "detailed_quality_report.txt"), "w")
        detailed_quality_report.write(formatted_detailed_data)
        detailed_quality_report.close()
        
    histogram_data =\
     format_histograms_two_categories(*make_histograms(seq_length_data['raw'],
     seq_length_data['processed']))
    hist_file = open(join(output_dir, "histograms.txt"), "w")
    hist_file.write(histogram_data)
    hist_file.close()

    return
    
# Quality filtering functions ******

def process_seqs(fasta_files,
                 qual_files,
                 ids_primers,
                 ids_rev_primers,
                 min_seq_len = 200,
                 max_seq_len = 1000,
                 min_qual_score = 25,
                 retain_primer = False,
                 max_ambig = 6,
                 suppress_ambig_check = False,
                 max_homopolymer = 6,
                 suppress_homopolymer_check = False,
                 max_primer_mismatch = 0,
                 output_dir = ".",
                 qual_score_window = 0,
                 discard_bad_windows = False,
                 suppress_primer_check = False,
                 reverse_primers = "disable",
                 reverse_primer_mismatches = 0,
                 record_qual_scores = False,
                 truncate_ambi_bases = False,
                 suppress_sampleid_check = False,
                 enable_all_checks = False,
                 local_align_forward_primer = False):
    """ Iterates through fasta/qual files, quality checks, writes seqs
    
    fasta_files: set of input fasta filepaths
    qual_files: set of input qual filepaths
    ids_primers: dict of SampleIDs:forward primers
    ids_rev_primers: dict of SampleIDs:reverse primers
    min_seq_len: Minimum sequence length
    max_seq_len: Max sequence length
    min_qual_score: Minimum overall qual score, also used for sliding window
    keep_primer: If True, forward primer is not removed from output
    max_ambig: Max ambiguous "N" bases allowed
    suppress_ambig_check: If True, disables ambiguous base check
    max_homopolymer: Max homopolymer length allowed
    suppress_homopolymer_check: If True, disabled homopolymer check
    max_primer_mismatch: Maximum allowed forward primer mismatches
    output_dir: output directory
    qual_score_window: Sliding qual window check, 0 disables, any other size
     is the size of the window
    discard_bad_windows: If True, and low qual sliding window detected, 
     discards the sequence rather than writing
    suppress_primer_check: If True, suppresses forward primer mismatch and
     removal check
    reverse_primers: Can be 'disable', 'truncate_only', and 'truncate_remove'.
     Uses ReversePrimer from mapping data to do local alignment of reverse
     primer, and will truncate at that point in the sequence. If 
     truncate_remove option is enabled, will discard reads where reverse
     primer can not be found.
    reverse_primer_mismatches: Allowed mismatches in reverse primer.
    record_qual_scores: If True, will write matching quality score sequence to
     trimmed fasta file.
    truncate_ambi_bases: If True, will ignore max_homopolymer value, and instead
     will truncate at first ambiguous "N" base encountered.
    suppress_sampleid_check: If True, will not test for unique SampleIDs in 
     the mapping file, and will use all forward and (if present) reverse 
     primers in the mapping file during tests. 
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    local_align_forward_primer: If True, will attempt to do local alignment to
     find the forward primer rather than slicing initial bases.
    """
    
    # initialize log data
    final_log_data = {'seq_counts':0,
                      'seqs_written':0,
                      'fasta_files':fasta_files,
                      'qual_files':qual_files,
                      'below_min_seq_len':0,
                      'max_seq_len':0,
                      'min_ave_qual_score':0,
                      'max_ambig':0,
                      'ambig_base_truncation':0,
                      'too_short_after_ambig_truncation':0,
                      'max_homopolymer':0,
                      'max_primer_mismatch':0,
                      'low_qual_window_found':0,
                      'too_short_after_window_truncation':0,
                      'low_qual_window_discarded':0,
                      'rev_primers_found':0,
                      'max_rev_primer_mismatch':0,
                      'seqs_discarded_no_rev_primer':0,
                      'too_short_after_revprimer_truncation':0,
                      'seq_ids_not_in_mapping':[]
                      }
    
    detailed_quality_data = {'seq_order':[]}
    detailed_quality_data['all_tests'] = ['below_min_seq_len',
    'exceeds_max_seq_len', 'exceeds_max_primer_mismatch',
    'rev_primer_found', 'exceeds_max_rev_primer_mismatch',
    'seqs_discarded_for_no_rev_primer', 'too_short_after_revprimer_truncation',
    'truncation_for_ambig_base', 'too_short_after_ambig_truncation',
    'exceeds_max_ambig', 'exceeds_max_homopolymer', 'below_min_ave_qual_score',
    'low_qual_window_found', 'discarded_for_low_qual_window',
    'too_short_after_window_truncation', 'seq_id_not_in_mapping']
        
    seq_length_data = {'raw':[], 'processed':[]}

    # empty qual seq for case of no qual file use
    qual_seq=""
    
    fasta_f = [qiime_open(fasta_f) for fasta_f in fasta_files]
    qual_f = [qiime_open(qual_f) for qual_f in qual_files]
    
    quality_filtered_seqs_f = open(join(output_dir,
     "quality_filtered_seqs.fna.incomplete"), "w")
    if qual_files and record_qual_scores:
        quality_filtered_seqs_q = open(join(output_dir,
         "quality_filtered_seqs.qual.incomplete"), "w")
         
    # Possibly use empty array and avoid second for loop for fasta only situation
         
    if qual_files:
        for curr_fasta, curr_qual in zip(fasta_f, qual_f):
            for fasta_data, qual_data in izip(MinimalFastaParser(curr_fasta),
             MinimalQualParser(curr_qual, full_header=True)):

                final_log_data['seq_counts'] += 1
                fasta_label, fasta_seq = fasta_data
                qual_label, qual_seq = qual_data
                
                seq_length_data['raw'].append(len(fasta_seq))
                
                # Check that labels, fasta/qual length equal
                if fasta_label != qual_label:
                    raise ValueError,("Found fasta and qual labels that do "
                     "not match: fasta label %s, qual label %s" % (fasta_label,
                     qual_label))
                if len(fasta_seq) != len(qual_seq):
                    raise ValueError,("Found fasta and qual seqs that do "
                     "not have matching lengths: fasta label "
                     "%s, qual label %s" % (fasta_label, qual_label))
                final_fasta, final_qual, final_log_data, detailed_quality_data=\
                 quality_filter_seq(fasta_label, fasta_seq.upper(), qual_seq,
                 final_log_data, detailed_quality_data, ids_primers,
                 ids_rev_primers, min_seq_len, max_seq_len, min_qual_score,
                 retain_primer, max_ambig, suppress_ambig_check,
                 max_homopolymer, suppress_homopolymer_check,
                 max_primer_mismatch, qual_score_window, discard_bad_windows,
                 suppress_primer_check, reverse_primers,
                 reverse_primer_mismatches, truncate_ambi_bases,
                 suppress_sampleid_check, enable_all_checks,
                 local_align_forward_primer)
                 
                if final_fasta:
                    final_log_data['seqs_written'] += 1
                    quality_filtered_seqs_f.write(">%s\n%s\n" %\
                     (fasta_label, final_fasta))
                    seq_length_data['processed'].append(len(final_fasta))
                    if record_qual_scores:
                        write_qual_line(quality_filtered_seqs_q, qual_label,
                         final_qual)
                                 
    else:
        for curr_fasta in fasta_f:
            for fasta_label, fasta_seq in MinimalFastaParser(curr_fasta):
                
                final_log_data['seq_counts'] += 1
                
                seq_length_data['raw'].append(len(fasta_seq))
                
                final_fasta, final_qual, final_log_data, detailed_quality_data=\
                 quality_filter_seq(fasta_label, fasta_seq.upper(), qual_seq,
                 final_log_data, detailed_quality_data, ids_primers,
                 ids_rev_primers, min_seq_len, max_seq_len, min_qual_score,
                 retain_primer, max_ambig, suppress_ambig_check,
                 max_homopolymer, suppress_homopolymer_check,
                 max_primer_mismatch, qual_score_window, discard_bad_windows,
                 suppress_primer_check, reverse_primers,
                 reverse_primer_mismatches, truncate_ambi_bases,
                 suppress_sampleid_check, enable_all_checks,
                 local_align_forward_primer)
                 
                
                if final_fasta:
                    final_log_data['seqs_written'] += 1
                    quality_filtered_seqs_f.write(">%s\n%s\n" %\
                     (fasta_label, final_fasta))
                    seq_length_data['processed'].append(len(final_fasta))
         
    # Rename .incomplete files to .fna/.qual files
    rename(quality_filtered_seqs_f.name,join(output_dir, "seqs.fna"))
    if qual_files and record_qual_scores:
        rename(quality_filtered_seqs_q.name, join(output_dir, "seqs.qual"))
        
    return final_log_data, detailed_quality_data, seq_length_data
    
def quality_filter_seq(fasta_label,
                       fasta_seq,
                       qual_seq,
                       final_log_data,
                       detailed_quality_data,
                       ids_primers,
                       ids_rev_primers,
                       min_seq_len = 200,
                       max_seq_len = 1000,
                       min_qual_score = 25,
                       retain_primer = False,
                       max_ambig = 6,
                       suppress_ambig_check = False,
                       max_homopolymer = 6,
                       suppress_homopolymer_check = False,
                       max_primer_mismatch = 0,
                       qual_score_window = 0,
                       discard_bad_windows = False,
                       suppress_primer_check = False,
                       reverse_primers = "disable",
                       reverse_primer_mismatches = 0,
                       truncate_ambi_bases = False,
                       suppress_sampleid_check = False,
                       enable_all_checks = False,
                       local_align_forward_primer = False):
    """ Does filtering checks for sequence/qual data passed
    
    Returns truncated fasta, qual (empty string if discarded) and added log 
    data, detailed log data if specified
    
    fasta_label: Fasta sequence ID, used for detailed logging, primer matching
    fasta_seq: Fasta sequence to quality check, string of nucleotides
    qual_seq: Associated quality scores for fasta sequence, array of ints
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: list of lines containing detailed per seq data
    ids_primers: dict of SampleIDs:forward primers
    ids_rev_primers: dict of SampleIDs:reverse primers
    min_seq_len: Minimum sequence length
    max_seq_len: Max sequence length
    min_qual_score: Minimum overall qual score, also used for sliding window
    keep_primer: If True, forward primer is not removed from output
    max_ambig: Max ambiguous "N" bases allowed
    suppress_ambig_check: If True, disables ambiguous base check
    max_homopolymer: Max homopolymer length allowed
    suppress_homopolymer_check: If True, disabled homopolymer check
    max_primer_mismatch: Maximum allowed forward primer mismatches
    qual_score_window: Sliding qual window check, 0 disables, any other size
     is the size of the window
    discard_bad_windows: If True, and low qual sliding window detected, 
     discards the sequence rather than writing
    suppress_primer_check: If True, suppresses forward primer mismatch and
     removal check
    reverse_primers: Can be 'disable', 'truncate_only', and 'truncate_remove'.
     Uses ReversePrimer from mapping data to do local alignment of reverse
     primer, and will truncate at that point in the sequence. If 
     truncate_remove option is enabled, will discard reads where reverse
     primer can not be found.
    reverse_primer_mismatches: Allowed mismatches in reverse primer.
    record_qual_scores: If True, will write matching quality score sequence to
     trimmed fasta file.
    truncate_ambi_bases: If True, will ignore max_homopolymer value, and instead
     will truncate at first ambiguous "N" base encountered.
    suppress_sampleid_check: If True, will not test for unique SampleIDs in 
     the mapping file, and will use all forward and (if present) reverse 
     primers in the mapping file during tests. 
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    local_align_forward_primer: If True, will attempt to do local alignment to
     find the forward primer rather than slicing initial bases.
    """
    
    if enable_all_checks:
        detailed_quality_data['seq_order'].append(fasta_label) 
    failures = False
    detailed_quality_data[fasta_label] = {}
    
    if not suppress_primer_check:
        failed, fasta_seq, qual_seq, final_log_data =\
         check_primer_mismatch(fasta_label, fasta_seq, qual_seq, final_log_data,
          detailed_quality_data[fasta_label], ids_primers, retain_primer,
          max_primer_mismatch, suppress_sampleid_check, enable_all_checks,
          local_align_forward_primer)
        failures = failures | failed
        if failed and not enable_all_checks:
            return "", "", final_log_data, detailed_quality_data
    
    failed, final_log_data = check_min_seq_length(fasta_seq, final_log_data,
     detailed_quality_data[fasta_label], min_seq_len, enable_all_checks)
    failures = failures | failed
    if failed and not enable_all_checks:
        return "", "", final_log_data, detailed_quality_data
    
    failed, final_log_data = check_max_seq_length(fasta_seq, final_log_data,
      detailed_quality_data[fasta_label], max_seq_len, enable_all_checks)
    failures = failures | failed
    if failed and not enable_all_checks:
        return "", "", final_log_data, detailed_quality_data
        
    if reverse_primers in ['truncate_only', 'truncate_remove']:
        failed, fasta_seq, qual_seq, final_log_data =\
         check_rev_primer_mismatch(fasta_label, fasta_seq, qual_seq,
          final_log_data, detailed_quality_data[fasta_label], ids_rev_primers,
          reverse_primer_mismatches, reverse_primers,
          suppress_sampleid_check, enable_all_checks)
        failures = failures | failed
        if failed and not enable_all_checks:
            return "", "", final_log_data, detailed_quality_data
            
        # Have to check if sequence too short after truncation
        failed, final_log_data = check_min_seq_length(fasta_seq,
         final_log_data, detailed_quality_data[fasta_label], min_seq_len,
         enable_all_checks, 'too_short_after_revprimer_truncation')
        failures = failures | failed
        if failed and not enable_all_checks:
            return "", "", final_log_data, detailed_quality_data
    
    if truncate_ambi_bases:
        fasta_seq, qual_seq, final_log_data, truncated=\
         truncate_ambig_bases(fasta_seq, qual_seq, final_log_data,
         detailed_quality_data[fasta_label], enable_all_checks)
         
        # Have to check if sequence too short after truncation
        if truncated:
            failed, final_log_data = check_min_seq_length(fasta_seq,
             final_log_data, detailed_quality_data[fasta_label], min_seq_len,
             enable_all_checks, 'too_short_after_ambig_truncation')
            failures = failures | failed
            if failed and not enable_all_checks:
                return "", "", final_log_data, detailed_quality_data
    
    if qual_score_window:
        failed, fasta_seq, qual_seq, final_log_data =\
         check_sliding_qual_window(fasta_seq, qual_seq, final_log_data,
         detailed_quality_data[fasta_label], qual_score_window, min_qual_score,
         discard_bad_windows, enable_all_checks)
        failures = failures | failed
        if failed and not enable_all_checks:
            return "", "", final_log_data, detailed_quality_data
         
        # Have to check if sequence too short after truncation
        failed, final_log_data = check_min_seq_length(fasta_seq,
         final_log_data, detailed_quality_data[fasta_label], min_seq_len,
         enable_all_checks, 'too_short_after_window_truncation')
        failures = failures | failed
        if failed and not enable_all_checks:
            return "", "", final_log_data, detailed_quality_data
    
    if len(qual_seq) > 0:
        failed, final_log_data = check_average_quality(qual_seq,
         final_log_data, detailed_quality_data[fasta_label], min_qual_score,
         enable_all_checks)
        failures = failures | failed
        if failed and not enable_all_checks:
            return "", "", final_log_data, detailed_quality_data
    
    if not suppress_ambig_check:
        failed, final_log_data = check_ambig_count(fasta_seq,
         final_log_data, detailed_quality_data[fasta_label], max_ambig,
         enable_all_checks)
        failures = failures | failed
        if failed and not enable_all_checks:
            return "", "", final_log_data, detailed_quality_data
            
    if not suppress_homopolymer_check:
        failed, final_log_data = check_homopolymers(fasta_seq,
         final_log_data, detailed_quality_data[fasta_label], max_homopolymer,
         enable_all_checks)
        failures = failures | failed
        if failed and not enable_all_checks:
            return "", "", final_log_data, detailed_quality_data
            
    if enable_all_checks:
        # If any failures, don't write the sequence
        if failures:
            return "", "", final_log_data, detailed_quality_data
            
    # Should only return this if everything passes filters
    return fasta_seq, qual_seq, final_log_data, detailed_quality_data
    

def check_primer_mismatch(fasta_label,
                          fasta_seq,
                          qual_seq,
                          final_log_data,
                          detailed_quality_data,
                          ids_primers,
                          retain_primer = False,
                          max_primer_mismatch = 0,
                          suppress_sampleid_check = False,
                          enable_all_checks = False,
                          local_align_forward_primer = False):
    """ Checks for max allowed primer mismatches, strips primer sequence
    
    fasta_label: Fasta sequence ID, used for detailed logging, primer matching
    fasta_seq: fasta sequence, string of nucleotide characters
    qual_seq: quality scores, array of ints
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: current dictionary from detailed log data, keyed
     to the current label ID being checked.
    ids_primers: dictionary of SampleID:list of primers
    retain_primer: If True, do not remove the primer from the beginning of the 
     read.
    max_primer_mismatch: Maximum primer mismatches allowed.
    suppress_sampleid_check: If True, will not test for unique SampleIDs in 
     the mapping file, and will use all forward and (if present) reverse 
     primers in the mapping file during tests.
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    local_align_forward_primer: If True, will attempt to do local alignment to
     find the forward primer rather than slicing initial bases
    """ 
    
    # If SampleID doesn't match known primers, don't proceed with primer check
    # and log the sequence label
    if not suppress_sampleid_check:
        seq_label = fasta_label.split('_')[0]
        if seq_label not in ids_primers.keys():
            final_log_data['seq_ids_not_in_mapping'].append(fasta_label)
            if enable_all_checks:
                detailed_quality_data['seq_id_not_in_mapping'] = 1
            return False, fasta_seq, qual_seq, final_log_data
        else:
            primers = ids_primers[seq_label]
            if enable_all_checks:
                detailed_quality_data['seq_id_not_in_mapping'] = 0
    else:
        primers = ids_primers['all_primers']
    
    failed = True        

    if local_align_forward_primer:
        for primer in primers:
            mismatches, hit_start = local_align_primer_seq(primer, fasta_seq)
            if mismatches <= max_primer_mismatch:
                if enable_all_checks:
                    detailed_quality_data['exceeds_max_primer_mismatch'] = 0
                fasta_seq = fasta_seq[hit_start + len(primer):]
                qual_seq = qual_seq[hit_start + len(primer):]
                failed = False
                break
    else:
        for primer in primers:
            exceeds_mismatch = count_mismatches(fasta_seq, primer,
             max_primer_mismatch)
            if not exceeds_mismatch:
                if enable_all_checks:
                    detailed_quality_data['exceeds_max_primer_mismatch'] = 0
                if not retain_primer:
                    fasta_seq = fasta_seq[len(primer):]
                    qual_seq = qual_seq[len(primer):]
                failed = False
                break
    
    if failed:
        final_log_data['max_primer_mismatch'] += 1
        if enable_all_checks:
            detailed_quality_data['exceeds_max_primer_mismatch'] = 1

    return failed, fasta_seq, qual_seq, final_log_data
    
def check_rev_primer_mismatch(fasta_label,
                              fasta_seq,
                              qual_seq,
                              final_log_data,
                              detailed_quality_data,
                              ids_rev_primers,
                              reverse_primer_mismatches = 0,
                              reverse_primers = "truncate_only",
                              suppress_sampleid_check = False,
                              enable_all_checks = False):
    """ Finds, checks reverse primer for mismatches
    
    fasta_label: Fasta sequence ID, used for detailed logging, primer matching
    fasta_seq: fasta sequence, string of nucleotide characters
    qual_seq: quality scores, array of ints
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: current dictionary from detailed log data, keyed
     to the current label ID being checked.
    ids_rev_primers: dictionary of SampleID:list of reverse primers
    reverse_primer_mismatches: Maximum primer mismatches allowed.
    reverse_primers: Can be 'truncate_only' or 'truncate_remove'. truncate_only
     will try to find the primer, and truncate at that point in the sequence.
     If truncate_remove option is enabled, will discard reads where reverse
     primer can not be found.
    suppress_sampleid_check: If True, will not test for unique SampleIDs in 
     the mapping file, and will use all forward and (if present) reverse 
     primers in the mapping file during tests.
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    """
    
    # If SampleID doesn't match known primers, don't proceed with primer check
    # and log the sequence label
    if not suppress_sampleid_check:
        curr_ids_not_found = set(final_log_data['seq_ids_not_in_mapping'])
        seq_label = fasta_label.split('_')[0]
        if seq_label not in ids_rev_primers.keys():
            # Doing this check to avoid double logging for both forward and
            # reverse primer checks
            if fasta_label not in curr_ids_not_found:
                final_log_data['seq_ids_not_in_mapping'].append(fasta_label)
            if enable_all_checks:
                detailed_quality_data['seq_id_not_in_mapping'] = 1
            return False, fasta_seq, qual_seq, final_log_data
        else:
            primers = ids_rev_primers[seq_label]
            if enable_all_checks:
                detailed_quality_data['seq_id_not_in_mapping'] = 0
    else:
        primers = ids_rev_primers['rev_primers']
    
    failed = True        
    for primer in primers:
        mismatches, hit_start = local_align_primer_seq(primer, fasta_seq)
        if mismatches <= reverse_primer_mismatches:
            if enable_all_checks:
                detailed_quality_data['exceeds_max_rev_primer_mismatch'] = 0
                detailed_quality_data['rev_primer_found'] = 1
                if reverse_primers == 'truncate_remove':
                    detailed_quality_data['seqs_discarded_for_no_rev_primer']=0
            final_log_data['rev_primers_found'] += 1
            fasta_seq = fasta_seq[0:hit_start]
            qual_seq = qual_seq[0:hit_start]
            failed = False
            break
    
    if failed:
        final_log_data['max_rev_primer_mismatch'] += 1

        if enable_all_checks:
            detailed_quality_data['exceeds_max_rev_primer_mismatch'] = 1
            detailed_quality_data['rev_primer_found'] = 0
            
        # Only should fail and not write sequence if truncate_remove is True
        if reverse_primers == 'truncate_remove':
            failed = True
            final_log_data['seqs_discarded_no_rev_primer'] += 1
            if enable_all_checks:
                detailed_quality_data['seqs_discarded_for_no_rev_primer'] = 1
        else:
            failed = False

    return failed, fasta_seq, qual_seq, final_log_data
    
    
def check_min_seq_length(fasta_seq,
                         final_log_data,
                         detailed_quality_data,
                         min_seq_len = 200,
                         enable_all_checks = False,
                         curr_log_key = 'below_min_seq_len'
                         ):
    """ Checks minimum sequence length, logs details if fails check
    
    fasta_seq: Fasta sequence to quality check, string of nucleotides
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: current dictionary from detailed log data, keyed
     to the current label ID being checked.
    min_seq_len: Minimum sequence length
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    """
    
    if len(fasta_seq) < min_seq_len:
        failed = True
        final_log_data[curr_log_key] += 1
        if enable_all_checks:
            detailed_quality_data[curr_log_key] = 1
    else:
        failed = False
        if enable_all_checks:
            detailed_quality_data[curr_log_key] = 0

    return failed, final_log_data
    
def check_max_seq_length(fasta_seq,
                         final_log_data,
                         detailed_quality_data,
                         max_seq_len = 1000,
                         enable_all_checks = False
                         ):
    """ Checks minimum sequence length, logs details if fails check
    
    fasta_seq: Fasta sequence to quality check, string of nucleotides
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: current dictionary from detailed log data, keyed
     to the current label ID being checked.
    max_seq_len: Maximum sequence length allowed
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    """
    if len(fasta_seq) > max_seq_len:
        failed = True
        final_log_data['max_seq_len'] += 1
        if enable_all_checks:
            detailed_quality_data['exceeds_max_seq_len'] = 1
    else:
        failed = False
        if enable_all_checks:
            detailed_quality_data['exceeds_max_seq_len'] = 0

    return failed, final_log_data
    
def truncate_ambig_bases(fasta_seq,
                         qual_seq,
                         final_log_data,
                         detailed_quality_data,
                         enable_all_checks = False,
                         ):
    """ Truncates fasta and qual seq at first ambiguous base position if found
    
    fasta_seq: Fasta sequence to quality check, string of nucleotides
    qual_seq: Quality sequence, array of ints.
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: current dictionary from detailed log data, keyed
     to the current label ID being checked.
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    """
    ambig_bases = "RYMKWSBDHVN"
    ambig_trunc_lens = []
    
    for curr_base in ambig_bases:
        try:
            trunc_ix = fasta_seq.index(curr_base)
            ambig_trunc_lens.append(trunc_ix)
        except ValueError:
            continue
            
    if ambig_trunc_lens:
        min_len_trunc = min(ambig_trunc_lens)
        fasta_seq = fasta_seq[0:min_len_trunc]
        qual_seq = qual_seq[0:min_len_trunc]
        final_log_data['ambig_base_truncation'] += 1
        if enable_all_checks:
            detailed_quality_data['truncation_for_ambig_base'] = 1
        truncated = True
    else:
        if enable_all_checks:
            detailed_quality_data['truncation_for_ambig_base'] = 0
        truncated = False

    return fasta_seq, qual_seq, final_log_data, truncated
    
def check_ambig_count(fasta_seq,
                      final_log_data,
                      detailed_quality_data,
                      max_ambig = 6,
                      enable_all_checks = False):
    """ Checks for max allowed ambiguous characters in sequence
    
    fasta_seq: Fasta sequence to quality check, string of nucleotides
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: current dictionary from detailed log data, keyed
     to the current label ID being checked.
    max_ambig: Maximum allowed ambiguous characters
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    """
    
    ambig_counts = count_ambig(fasta_seq)
    if ambig_counts > max_ambig:
        failed = True
        final_log_data['max_ambig'] += 1
        if enable_all_checks:
            detailed_quality_data['exceeds_max_ambig'] = 1
    else:
        failed = False
        if enable_all_checks:
            detailed_quality_data['exceeds_max_ambig'] = 0

    return failed, final_log_data
    
def check_sliding_qual_window(fasta_seq,
                              qual_seq,
                              final_log_data,
                              detailed_quality_data,
                              qual_score_window = 50,
                              min_qual_score = 25,
                              discard_bad_windows = False,
                              enable_all_checks = False):
    """ Checks sliding window of quality scores, truncates at first window found
    
    fasta_seq: dna sequence, string of nucleotides
    qual_seq: quality scores, array of ints
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: current dictionary from detailed log data, keyed
     to the current label ID being checked.
    qual_score_window: size of quality score window to check
    min_qual_score: Minimum quality score allowed.
    discard_bad_windows: If True, will discard sequence where low quality
     window found.
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    """
    
    found_window, index = check_window_qual_scores(qual_seq, qual_score_window,
     min_qual_score)
     
    if found_window:
        final_log_data['low_qual_window_found'] += 1
        if enable_all_checks:
            detailed_quality_data['low_qual_window_found'] = 1
        if discard_bad_windows:
            failed = True
            final_log_data['low_qual_window_discarded'] += 1
            if enable_all_checks:
                detailed_quality_data['discarded_for_low_qual_window'] = 1
        else:
            failed = False
            fasta_seq = fasta_seq[0:index]
            qual_seq = qual_seq[0:index]

    else:
        failed = False
        if enable_all_checks:
            detailed_quality_data['low_qual_window_found'] = 0
    
    return failed, fasta_seq, qual_seq, final_log_data
    
def check_average_quality(qual_seq,
                          final_log_data,
                          detailed_quality_data,
                          min_qual_score = 25,
                          enable_all_checks = False):
    """ Checks for max allowed ambiguous characters in sequence
    
    qual_seq: quality scores, array of ints
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: current dictionary from detailed log data, keyed
     to the current label ID being checked.
    min_qual_score: Minimum quality score allowed.
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    """ 
    if mean(qual_seq) < min_qual_score:
        failed = True
        final_log_data['min_ave_qual_score'] += 1
        if enable_all_checks:
            detailed_quality_data['below_min_ave_qual_score'] = 1
    else:
        failed = False
        if enable_all_checks:
            detailed_quality_data['below_min_ave_qual_score'] = 0

    return failed, final_log_data
    
def check_homopolymers(fasta_seq,
                       final_log_data,
                       detailed_quality_data,
                       max_homopolymer = 6,
                       enable_all_checks = False):
    """ Checks for homopolymer runs that exceed max allowed length
    
    fasta_seq: Fasta sequence to quality check, string of nucleotides
    final_log_data: Dict that contains details of quality filter counts
    detailed_quality_data: current dictionary from detailed log data, keyed
     to the current label ID being checked.
    max_homopolymer: Maximum allowed homopolymer run.
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    """
    
    homopolymer_found = seq_exceeds_homopolymers(fasta_seq, max_homopolymer)
    
    if homopolymer_found:
        failed = True
        final_log_data['max_homopolymer'] += 1
        if enable_all_checks:
            detailed_quality_data['exceeds_max_homopolymer'] = 1
    else:
        failed = False
        if enable_all_checks:
            detailed_quality_data['exceeds_max_homopolymer'] = 0
    
    return failed, final_log_data
        
# End quality filtering functions ***********  

# Mapping file data processing **************
    
def check_map(mapping_data,
              suppress_primer_check = False,
              reverse_primers = 'disable',
              suppress_sampleid_check = False):
    """ Gets header, mapping data, halts execution if there are errors
    
    mapping_file:  list of lines of metadata mapping file
    suppress_primer_check: If True, suppresses forward primer mismatch and
     removal check
    reverse_primers: Can be 'disable', 'truncate_only', and 'truncate_remove'.
     Uses ReversePrimer from mapping data to do local alignment of reverse
     primer, and will truncate at that point in the sequence. If 
     truncate_remove option is enabled, will discard reads where reverse
     primer can not be found.
    suppress_sampleid_check: If True, will not test for unique SampleIDs in 
     the mapping file, and will use all forward and (if present) reverse 
     primers in the mapping file during tests.
    """
    
    header, mapping_data, run_description, errors, warnings= \
        process_id_map(mapping_data, disable_primer_check=suppress_primer_check)
        
    halt_for_errors = False
    sample_id_col = 0
    primer_col = 2
    if reverse_primers != 'disable':
        rev_primer_test = True
        try:
            rev_primer_col = header.index("ReversePrimer")
        except ValueError:
            raise ValueError,('ReversePrimer not found in mapping file header,'
             ' please checking mapping file and add ReversePrimer column and '
             'values if necessary.')
    else:
        rev_primer_test = False
    
    for curr_error in errors:
        # Check for bad primer sequence based on column
        curr_col = int(curr_error.split('\t')[1].split(',')[1])
        if curr_col == primer_col:
            halt_for_errors = True
        if rev_primer_test:
            if curr_col == rev_primer_col:
                halt_for_errors = True
        # Check for duplicate SampleIDs
        if not suppress_sampleid_check:
            if curr_col == sample_id_col:
                if "Duplicate SampleID" in curr_error:
                    halt_for_errors = True
        # Check for invalid headers
        if "Found header field" in curr_error:
            halt_for_errors = True
        
    """This should only halt for duplicate SampleIDs (unless check suppressed),
    missing primers, invalid characters in primers, or errors in mapping
    header."""
    if halt_for_errors:
        raise ValueError,("Errors found in mapping file, please check "+\
         "mapping file with check_id_map.py. See suppress_primer_check and "+\
         "suppress_sampleid_check options if primers or SampleIDs should "+\
         "not be checked.")
    
    return header, mapping_data
    
def get_ids_primers(header,
                    mapping_data,
                    reverse_primers = "disable",
                    suppress_sampleid_check = False):
    """ Creates dict(s) of SampleIDs:primers (forward and reverse)
    
    header:  list of strings of header data.
    mapping_data:  list of lists of mapping data
    reverse_primers: Can be 'disable', 'truncate_only', and 'truncate_remove'.
     Uses ReversePrimer from mapping data to do local alignment of reverse
     primer, and will truncate at that point in the sequence. If 
     truncate_remove option is enabled, will discard reads where reverse
     primer can not be found.
    suppress_sampleid_check: If True, will not test for unique SampleIDs in 
     the mapping file, and will use all forward and (if present) reverse 
     primers in the mapping file during tests.
    
    Dict is of format SampleID:primer sequence
    """
    sample_id_ix = header.index("SampleID")
    primer_ix = header.index("LinkerPrimerSequence")
    if reverse_primers != "disable":
        rev_primer_ix = header.index("ReversePrimer")
    
    ids_primers = {}
    ids_rev_primers = {}
    all_primers = set([])
    all_rev_primers = set([])
       
    for line in mapping_data:
        all_forward_primers = set([])
        all_reverse_primers = set([])
        if suppress_sampleid_check:
            all_primers.update(
             [upper(primer).strip() for primer in (line[primer_ix].split(','))])
            if reverse_primers != "disable":
                all_rev_primers.update(
                 [DNA.rc(upper(primer)).strip()
                 for primer in (line[rev_primer_ix].split(','))])
            continue
        if reverse_primers != "disable":
            all_raw_reverse_primers = set(
             [DNA.rc(upper(primer)).strip() for 
             primer in line[rev_primer_ix].split(',')])
            for curr_primer in all_raw_reverse_primers:
                all_reverse_primers.update(expand_degeneracies([curr_primer]))
            ids_rev_primers[line[sample_id_ix]] = all_reverse_primers
        
        all_raw_forward_primers =\
         set([upper(primer).strip() for primer in line[primer_ix].split(',')])
        for curr_primer in all_raw_forward_primers:
            all_forward_primers.update(expand_degeneracies([curr_primer]))
        ids_primers[line[sample_id_ix]] = all_forward_primers
    
    if suppress_sampleid_check:
        final_all_primers = set([])
        final_all_rev_primers = set([])
        for curr_primer in set(all_primers):
            final_all_primers.update(expand_degeneracies([curr_primer]))
        for curr_primer in set(all_rev_primers):
            final_all_rev_primers.update(expand_degeneracies([curr_primer]))
        ids_primers["all_primers"] = final_all_primers
        ids_rev_primers["rev_primers"] = final_all_rev_primers
        
    return ids_primers, ids_rev_primers
    
# End mapping data processing ****

# Log formatting functions **********

def format_detailed_data(detailed_quality_data):
    """ Formats detailed quality data for logging
    
    detailed_quality_data: dictionary of fasta labels:{filters} with filters
     set to 1 if particular test failed or if check was positive for detection, 
     such as detecting a low quality window
    """
    
    formatted_detailed_data = ""
    # Variable number of tests, so no static header, generate header from
    # tests recorded.
    
    initial_line = "# Each line contains the fasta label followed by a 1 or 0 for the results of each filter/check\n"
    
    filters = detailed_quality_data['all_tests']

    header = "Sequence label\t%s\n" % "\t".join(filters)
    
    formatted_detailed_data += initial_line
    formatted_detailed_data += header
    
    # Step through this way so original order of sequence maintained
    for curr_label in detailed_quality_data['seq_order']:
        curr_line = curr_label
        for curr_filter in filters:
            try:
                curr_line += "\t%s" %\
                 detailed_quality_data[curr_label][curr_filter]
            except KeyError:
                curr_line += "\tNA"
        curr_line += "\n"
        formatted_detailed_data += curr_line
            
    return formatted_detailed_data
    
def format_log_data(final_log_data,
                    seq_length_data,
                    min_seq_len,
                    max_seq_len,
                    min_qual_score,
                    retain_primer,
                    max_ambig,
                    suppress_ambig_check,
                    max_homopolymer,
                    suppress_homopolymer_check,
                    max_primer_mismatch,
                    output_dir,
                    qual_score_window,
                    discard_bad_windows,
                    suppress_primer_check,
                    reverse_primers,
                    reverse_primer_mismatches,
                    record_qual_scores,
                    truncate_ambi_bases,
                    suppress_sampleid_check,
                    enable_all_checks):
    """ Formats summary log data

    final_log_data: Dictionary of counts of failures
    seq_length_data: 2d list of raw/processed sequence lengths
    min_seq_len: Minimum sequence length
    max_seq_len: Max sequence length
    min_qual_score: Minimum overall qual score, also used for sliding window
    retain_primer: If True, forward primer is not removed from output
    max_ambig: Max ambiguous "N" bases allowed
    suppress_ambig_check: If True, disables ambiguous base check
    max_homopolymer: Max homopolymer length allowed
    suppress_homopolymer_check: If True, disabled homopolymer check
    max_primer_mismatch: Maximum allowed forward primer mismatches
    output_dir: output directory
    qual_score_window: Sliding qual window check, 0 disables, any other size
     is the size of the window
    discard_bad_windows: If True, and low qual sliding window detected, 
     discards the sequence rather than writing
    suppress_primer_check: If True, suppresses forward primer mismatch and
     removal check
    reverse_primers: Can be 'disable', 'truncate_only', and 'truncate_remove'.
     Uses ReversePrimer from mapping data to do local alignment of reverse
     primer, and will truncate at that point in the sequence. If 
     truncate_remove option is enabled, will discard reads where reverse
     primer can not be found.
    reverse_primer_mismatches: Allowed mismatches in reverse primer.
    record_qual_scores: If True, will write matching quality score sequence to
     trimmed fasta file.
    truncate_ambi_bases: If True, will ignore max_homopolymer value, and instead
     will truncate at first ambiguous "N" base encountered.
    suppress_sampleid_check: If True, will not test for unique SampleIDs in 
     the mapping file, and will use all forward and (if present) reverse 
     primers in the mapping file during tests. 
    enable_all_checks: If True, will do all quality tests rather than halting
     on the first failure, and will write detailed per-sequence report on
     failures.
    """
    
    formatted_log_data = ["# Log details for quality filtering\n\n"]
    
    fastas = ",".join(list(final_log_data['fasta_files']))
    formatted_log_data.append("Input fasta file(s):\t%s\n" % (fastas))

    if final_log_data['qual_files']:
        quals = ",".join(list(final_log_data['qual_files']))
    else:
        quals = "NA"
    formatted_log_data.append("Input qual file(s):\t%s\n" % (quals))
    
    formatted_log_data.append("Output directory:\t%s\n\n" % output_dir)

    formatted_log_data.append("Parameter Settings\n\n")
    formatted_log_data.append("Minimum sequence lengths:\t%d\n" % min_seq_len)
    formatted_log_data.append("Maximum sequence lengths:\t%d\n" % max_seq_len)
    formatted_log_data.append("Minimum quality score:\t%d\n" % min_qual_score)
    formatted_log_data.append("Retain primer sequence:\t%s\n" % retain_primer)
    formatted_log_data.append("Maximum ambiguous bases allowed:\t%d\n"
     % max_ambig)
    formatted_log_data.append("Suppress ambiguous base check:\t%s\n" %\
     suppress_ambig_check)
    formatted_log_data.append("Truncate at first ambiguous base:\t%s\n" %\
     truncate_ambi_bases)
    formatted_log_data.append("Maximum homopolymer length allowed:\t%d\n" %\
     max_homopolymer)
    formatted_log_data.append("Suppress homopolymer check:\t%s\n" %\
     suppress_homopolymer_check)
    formatted_log_data.append("Suppress primer check:\t%s\n" %\
     suppress_primer_check)
    formatted_log_data.append("Maximum allowed primer mismatches:\t%d\n" %\
     max_primer_mismatch)
    formatted_log_data.append("Sliding quality window check setting:\t%d\n" %\
     qual_score_window)
    formatted_log_data.append("Discard sequences with low quality windows:\t%s\n" %\
     discard_bad_windows)
    formatted_log_data.append("Reverse primer settings:\t%s\n" %\
     reverse_primers)
    formatted_log_data.append("Reverse primer maximum mismatches:\t%d\n" %\
     reverse_primer_mismatches)
    formatted_log_data.append("Suppress SampleID matching:\t%s\n" %\
     suppress_sampleid_check)
    formatted_log_data.append("Enabled detailed logging of all filters:\t%s\n" %\
     enable_all_checks)
    formatted_log_data.append("Record matching quality scores:\t%s\n\n" %\
     record_qual_scores)
     
    formatted_log_data.append("Input sequence count:\t%d\n" %\
     final_log_data['seq_counts'])
    formatted_log_data.append("Seqs written:\t%d\n" %\
     final_log_data['seqs_written'])
    formatted_log_data.append("Percent of seqs written\t%3.2f\n" %\
     (int(final_log_data['seqs_written'])/int(final_log_data['seq_counts'])))
    # Handle case of no sequence data
    try:
        formatted_log_data.append("Raw min/max/mean sequence lengths:"+\
         "\t%3.2f/%3.2f/%3.2f\n" % (min(seq_length_data['raw']),
         max(seq_length_data['raw']), mean(seq_length_data['raw'])))
    except ValueError:
        formatted_log_data.append("Raw min/max/mean sequence lengths: NA\n")
    try:
        formatted_log_data.append("Processed min/max/mean sequence lengths:"+\
         "\t%3.2f/%3.2f/%3.2f\n\n" % (min(seq_length_data['processed']),
         max(seq_length_data['processed']), mean(seq_length_data['processed'])))
    except ValueError:
        formatted_log_data.append("Processed min/max/mean sequence lengths: NA\n")
     
    formatted_log_data.append("Filters that resulted in discarded reads\n")
    formatted_log_data.append("Sequences below minimum length %d:\t%d\n" %\
     (min_seq_len, final_log_data['below_min_seq_len']))
    formatted_log_data.append("Sequences above maximum length %d:\t%d\n" %\
     (max_seq_len, final_log_data['max_seq_len']))
    formatted_log_data.append("Sequences below minimum average quality score %d:\t%d\n" %\
     (min_qual_score, final_log_data['min_ave_qual_score']))
    formatted_log_data.append("Sequences with homopolymers longer than %d:\t%d\n" %\
     (max_homopolymer, final_log_data['max_homopolymer']))
    formatted_log_data.append("Sequences exceeding %d ambiguous bases:\t%d\n" %\
     (max_ambig, final_log_data['max_ambig']))
    formatted_log_data.append("Exceeded maximum primer mismatch of %d:\t%d\n" %\
     (max_primer_mismatch, final_log_data['max_primer_mismatch']))
    formatted_log_data.append("Sequences too short after reverse primer truncation:\t%d\n" %\
     final_log_data['too_short_after_revprimer_truncation'])
    formatted_log_data.append("Sequences too short after sliding quality window truncation:\t%d\n" %\
     final_log_data['too_short_after_window_truncation'])
    formatted_log_data.append("Sequences too short after ambiguous base truncation:\t%d\n" %\
     final_log_data['too_short_after_ambig_truncation'])
    formatted_log_data.append("Sequences discarded for low quality window detection:\t%d\n" %\
     final_log_data['low_qual_window_discarded'])
    formatted_log_data.append("Sequences discarded for failing to detect reverse primer:\t%d\n\n" %\
     final_log_data['seqs_discarded_no_rev_primer'])
      
    formatted_log_data.append("Other sequence details:\n\n")
    formatted_log_data.append("Sequences with detected reverse primer:\t%d\n" %\
     final_log_data['rev_primers_found'])
    formatted_log_data.append("Sequences with detected low quality window:\t%d\n" %\
     final_log_data['low_qual_window_found'])
    formatted_log_data.append("Sequences with detected ambiguous base:\t%d\n\n" %\
     final_log_data['ambig_base_truncation'])
      
    if final_log_data['seq_ids_not_in_mapping']:
        formatted_log_data.append("Fasta labels not matching SampleIDs in mapping file:\n")
        for curr_id in final_log_data['seq_ids_not_in_mapping']:
            formatted_log_data.append("%s\n" % curr_id)

    combined_log_data = "".join(formatted_log_data)
    
    return combined_log_data
# End log formatting functions *******
    
    
# File IO functions *****************


def write_qual_line(demultiplexed_qual_f,
                    label_line,
                    qual_seq,
                    qual_line_size = 60):
    """ Writes quality score sequence out in proper format
    
    demultiplexed_qual_f:  open file object to write label/qual scores to.
    label_line:  qual label to write
    qual_seq:  current qual sequence, list of scores.
    """
    current_qual_scores_lines = []  
    # Quality score format is a string of 60 base calls, followed by a 
    # newline, until the last N bases are written
    for slice in range(0, len(qual_seq), qual_line_size):

        current_segment = qual_seq[slice:slice + qual_line_size]
        current_qual_scores_lines.append(" ".join(map(str, current_segment)))
    
    demultiplexed_qual_f.write(">%s\n" % label_line)
    demultiplexed_qual_f.write('\n'.join(current_qual_scores_lines))
    demultiplexed_qual_f.write('\n')
    
# End IO functions *************

# Misc functions ***********

def count_mismatches(seq, primer, max_mm):
    """Counts mismatches, primer should be <= length of the seq.
    
    seq: fasta sequence, string of nucleotides
    primer: primer sequence, string of nucleotides
    max_mm: Maximum allowed mismatches
    """
    mm = 0
    for i in range(min(len(primer), len(seq))):
        if seq[i] != primer[i]:
            mm += 1
            if mm > max_mm:
                return True
    return False

def expand_degeneracies(raw_primers):
    """ Returns all non-degenerate versions of a given primer sequence """
    
    expanded_primers=[]
    
    for raw_primer in raw_primers:
        primers=SequenceGenerator(template=raw_primer.strip(),
         alphabet=IUPAC_DNA)
        
        for primer in primers:
            expanded_primers.append(primer)
        
    return expanded_primers

def make_histograms(raw_lengths, post_lengths, binwidth=10):
    """Makes histogram data for raw and post lengths"""
    if post_lengths:
        min_len = min([min(post_lengths), min(raw_lengths)])
    else:
        min_len = min(raw_lengths)
    max_len = max(raw_lengths)
    hist_floor = int((min_len/binwidth)*binwidth)
    ceil = int(((max_len/binwidth)+2)*binwidth)
    bins = arange(hist_floor, ceil, binwidth)
    raw_hist = histogram(raw_lengths,bins)[0]
    post_hist, bin_edges = histogram(post_lengths,bins)
    return raw_hist, post_hist, bin_edges
    
def count_ambig(curr_seq, ambig_chars="RYMKWSBDHVN"):
    """Counts non-ambiguous characters in seq
    
    curr_seq: fasta sequence
    ambig_chars: nucleotide characters"""
    total = 0
    for ambig_char in ambig_chars:
        total += curr_seq.count(ambig_char)
    return total

    
def pair_hmm_align_unaligned_seqs(seqs,moltype=DNA,params={}):
    """
        Checks parameters for pairwise alignment, returns alignment.
        
        Code from Greg Caporaso.
    """
    
    seqs = LoadSeqs(data=seqs,moltype=moltype,aligned=False)
    try:
        s1, s2 = seqs.values()
    except ValueError:
        raise ValueError,\
         "Pairwise aligning of seqs requires exactly two seqs."
    
    try:
        gap_open = params['gap_open']
    except KeyError:
        gap_open = 5
    try:
        gap_extend = params['gap_extend']
    except KeyError:
        gap_extend = 2
    try:
        score_matrix = params['score_matrix']
    except KeyError:
        score_matrix = make_dna_scoring_dict(\
         match=1,transition=-1,transversion=-1)
    
    return local_pairwise(s1,s2,score_matrix,gap_open,gap_extend)
    
def local_align_primer_seq(primer,sequence):
    """Perform local alignment of primer and sequence
    
        primer: Input primer sequence
        sequence: target sequence to test primer against
        
        Returns the number of mismatches,
         and the start position in sequence of the hit.
         
        Modified from code written by Greg Caporaso.
    """
    # The scoring function which can be passed to cogent.alignment.algorithms.sw_align
    
    query_primer = primer

    query_sequence = str(sequence)
     
    # Get alignment object from primer, target sequence
    alignment = pair_hmm_align_unaligned_seqs([query_primer,query_sequence])
    
    # Extract sequence of primer, target site, may have gaps in insertions
    # or deletions have occurred.
    primer_hit = str(alignment.Seqs[0])
    target_hit = str(alignment.Seqs[1])
    
    # Count insertions and deletions
    insertions = primer_hit.count('-')
    deletions = target_hit.count('-')
    
    mismatches = 0
    for i in range(len(target_hit)):
        if target_hit[i] != primer_hit[i] and\
         target_hit[i] != '-' and primer_hit[i] != '-': 
            mismatches += 1
    try:
        hit_start = query_sequence.index(target_hit.replace('-',''))
    except ValueError:
        raise ValueError,('substring not found, query string '
         '%s, target_hit %s' % (query_sequence, target_hit))
    
    # sum total mismatches
    mismatch_count = insertions + deletions + mismatches
    
    return mismatch_count, hit_start
    
def check_window_qual_scores(qual_scores, window=50, min_average=25):
    """Check that all windows have ave qual score > threshold.
    
    qual_scores: list of quality scores
    window: size of window in base pairs
    min_average: minimum average quality score for window
    
    If does not find a low quality window, will return False, last index of
    window, if a low quality window is found, this will return True, first
    index of the low quality window."""

    
    # Code from Jens Reeder, added 1-13-2010
    l = len(qual_scores)

    window = min(window, l)
    if (window == 0):
        return True
    #initialize with sum of first window
    window_score = sum(qual_scores[:window])   
    idx = 0
    while (window_score/float(window) >= min_average
           and idx < l-window):
            #'Move' window
            window_score += qual_scores[idx+window] - qual_scores[idx]
            idx += 1
    if (idx == l-window):
        # we processed all qual_scores, must be good
        # Return index for truncation purposes
        return False, idx
    else:
        return True, idx
        
def seq_exceeds_homopolymers(curr_seq, max_len=6):
    """Returns False if primer contains any homopolymer > allowed length"""
    for base in 'ATGC':
        curr = base * (max_len+1)
        if curr in curr_seq:
            return True
    return False

# End misc functions **********