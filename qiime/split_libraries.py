#!/usr/bin/env python
#file split_libraries.py
"""Performs preprocessing steps for barcoded library analysis, e.g. 454.

Specifically, does the quality-filtering step (using several criteria) and
renames each read with the appropriate library id.
"""
__author__ = "Rob Knight and Micah Hamady"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight", "Micah Hamady", "Greg Caporaso", "Kyle Bittinger"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

import re
from cogent.parse.fasta import FastaFinder, MinimalFastaParser
from numpy import array, mean, arange, histogram
from numpy import __version__ as numpy_version
from qiime.check_id_map import parse_id_map
from qiime.barcode import correct_barcode
from gzip import GzipFile
from optparse import OptionParser
from os import mkdir, stat
from collections import defaultdict
from hamming import decode_barcode_8
from golay import decode as decode_golay_12

## Including new=True in the histogram() call is necessary to 
## get the correct result in versions prior to NumPy 1.2.0,
## but the new keyword will be removed in NumPy 1.4. In NumPy 1.2.0 
## or later, new=True raises a Warning regarding 
## deprecation of new. One solution to this is to check the 
## numpy version, and if it's less than 1.2.0, overwrite histogram 
## with new=True as the default. This avoids a deprecation warning message 
## in versions 1.2.0 through 1.3.*, and a try/except to handle errors from 
## versions 1.4.0 or later. 
numpy_version = re.split("[^\d]", numpy_version)
numpy_version = tuple([int(i) for i in numpy_version if i.isdigit()])
if numpy_version < (1,3,0):
    numpy_histogram = histogram
    def histogram(a, bins=10, range=None, normed=False, weights=None):
        return numpy_histogram(a,bins=bins,range=range,\
         normed=normed,weights=weights,new=True)

IUPAC_DNA = {'U':'T','T':'T','C':'C','A':'A','G':'G', 'R':'AG','Y':'TC',
             'W':'TA','S':'CG','M':'CA','K':'TG','B':'TCG','D':'TAG','H':'TCA',
             'V':'CAG', 'N':'TCAG' }
     
STANDARD_BACTERIAL_PRIMER = "CATGCTGCCTCCCGTAGGAGT" #R357
STANDARD_FUNGAL_PRIMER = "ACTCTGGACCTGGTGAGTTTC" #need region info from Noah

# Supported barcode types - need to adapt these functions to ignore list
# of valid barcodes that the generic decoder requires
BARCODE_TYPES = { "golay_12":(12, lambda bc, bcodes: decode_golay_12(bc)), 
        "hamming_8":(8, lambda bc, bcodes: decode_barcode_8(bc))}

def generate_possibilities(seq_str):
    """Generate non-degenerate patterns for given degenerate pattern"""
    expansions = [IUPAC_DNA[x] for x in seq_str]
    all = []
    for expan in expansions:
        if not expan: return []
        if all:
            all = [i + j for i in all for j in expan]
        else:
            all = [j for j in expan]
    return all

def get_infile(filename):
    """Returns filehandle, allowing gzip input."""
    if filename.endswith(".gz"):
        fin = GzipFile(filename, "rb")     
    else:
        fin = open(filename, "U") 
    return fin

def count_mismatches(seq1, seq2, max_mm):
    """Counts mismatches, assumes lengths equal.
    """
    mm = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mm += 1
            if mm > max_mm:
                return mm
    return mm

def ok_mm_primer(primer_seq, all_primer_seqs, primer_mm):
    """Check if primer_seq matches any primer within max mismatches.
    
    TODO: if we start using highly degenerate primers, should refactor using
    faster algorithm.
    """ 
    for curr_pat in all_primer_seqs:
        curr_mm = count_mismatches(primer_seq, curr_pat, primer_mm)
        if count_mismatches(primer_seq, curr_pat, primer_mm) <= primer_mm:
            return True
    return False

def check_map(infile, has_barcodes=True):
    """Check mapping file and extract list of valid barcodes """
    hds, id_map, dsp, run_description, errors, warnings = \
        parse_id_map(infile, has_barcodes)
    barcode_to_sample_id = {} 
    for sample_id, sample in id_map.items():
        barcode_to_sample_id[sample['BarcodeSequence']] = sample_id
    return hds, id_map, barcode_to_sample_id, warnings, errors

def get_primer_seqs(primer_seq_pats, min_primer_len=5):
    """Load primer sequences from primer patterns.
    
    Question: should it be possible to add support for unequal length primers?
    """
    primer_seq_len = None 
    all_primer_seqs = set([])
    for curr_pat in primer_seq_pats:
        curr_pat_len = len(curr_pat)
        if curr_pat_len < min_primer_len:
            raise ValueError, "Primer sequence pattern must be >=%s nt" % \
                min_primer_len
        if primer_seq_len is None:
            primer_seq_len = curr_pat_len 
        elif curr_pat_len != primer_seq_len:
            raise ValueError, "Primer seq is wrong length"
        for curr_seq in generate_possibilities(curr_pat):
            all_primer_seqs.add(curr_seq)
    return all_primer_seqs, primer_seq_len

def fasta_ids(fasta_files, verbose=False):
    """ Returns list of ids in FASTA files """
    all_ids = set([])
    for fasta_in in fasta_files:
        for label, seq in MinimalFastaParser(fasta_in):
            rid = label.split()[0]
            if rid in all_ids:
                raise ValueError, \
                    "Duplicate ID found in FASTA/qual file: %s" % label      
            all_ids.add(rid)
    return all_ids

def qual_score(infile):
    """Load quality scores."""
    id_to_qual = {}
    for rec in FastaFinder(infile):
        curr_id = rec[0][1:]
        curr_qual = ' '.join(rec[1:])
        parts = array(map(int, curr_qual.split()))
        curr_pid = curr_id.split()[0]
        id_to_qual[curr_pid] = parts
    return id_to_qual

def qual_scores(qual_files):
    """ Load qual scores into dict of {id:qual_scores}.
    
    No filtering is performed at this step.
    """
    qual_mappings = {}
    for qual_file in qual_files:
        qual_mappings.update(qual_score(qual_file))
    return qual_mappings

def count_ambig(curr_seq, valid_chars='ATCG'):
    """Counts non-standard characters in seq"""
    up_seq  = curr_seq.upper()
    total = 0
    for vchar in valid_chars:
        total += up_seq.count(vchar)
    return len(curr_seq) - total

def split_seq(curr_seq, barcode_len, primer_seq_len):
    """Split sequence into parts barcode, primer and remainder"""
    curr_barcode = curr_seq[0:barcode_len] 
    rest_of_seq = curr_seq[barcode_len:]
    primer_seq = rest_of_seq[0:primer_seq_len]
    rest_of_seq = rest_of_seq[primer_seq_len:] 
    return curr_barcode, primer_seq, rest_of_seq 

def primer_exceeds_mismatches(primer_seq, all_primer_seqs, max_primer_mm):
    """Returns True if primer exceeds allowed mismatches"""   
    if primer_seq not in all_primer_seqs:
        if not ok_mm_primer(primer_seq, all_primer_seqs, max_primer_mm):
            return True
    return False

def seq_exceeds_homopolymers(curr_seq, max_len=6):
    """Returns False if primer contains any homopolymer > allowed length"""
    for base in 'ATGC':
        curr = base * (max_len+1)
        if curr in curr_seq:
            return True
    return False

def check_barcode(curr_barcode, barcode_type, valid_map):
    """Return whether barcode is valid, and attempt at correction."""
    if curr_barcode in valid_map:
        return False, curr_barcode
    if barcode_type in BARCODE_TYPES:
        expect_len, curr_bc_fun  = BARCODE_TYPES[barcode_type]
        barcode, num_errors = curr_bc_fun(curr_barcode, valid_map)
        return num_errors, barcode
    else:
        return True, None

def make_histograms(pre_lengths, post_lengths, binwidth=10):
    """Makes histogram data for pre and post lengths"""
    min_len = min(pre_lengths)
    max_len = max(pre_lengths)
    floor = (min_len/binwidth)*binwidth
    ceil = ((max_len/binwidth)+2)*binwidth
    bins = arange(floor, ceil, binwidth)
    pre_hist = histogram(pre_lengths,bins)[0]
    post_hist, bin_edges = histogram(post_lengths,bins)
    return pre_hist, post_hist, bin_edges

def format_histograms(pre_hist, post_hist, bin_edges):
    """Returns text-formatted histogram."""
    lines = []
    lines.append('Length\tBefore\tAfter')
    for edge, pre, post in zip(bin_edges, pre_hist, post_hist):
        lines.append('\t'.join(map(str, [edge, pre, post])))
    return '\n'.join(lines)

class SeqQualBad(object):
    """Checks if a seq and qual score are bad, saving ids that are bad."""

    def __init__(self, name, f):
        """New SeqQualBad keeps track of failed ids."""
        self.FailedIds = []
        self.Name = name
        self.F = f

    def __call__(self, id_, seq, qual):
        """SeqQualBad called on id, seq and qual returns bool.

        Note: saves failed ids in self.FailedIds."""
        result = self.F(id_, seq, qual)
        if result:
            self.FailedIds.append(id_)
        return result

    def __str__(self):
        """SeqQualBad str returns tab-delimited output of counts."""
        return "%s\t%s" % (self.Name, len(self.FailedIds))


def qual_missing(id_, seq, qual):
    """returns True if qual is None"""
    return qual is None

QualMissing = SeqQualBad('Missing Qual Score', qual_missing)

def get_seq_lengths(seq_lengths, bc_counts):
    """Convenience wrapper for getting lengths of good and bad seqs"""
    all_seq_lengths = seq_lengths.values()
    all_seq_ids = set(seq_lengths.keys())
    bad_seq_ids = set(bc_counts[None]).union(set(bc_counts['#FAILED']))
    good_seq_ids = all_seq_ids - bad_seq_ids
    good_seq_lengths = map(seq_lengths.__getitem__, good_seq_ids)
    return all_seq_lengths, good_seq_lengths
 
def check_seqs(fasta_out, fasta_files, starting_ix, valid_map, qual_mappings, 
    filters, barcode_len, primer_seq_len, keep_primer, barcode_type, 
    max_bc_errors):
    """Checks fasta-format sequences and qual files for validity."""
    seq_lengths = {}
    bc_counts = defaultdict(list)
    curr_ix = starting_ix
    corr_ct = 0 #count of corrected barcodes
    
    for fasta_in in fasta_files:
        for curr_id, curr_seq in MinimalFastaParser(fasta_in):
            curr_rid = curr_id.split()[0]
            curr_len = len(curr_seq)
            curr_qual = qual_mappings.get(curr_rid, None)
            seq_lengths[curr_rid] = curr_len
            failed = False
            for f in filters:
                failed = failed or f(curr_rid, curr_seq, curr_qual)
            if failed:  #if we failed any of the checks, bail out here
                bc_counts['#FAILED'].append(curr_rid)
                continue
            # split seqs
            cbc, cpr, cres = split_seq(curr_seq, barcode_len, primer_seq_len)
            # get current barcode
            try:
                bc_diffs, curr_bc = \
                    check_barcode(cbc, barcode_type, valid_map.keys())
                if bc_diffs > max_bc_errors:
                    raise ValueError, "Too many errors in barcode"
                corr_ct += bool(bc_diffs)
            except Exception, e:
                bc_counts[None].append(curr_rid)
                continue
            bc_counts[curr_bc].append(curr_rid)
            curr_samp_id = valid_map.get(curr_bc, 'Unassigned')
            new_id = "%s_%d" % (curr_samp_id, curr_ix)
            # check if writing out primer
            write_seq = cres
            if keep_primer:
                write_seq = cpr + cres
            fasta_out.write(">%s %s orig_bc=%s new_bc=%s bc_diffs=%s\n%s\n" % 
                (new_id, curr_rid, cbc, curr_bc, int(bc_diffs), write_seq))
            curr_ix += 1
    log_out = format_log(bc_counts, corr_ct, seq_lengths, valid_map, filters)
    all_seq_lengths, good_seq_lengths = get_seq_lengths(seq_lengths, bc_counts)
    return log_out, all_seq_lengths, good_seq_lengths

def format_log(bc_counts, corr_ct, seq_lengths, valid_map, filters):
    """Makes log lines"""
    log_out = []
    all_seq_lengths, good_seq_lengths = get_seq_lengths(seq_lengths, bc_counts)
    log_out.append("Number raw seqs\t%d" % len(seq_lengths)) 
    log_out.append("Number seqs written\t%d" % len(good_seq_lengths)) 
    for f in filters:
        log_out.append(str(f))
    log_out.append("Raw len min/max/avg\t%.1f/%.1f/%.1f" % 
        (min(all_seq_lengths), max(all_seq_lengths), mean(all_seq_lengths)))
    
    if good_seq_lengths:
        log_out.append("Wrote len min/max/avg\t%.1f/%.1f/%.1f" % 
        (min(good_seq_lengths), max(good_seq_lengths), mean(good_seq_lengths))) 
    
    #figure out which barcodes we got that didn't come from valid samples
    valid_bc_nomap = set(bc_counts) - set(valid_map) - set([None,'#FAILED'])
    valid_bc_nomap_counts = [(len(bc_counts[b]),b) for b in valid_bc_nomap]

    log_out.append("BC ok, not in mapping\t%d" % len(valid_bc_nomap_counts)) 
    log_out.append("BC corrected/not\t%d/%d" % 
            (corr_ct, len(bc_counts[None]))) 
    if valid_bc_nomap:
        log_out.append("BC ok (missed)\tCount")
        for count, bc in reversed(sorted(valid_bc_nomap_counts)):
            log_out.append("%s\t%d" % (bc, count)) 
   
    sample_cts = [(len(bc_counts[bc]), bc, sample_id) for bc, sample_id 
            in valid_map.items()]
    if sample_cts:
        log_out.append("Sample\tSequence Count\tBarcode")
        for count, bc, sample_id in reversed(sorted(sample_cts)):
            log_out.append("%s\t%d\t%s" % (sample_id, count, bc))

        filtered_sample_cts = [s[0] for s in sample_cts if s[0]]
        if filtered_sample_cts:
            log_out.append("Num Samples\t%d" % len(filtered_sample_cts))
            log_out.append("Sample ct min/max/mean: %d / %d / %.2f" % (
            min(filtered_sample_cts), max(filtered_sample_cts), 
            mean(filtered_sample_cts)))
    return log_out

def preprocess(fasta_files, qual_files, mapping_file, 
    primer_seq_pats=STANDARD_BACTERIAL_PRIMER, 
    barcode_type="golay_12",
    min_seq_len=200, max_seq_len=1000, min_qual_score=25, starting_ix=1,
    keep_primer=True, max_ambig=0, max_primer_mm=1, trim_seq_len=True,
    dir_prefix='.', max_bc_errors=2, max_homopolymer=4):
    """
    Preprocess barcoded libraries, e.g. from 454.

    Parameters:

    fasta_files: list of raw 454 fasta files, fasta format.
    
    qual_files: list of raw 454 qual file(s)
    
    mapping_file: mapping file with BarcodeSequence column containing valid 
    barcodes used in the 454 run 

    primer_seq_pats: list of valid primer sequences, can be degenerate.

    barcode_type: type of barcode, e.g. golay_12. Should appear in list of
    known barcode types.

    min_seq_len: minimum sequence length to allow.

    max_seq_len: maximum sequence length to allow.

    min_qual_score: minimum average qual score considered acceptaable.
    
    starting_ix: integer to start sample sequence numbering at.

    keep_primer: when True, will keep primer sequence, otherwise will strip it 

    max_ambig: maximum number of ambiguous bases to allow in the read.

    max_primer_mm: maximum number of primer mismatches to allow.

    trim_seq_len: if True (default), calculates lengths after trimming.

    dir_prefix: prefix of directories to write files into.

    max_bc_errors: maximum number of barcode errors to allow in output seqs

    Result:
    in dir_prefix, writes the following files:
    id_map.xls: 2-column tab-delimited text format orig_id:new_id
    error_map.xls: 2-column tab-delimited text format orig_id:fail_reasons
    seqs.fasta: sequences with new ids lib_index in fasta format
    lengths.xls: histograms of unfiltered and filtered lengths, resolution 10 bp
    """
    if max_seq_len < 10:
        raise ValueError, "Max sequence must be >= 10"
    if min_seq_len >= max_seq_len:
        raise ValueError, "Min len cannot be >= max len"
    if min_qual_score < 0:
        raise ValueError, "Min qual score must be > 0"
    if starting_ix < 1:
        raise ValueError, "Starting index must be > 0."
    if max_ambig < 0:
        raise ValueError, "Max ambig chars must be >= 0."
    if max_primer_mm < 0:
        raise ValueError, "Max primer mismatches must be >= 0."
    if min_qual_score < 5:
        raise ValueError, "Min qual score must be >= 0."

    try:
        stat(dir_prefix)
    except OSError:
        mkdir(dir_prefix)

    # Generate primer sequence patterns
    all_primer_seqs, primer_seq_len = \
        get_primer_seqs(primer_seq_pats.split(','))

    # Check mapping file and get barcode mapping 
    map_file = open(mapping_file, 'U')
    headers, id_map, valid_map, warnings, errors = \
        check_map(map_file)
    map_file.close()
    if errors:
        raise ValueError, "Invalid mapping file. Validate with MapCheck first: %s" % "\n".join(errors)

    # Check barcode type
    if barcode_type not in BARCODE_TYPES:
        try:
            barcode_len, barcode_fun = int(barcode_type), correct_barcode
        except ValueError:
            raise ValueError, "Unsupported barcode type: %s" % barcode_type
    else:
        barcode_len, barcode_fun = BARCODE_TYPES[barcode_type]

    fasta_files = map(get_infile, fasta_files)
    qual_files = map(get_infile, qual_files)
    
    # Check fasta files valid format, no duplicate ids
    # and ids match between fasta and qual files
    all_fasta_ids = fasta_ids(fasta_files)
    all_qual_ids = fasta_ids(qual_files)
    if qual_files and (len(all_fasta_ids) != len(all_qual_ids)):
        f_ids = all_fasta_ids.difference(all_qual_ids)
        q_ids = all_qual_ids.difference(all_fasta_ids)
        raise ValueError, "Found %d ids in fasta file not in qual file, %d ids in qual file not in fasta"  % (len(f_ids), len(q_ids))

    for f in fasta_files:
        f.seek(0)
    if qual_files:
        for q in qual_files:
            q.seek(0)
        # Load quality scores 
        qual_mappings = qual_scores(qual_files)
        for q in qual_files:
            q.close()
    else:
        qual_mappings = {}

    #make filters
    filters = []
    if qual_mappings:
        filters.append(QualMissing)
        filters.append(SeqQualBad(
            'Mean qual score below minimum of %s' % min_qual_score, 
            lambda id_, seq, qual: mean(qual) < min_qual_score))
    #seq len filter depends on whether we're including the barcode
    if trim_seq_len:
        trim = barcode_len + primer_seq_len
        filters.append(SeqQualBad(
            'Length outside bounds of %s and %s' % (min_seq_len,max_seq_len),
            lambda id_, seq, qual: \
                not (min_seq_len<=len(seq)-trim<= max_seq_len)))
    else:
        filters.append(SeqQualBad(
            'Length outside bounds of %s and %s' % (min_seq_len,max_seq_len),
            lambda id_, seq, qual: not (min_seq_len<=len(seq)<= max_seq_len)))
    filters.append(SeqQualBad(
        'Num ambiguous bases exceeds limit of %s' % max_ambig,
        lambda id_, seq, qual: count_ambig(seq) > max_ambig))
    filters.append(SeqQualBad(
        'Num mismatches in primer exceeds limit of %s' % max_primer_mm,
        lambda id_, seq, qual: primer_exceeds_mismatches(
            seq[barcode_len:barcode_len+primer_seq_len], all_primer_seqs,
            max_primer_mm)))
    filters.append(SeqQualBad(
        'Max homopolymer run exceeds limit of %s' % max_homopolymer,
        lambda id_, seq, qual: seq_exceeds_homopolymers(
            seq[barcode_len+primer_seq_len:], max_homopolymer)))

    # Check seqs and write out
    fasta_out = open(dir_prefix + '/' + 'seqs.fna', 'w+')
    log_stats, pre_lens, post_lens = check_seqs(fasta_out, fasta_files, 
        starting_ix, valid_map, qual_mappings, filters, barcode_len,
        primer_seq_len, keep_primer, barcode_type, max_bc_errors)

    # Write log file
    log_file = open(dir_prefix + '/' + "split_library_log.txt", 'w+')
    log_file.write('\n'.join(log_stats))
    log_file.close()

    # Write sequence distros here
    histogram_file = open(dir_prefix + '/' + 'histograms.txt', 'w+')
    histogram_file.write(format_histograms
        (*make_histograms(pre_lens, post_lens)))
    histogram_file.close()


usage_str = \
"""usage: %prog [options] {-i FASTA_FNAMES -q QUAL_FNAMES -m MAP_FNAME}

[] indicates optional input (order unimportant) 
{} indicates required input (order unimportant) 

FASTA_FNAMES: Comma-delimited paths of sequence files, in FASTA format.
QUAL_FNAMES: Comma-delimited paths of quality files, in FASTA-like
  format.
MAP_FNAME: Path to tab-delimited mapping file.  Must contain a header
  line indicating SampleID in the first column and BarcodeSequence in
  the second.

Example usage:

Process sequences from two files (a.fna, b.fna) using the sample-to-
barcode associations in samples.txt.  The quality files for the
sequences (a.qual, b.qual) must also be provided.

python %prog -i a.fna,b.fna -q a.qual,b.qual -m samples.txt
"""


def make_cmd_parser():
    """Returns the command-line options."""
    parser = OptionParser()
    parser.usage = usage_str
    parser.version = 'Version: %prog ' +  __version__
    parser.add_option('-m', '--map', dest='map_fname', 
        help='name of mapping file [REQUIRED]')
    parser.add_option('-f', '--fasta', dest='fasta_fnames', 
        help='names of fasta files, comma-delimited [REQUIRED]')
    parser.add_option('-q', '--qual', dest='qual_fnames', 
        help='names of qual files, comma-delimited [REQUIRED]')
    parser.add_option('-p', '--primers', default=STANDARD_BACTERIAL_PRIMER,
        help='degen sequences of primers, comma-delimited [default: %default]')
    parser.add_option('-l', '--min-seq-length', dest='min_seq_len',
        type=int, default=200,
        help='minimum sequence length, in nucleotides [default: %default]')
    parser.add_option('-L', '--max-seq-length', dest='max_seq_len',
        type=int, default=1000,
        help='maximum sequence length, in nucleotides [default: %default]')
    parser.add_option('-t', '--trim-seq-length', dest='trim_seq_len',
        action='store_true',
        help='calculate sequence lengths after trimming primers and barcodes')
    parser.add_option('-Q', '--min-qual-score', type=int, default=25,
        help='min qual score allowed in read [default: %default]')
    parser.add_option('-k', '--keep-primer', action='store_true',
        help='do not remove primer from sequences')
    parser.add_option('-a', '--max-ambig', type=int, default=0,
        help='maximum number of ambiguous bases [default: %default]')
    parser.add_option('-H', '--max-homopolymer', type=int, default=6,
        help='maximum length of homopolymer run [default: %default]')
    parser.add_option('-M', '--max-primer-mismatch', dest='max_primer_mm',
        type=int, default=0,
        help='maximum number of primer mismatches [default: %default]')
    parser.add_option('-b', '--barcode-type', default='golay_12', 
        help=\
        'barcode type, e.g. 4 or hamming_8 or golay_12 [default: %default]')
    parser.add_option('-d', '--dir-prefix', default='.',
        help='directory prefix for output files [default: %default]')
    parser.add_option('-e', '--max-barcode-errors', dest='max_bc_errors',
        default=1.5, type=float,
        help='maximum number of errors in barcode [default: %default]')
    parser.add_option('-s', '--start-numbering-at', dest='start_index',
        default=1, type=int,
        help='seq id to use for the first sequence [default: %default]')
    options, args = parser.parse_args()

    required_options = [
        ('map_fname', '-m'),
        ('fasta_fnames', '-f'), 
        ('qual_fnames', '-q')]
    for attr_name, flag in required_options:
        if not getattr(options, attr_name):
            parser.error('Required option %s not found.  Must provide a '
                         'mapping file (-m), at least one fasta input file '
                         '(-i), and at least one qual file (-q).' % flag)
    return options


if __name__ == "__main__":
    from sys import exit, stderr
    options = make_cmd_parser()
    mapping_file = options.map_fname
    fasta_files = set(options.fasta_fnames.split(','))
    if options.qual_fnames:
        qual_files = set(options.qual_fnames.split(','))
    else:
        qual_files = set()

    for q in qual_files:
        if not q.endswith('qual'):
            stderr.write(
            "Qual file does not end with .qual: is it really a qual file?\n%s\n" 
            % q)

    for f in fasta_files:
        if not (f.endswith('fasta') or f.endswith('fna')):
            stderr.write(
            "Fasta file does not end with .fna: is it really a seq file?\n%s\n" 
            % f)

    preprocess(fasta_files, qual_files, mapping_file,
    primer_seq_pats = options.primers,
    barcode_type=options.barcode_type,
    starting_ix = options.start_index,
    min_seq_len = options.min_seq_len,
    max_seq_len = options.max_seq_len, 
    min_qual_score=options.min_qual_score,
    keep_primer=options.keep_primer,
    max_ambig=options.max_ambig,
    max_primer_mm=options.max_primer_mm,
    trim_seq_len=options.trim_seq_len,
    dir_prefix=options.dir_prefix,
    max_bc_errors = options.max_bc_errors,
    max_homopolymer = options.max_homopolymer,
    )
 
