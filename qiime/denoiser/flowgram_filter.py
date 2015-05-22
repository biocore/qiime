#!/usr/bin/env python

"""Functions to filter and crop flowgrams"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from re import compile, search
from itertools import imap
from collections import defaultdict
from tempfile import mkstemp
from os import close

from qiime.util import FileFormatError
from skbio.parse.sequences import parse_fasta
from bfillings.denoiser import Flowgram

from qiime.denoiser.utils import cat_sff_files, write_sff_header

# The default key sequence of the 454 machine
DEFAULT_KEYSEQ = "TCAG"


def filter_sff_file(flowgrams, header, filter_list, out_fh):
    """Filters all flowgrams in handle with filter.

    flowgrams: a list of flowgrams (or something similar)

    header: the header for the flowgrams

    filter_list: list of filters to be applied on sff.txt file

    out_fh: output file handle

    returns: number of flowgrams in filtered out file
    """

    write_sff_header(header, out_fh)

    l = 0
    for f in flowgrams:
        passed = True
        for filter in filter_list:
            passed = passed and filter(f)
            if not passed:
                # bail out
                break
        if (passed):
            out_fh.write(f.createFlowHeader() + "\n")
            l += 1
    return l


def within_length(flowgram, minlength=0, maxlength=400):
    """Checks if the (quality trimmed) seq of flowgram is within a specified length.

    flowgram: flowgram to check

    minlenght: minimal required length

    maxlenght: maximal allowed length
    """
    seq = flowgram.toSeq()
    l = len(seq)
    return (l >= minlength and l <= maxlength)


def check_ambigous(flowgram, max_allowed=4):
    """Check if there is a stretch of no-signals.

    flowgram: the flowgram to check

    max_allowed: number of consecutive signals below 0.5. This should be set to three,
                 or higher. 3 would mean no ambigous base called, 4 at most one N.
    """

    if max_allowed < 3:
        raise ValueError(
            "Error in calling check_ambigous. max_allowed should be at least 3")
    count = 0
    for signal in flowgram.flowgram:
        if signal < 0.5:
            count += 1
            if count > max_allowed:
                return True
        else:
            count = 0
    return False


def truncate_flowgrams_in_SFF(
        flowgrams, header, outhandle=None, outdir="/tmp/",
        barcode_mapping=None, primer=None,
        allow_num_ambigous=4):
    """Truncate flowgrams at low quality 3' end and strip key+primers.

    flowgrams: a list of flowgrams (or something similar)

    header: the header for the flowgrams

    outhandle: output file handle, can be None

    outdir: directory where random file will be created if outhandle is None

    barcode_mapping: dictionary mapping of read ids to barcode seqs.
                     The barcode seq will be truncated of the 5' end of the read

    primer: primer sequence that will be truncated of the 5' end of the read

    allow_num_ambigous: number of 'N' allowed in flowgram
    """
    out_filename = ""
    if not outhandle:
        fd, out_filename = mkstemp(dir=outdir, prefix="trunc_sff",
                                  suffix=".sff.txt")
        close(fd)
        outhandle = open(out_filename, "w")

    write_sff_header(header, outhandle)

    l = 0
    for f in flowgrams:
        qual_trimmed_flowgram = f.getQualityTrimmedFlowgram()

        if barcode_mapping:
            if f.Name in barcode_mapping:
                trunc_flowgram = qual_trimmed_flowgram.getPrimerTrimmedFlowgram(
                    primerseq=DEFAULT_KEYSEQ + barcode_mapping[f.Name] + primer)
            else:
                continue
        else:
            prim = DEFAULT_KEYSEQ
            if primer:
                prim += primer
            trunc_flowgram = qual_trimmed_flowgram.getPrimerTrimmedFlowgram(
                primerseq=prim)

        if(trunc_flowgram is not None):
            outhandle.write(trunc_flowgram.createFlowHeader() + "\n")
            l += 1
    return (out_filename, l)


def cleanup_sff(flowgrams, header, outhandle=None, outdir="/tmp",
                min_length=150, max_length=400):
    """Clean a sff file and returns name of clean file and number of clean flowgrams.

    flowgrams: a list of flowgrams

    header: the header for the flowgrams

    outhandle: handle flowgrams will be written to if set, can be stdout

    outdir: if outhandle is not set, random file will be created in outdir

    min_length, max_length: hard sequence contraints, default are set for
                 GS FLX, increase for Titanium, decrease for GS20

    NOTE: It is strongly recommended to use a proper quality filtering as in
    QIIME's split_libraries.py . This function is intended as a last resort and
    should  hardly ever be used,
    """

    clean_filename = ""
    if not outhandle:
        fd, clean_filename = mkstemp(dir=outdir, prefix="cleanup_sff",
                                    suffix=".sff.txt")
        close(fd)
        outhandle = open(clean_filename, "w")

    l = filter_sff_file(
        flowgrams, header, [lambda f: within_length(f, min_length, max_length),
                            lambda f: f.hasProperKey()],
        outhandle)
    return (clean_filename, l)


def split_sff(sff_file_handles, map_file_handle, outdir="/tmp/"):
    """Splits a sff.txt file on barcode/mapping file."""

    try:
        (flowgrams, header) = cat_sff_files(sff_file_handles)
    except ValueError:
        # reading in the binary sff usually shows up as ValueError
        raise FileFormatError('Wrong flogram file format. Make sure you pass the sff.txt format ' +
                              'produced by sffinfo. The binary .sff will not work here.')

    (inverse_map, map_count) = build_inverse_barcode_map(
        parse_fasta(map_file_handle))

    filenames = []
    # we might have many barcodes and reach python open file limit
    # therefor we go the slow way and open and close files each time
    # First set up all files with the headers only
    for barcode_id in map_count.keys():
        fh = open(outdir + barcode_id, "w")
        write_sff_header(header, fh, map_count[barcode_id])
        fh.close()
        filenames.append(outdir + barcode_id)
    # Then direct each flowgram into its barcode file
    for f in flowgrams:
        if f.Name in inverse_map:
            barcode_id = inverse_map[f.Name]
            fh = open(outdir + barcode_id, "a")
            fh.write(f.createFlowHeader() + "\n")
    return filenames


def build_inverse_barcode_map(seqs):
    """Build a map from fasta header from split_libraries.

    seqs: a list of (label, seq)  pairs

    Returns: mapping of flowgram ID to sampleID and a sample count

    A fasta header from split_libraries looks like this
    >S160_1 E86FECS01DW5V4 orig_bc=CAGTACGATCTT new_bc=CAGTACGATCTT bc_diffs=0
    """
    inverse_map = {}
    map_count = defaultdict(int)
    for (label, seq) in seqs:
        (map_id, seq_id) = label.split()[:2]
        map_id = map_id.split("_")[0]
        inverse_map[seq_id] = map_id
        map_count[map_id] += 1

    return (inverse_map, map_count)


def extract_barcodes_from_mapping(labels):
    """extract barcodes from split_libraries fasta headers.

    Returns a dictionary{flowgram_id:barcode_seq}
    """
    barcodes = {}

    # use \w* to allow for non barcoded reads
    re = compile(
        r'(\w+) ([a-zA-Z0-9.]+) orig_bc=(\w*) new_bc=\w* bc_diffs=\d+')
    for label in labels:
        tmatch = search(re, label)
        flowgram_id = tmatch.group(2)
        barcode = tmatch.group(3)

        barcodes[flowgram_id] = barcode

    return barcodes
