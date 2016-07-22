#!/usr/bin/env python
# file trim_sff_primers.py: resets trim values in sff file based on primers.

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", 'Kyle Bittinger', "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"

"""Finds the technical read regions for each library, and resets the left trim.

Replaces the sff files with the trimmed versions.
"""

from itertools import imap
from os import walk, devnull, remove, rename, close
from os.path import splitext, join, exists
from shutil import move
from subprocess import check_call
from sys import stderr
from tempfile import mkstemp, TemporaryFile

from burrito.util import ApplicationNotFoundError
from cogent.parse.binary_sff import parse_binary_sff, write_binary_sff
from burrito.util import which

from qiime.parse import parse_mapping_file


def get_technical_lengths(input_map, debug=False):
    """Returns per-sample info on technical lengths.

    Note: KEY_SEQ, BARCODE and PRIMER fields are required. The LINKER
    field is optional.
    """
    if debug:
        print "Making debug output"
    body, header, comments = parse_mapping_file(input_map)
    if debug:
        print "HEADER:", header
    key_index = header.index('KEY_SEQ')
    bc_index = header.index('BARCODE')
    if 'LINKER' in header:
        linker_index = header.index('LINKER')
    else:
        linker_index = None
    primer_index = header.index('PRIMER')
    technical_lengths = {}
    for fields in body:
        curr_tech_len = len(fields[key_index]) + len(fields[bc_index]) + \
            len(fields[primer_index])
        if linker_index is not None:
            curr_tech_len += len(fields[linker_index])
        technical_lengths[fields[0]] = curr_tech_len
    if debug:
        print "Technical lengths:"
        print technical_lengths
    return technical_lengths


def get_per_lib_sff_fps(sff_dir):
    """Return a dict mapping library ID's to SFF filepaths in a directory.
    """
    for dirpath, dirnames, fnames in walk(sff_dir):
        for fname in fnames:
            if fname.endswith('.sff'):
                libname, _ = splitext(fname)
                yield libname, join(dirpath, fname)


def set_clip_qual_left(sff_data, technical_read_length):
    """Resets the value of clip_qual_left for each read in the SFF data.

    For efficiency, the original read objects are not copied, and are
    modified as a side effect of this function.
    """
    # TODO: Move to PyCogent
    header, reads = sff_data

    # sfftools use 1-based indexing
    clip_idx = technical_read_length + 1

    def adjust_read(read):
        read['clip_qual_left'] = clip_idx
        return read

    return header, imap(adjust_read, reads)


def set_sff_trimpoints(sff_dir, technical_lengths):
    """Set trimpoints to end of technical read for all SFF files in directory.
    """
    for lib_id, sff_fp in get_per_lib_sff_fps(sff_dir):
        try:
            readlength = technical_lengths[lib_id]
        except KeyError:
            continue
        sff_data = parse_binary_sff(open(sff_fp), True)
        clipped_header, clipped_reads = set_clip_qual_left(
            sff_data, readlength)

        fd, temp_fp = mkstemp(dir=sff_dir)
        close(fd)
        with open(temp_fp, 'w') as f:
            write_binary_sff(f, clipped_header, clipped_reads)

        move(temp_fp, sff_fp)


def set_sff_trimpoints_with_sfftools(
        sff_dir, technical_lengths, sffinfo_path='sffinfo', sfffile_path='sfffile',
        debug=False):
    """Set trimpoints to end of technical read for all SFF files in directory.

    This function essentially provides the reference implementation.
    It uses the official sfftools from Roche to process the SFF files.
    """
    if not (exists(sffinfo_path) or which(sffinfo_path)):
        raise ApplicationNotFoundError(
            'sffinfo executable not found. Is it installed and in your $PATH?')
    if not (exists(sfffile_path) or which(sfffile_path)):
        raise ApplicationNotFoundError(
            'sfffile executable not found. Is it installed and in your $PATH?')

    for lib_id, sff_fp in get_per_lib_sff_fps(sff_dir):
        try:
            readlength = technical_lengths[lib_id]
        except KeyError:
            continue

        sffinfo_args = [sffinfo_path, '-s', sff_fp]
        if debug:
            print "Running sffinfo command %s" % sffinfo_args
        sffinfo_output_file = TemporaryFile()
        check_call(sffinfo_args, stdout=sffinfo_output_file)
        sffinfo_output_file.seek(0)

        seqlengths = {}
        for line in sffinfo_output_file:
            if line.startswith('>'):
                fields = line[1:].split()
                seq_len = fields[1].split('=')[1]
                seqlengths[fields[0]] = seq_len

        trim_fp = sff_fp + '.trim'
        trim_file = open(trim_fp, 'w')
        for id_, length in seqlengths.items():
            curr_length = int(seqlengths[id_])
            # Sfftools use 1-based index
            left_trim = readlength + 1
            # Key sequence not included in FASTA length
            right_trim = curr_length + 4
            if curr_length > left_trim:
                trim_file.write(
                    "%s\t%s\t%s\n" % (id_, left_trim, right_trim))
            else:
                stderr.write(
                    'Rejected read %s with trim points %s and %s (orig '
                    'length %s)' % (id_, left_trim, curr_length, length))
        trim_file.close()

        trimmed_sff_fp = sff_fp + '.trimmed'
        sfffile_args = [
            sfffile_path, '-t', trim_fp, '-o', trimmed_sff_fp, sff_fp]
        if debug:
            print "Running sfffile command:", sfffile_args
        check_call(sfffile_args, stdout=open(devnull, 'w'))
        remove(sff_fp)
        rename(trimmed_sff_fp, sff_fp)
