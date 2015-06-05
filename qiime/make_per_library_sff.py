#!/usr/bin/env python

from __future__ import division

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"

import itertools
import os
import subprocess

from qiime.process_sff import (
    check_sfffile,
)
from cogent.parse.binary_sff import (
    parse_binary_sff, write_binary_sff,
)


def filter_sff_reads(sff_data, ids_to_keep=None, ids_to_remove=None):
    """Retain reads where the ID is in ids_to_keep but not in ids_to_remove.

    This function reproduces the behavior of the -i and -e options in
    Roche's sfffile program.
    """
    # TODO: Move to PyCogent
    header, reads = sff_data
    # Force evaluation of all reads.  We have no choice, since we need
    # the total number of reads to be returned with the header.
    # Another design choice would be to go back and correct the header
    # when we write the binary SFF file to disk -- maybe we'll switch
    # to that strategy in the future.
    if ids_to_keep is not None:
        reads = [r for r in reads if r['Name'] in ids_to_keep]
    if ids_to_remove is not None:
        reads = [r for r in reads if r['Name'] not in ids_to_remove]
    header['number_of_reads'] = len(reads)
    return header, reads


def parse_id_list(id_list_file):
    ids = set()
    for line in id_list_file:
        words = line.split()
        if words:
            ids.add(words[0].lstrip('>'))
    return ids


def combine_sff_data(*sff_datasets):
    combined_header = {'number_of_reads': 0}
    combined_reads = []

    for header, reads in sff_datasets:
        prev_num_reads = combined_header['number_of_reads']
        combined_header = header.copy()
        combined_header['number_of_reads'] += prev_num_reads

        combined_reads = itertools.chain(combined_reads, reads)

    combined_header['index_offset'] = 0
    combined_header['index_length'] = 0
    return combined_header, combined_reads


def make_per_library_sff(sff_fps, id_list_fp, debug=False):
    id_list_basepath, _ = os.path.splitext(id_list_fp)
    output_fp = id_list_basepath + '.sff'

    sff_datasets = [parse_binary_sff(open(fp), True) for fp in sff_fps]
    sff_data = combine_sff_data(*sff_datasets)
    ids = parse_id_list(open(id_list_fp))

    filtered_sff_data = filter_sff_reads(sff_data, ids_to_keep=ids)
    if debug:
        print 'Creating SFF file for %s' % id_list_fp
    write_binary_sff(open(output_fp, 'w'), *filtered_sff_data)


def make_per_library_sff_with_sfffile(
        sff_fps, id_list_fp, sfffile_path=None, debug=False):
    id_list_basepath, _ = os.path.splitext(id_list_fp)
    output_fp = id_list_basepath + '.sff'

    check_sfffile()
    args = ['sfffile', '-i', id_list_fp, '-o', output_fp] + sff_fps
    if debug:
        print args
    subprocess.check_call(
        args, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))


def make_per_library_sffs(
        sff_fps, id_list_dir, use_sfftools=False, sfffile_path=None, debug=False):
    for dirpath, dirnames, filenames in os.walk(id_list_dir):
        for filename in filenames:
            if filename.startswith('.'):
                continue
            id_list_fp = os.path.join(dirpath, filename)

            if use_sfftools:
                make_per_library_sff_with_sfffile(
                    sff_fps, id_list_fp, sfffile_path, debug)
            else:
                make_per_library_sff(sff_fps, id_list_fp, debug)
