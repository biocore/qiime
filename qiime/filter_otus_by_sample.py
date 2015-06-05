#!/usr/bin/env python
# file filter_otus_by_sample.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stombaugh"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"


from string import strip
import re

from skbio.alignment import SequenceCollection
from skbio.sequence import DNA

def filter_otus(otus, prefs):
    """filters the otus file based on which samples should be removed and
       determines which sequences to remove"""
    new_otus_list = []
    # create an list containing the seqs to remove

    for i in otus:
        new_otus = []
        for j in (otus[i]):
            sample_seq = j.split('_')
            if len(sample_seq) > 1:
                new_name = ''.join(sample_seq[:-1])
            else:
                new_name = sample_seq[0]
            is_sample = False
            for sample_id in prefs:
                if prefs[sample_id] == new_name:
                    is_sample = True
            if is_sample:
                pass
            else:
                new_otus.append('%s' % (j))
        if new_otus:
            new_otus_list.append((i, new_otus))
    return new_otus_list


def filter_aln_by_otus(aln, prefs):
    """filters the representative set of seqs based on which samples should
        be removed"""
    filtered_seqs = []
    removed_seqs = []
    for j in range(aln.sequence_count()):
        remove = False
        aln_name = aln[j].id
        stripped_aln_name = aln_name.split(' ')[0].split('_')
        if len(stripped_aln_name) > 1:
            new_aln_name = ''.join(stripped_aln_name[:-1])
        else:
            new_aln_name = stripped_aln_name[0]

        for sample_id in prefs:
            if prefs[sample_id] == new_aln_name:
                remove = True

        if remove:
            removed_seqs.append((aln_name, str(aln[aln_name])))
        else:
            filtered_seqs.append((aln_name, str(aln[aln_name])))

    return filtered_seqs, removed_seqs


def process_extract_samples(samples_to_extract):
    """Parses the samples to extract option from the command line"""
    prefs = {}

    if samples_to_extract:
        samples = samples_to_extract.strip().strip("'").split(',')

    for j, col in enumerate(samples):
        key = str(j)
        prefs[key] = {}
        prefs[key] = col

    return prefs


def filter_samples(prefs, data, dir_path='', filename=None):
    """processes the filtering of the otus file and representative seq set, then
        writes filtered otus and filtered representative seq set files"""

    aln = data['aln']
    otus = data['otus']

    # filter the otus file based on which samples to remove
    new_otus_list = filter_otus(otus, prefs)

    filtered_otus_output_filepath = '%s/%s_sfiltered_otus.txt' \
                                    % (dir_path, filename)
    filtered_otus_output_filepath = open(filtered_otus_output_filepath, 'w')

    # Write out a new otus file
    for key in (new_otus_list):
        filtered_otus_output_filepath.write(key[0])
        for j in key[1]:
            filtered_otus_output_filepath.write('\t' + str(j))
        filtered_otus_output_filepath.write('\n')
    filtered_otus_output_filepath.close()

    # filter seq set
    filtered_seqs, removed_seqs = filter_aln_by_otus(aln, prefs)

    # write a fasta containing list of sequences removed from
    # representative set
    if len(removed_seqs) > 0:
        removed_seqs = SequenceCollection.from_fasta_records(
                [(e[0], str(e[1])) for e in removed_seqs], DNA)
    else:
        raise ValueError(
            'No sequences were removed.  Did you specify the correct Sample ID?')
    output_filepath2 = '%s/%s_sremoved.fasta' % (dir_path, filename)
    output_file2 = open(output_filepath2, 'w')
    output_file2.write(removed_seqs.to_fasta())
    output_file2.close()

    # write a fasta containing the filtered representative seqs
    if len(filtered_seqs) > 0:
        filtered_seqs = SequenceCollection.from_fasta_records(
                [(e[0], str(e[1])) for e in filtered_seqs], DNA)
    else:
        raise ValueError(
            'No sequences were remaining in the fasta file.  Did you remove all Sample ID\'s?')

    output_filepath = '%s/%s_sfiltered.fasta' % (dir_path, filename)
    output_file = open(output_filepath, 'w')
    output_file.write(filtered_seqs.to_fasta())
    output_file.close()
