#!/usr/bin/env python
# File created Feb 1 2012

from __future__ import division

__author__ = "William Anton Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Anton Walters", "Emily TerAvest"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Anton Walters"
__email__ = "william.a.walters@gmail.com"

from collections import defaultdict
from os.path import split, join

from skbio.parse.sequences import parse_fasta
from skbio.io import RecordError
from cogent.parse.tree import DndParser
from qiime.check_id_map import process_id_map
from qiime.split_libraries import expand_degeneracies


def get_mapping_details(mapping_fp,
                        suppress_barcode_checks=False,
                        suppress_primer_checks=False):
    """ Returns SampleIDs, Barcodes, Primer seqs from mapping file

    mapping_fp: filepath to mapping file
    suppress_barcode_checks=If True, will skip getting barcodes from mapping
     file and searching for these in sequences.
    suppress_primer_checks=If True, will skip getting primers from mapping
     file and searching for these in sequences
    """

    mapping_f = open(mapping_fp, "U")

    # Only using the id_map and the errors from parsing the mapping file.
    hds, mapping_data, run_description, errors, warnings = \
        process_id_map(mapping_f)

    mapping_f.close()

    # Should raise errors for barcodes or primers unless suppressed, and
    # should raise errors for headers or duplicate SampleIDs in any case.
    loc_bcs = ",1"
    loc_primers = ",2"
    if errors:
        for curr_error in errors:
            # Halt when header has error
            if curr_error.startswith("Found header field"):
                raise ValueError('Error in mapping file, please validate '
                                 'mapping file with validate_mapping_file.py')
            elif curr_error.endswith(loc_bcs):
                # Halt for barcode errors unless suppressed
                if suppress_barcode_checks:
                    continue
                else:
                    raise ValueError('Error in mapping file, please validate '
                                     'mapping file with validate_mapping_file.py')
            elif curr_error.endswith(loc_primers):
                # Halt for primer errors unless suppressed
                if suppress_primer_checks:
                    continue
                else:
                    raise ValueError('Error in mapping file, please validate '
                                     'mapping file with validate_mapping_file.py')
            # Raise error on duplicate sample IDs
            elif curr_error.startswith("Duplicate SampleID"):
                raise ValueError('Error in mapping file, please validate '
                                 'mapping file with validate_mapping_file.py')

    # create dict of dicts with SampleID:{each header:mapping data}

    id_map = {}

    for curr_data in mapping_data:
        id_map[curr_data[0]] = {}

    for header in range(len(hds)):
        for curr_data in mapping_data:
            id_map[curr_data[0]][hds[header]] = curr_data[header]

    sample_ids = id_map.keys()

    barcode_seqs = []
    raw_linkerprimer_seqs = []

    for curr_id in id_map:
        if not suppress_barcode_checks:
            barcode_seqs.append(id_map[curr_id]['BarcodeSequence'])
        if not suppress_primer_checks:
            raw_linkerprimer_seqs.append(
                id_map[curr_id]['LinkerPrimerSequence'])

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
        for label, seq in parse_fasta(fasta_f):
            continue
    except RecordError:
        raise RecordError("Input fasta file not valid fasta format.  Error " +
                          "found at %s label and %s sequence " % (label, seq))

    fasta_f.close()


def get_fasta_labels(input_fasta_fp):
    """ Returns the fasta labels (text before whitespace) as a list

    input_fasta_fp: fasta filepath
    """

    fasta_labels = []

    fasta_f = open(input_fasta_fp, "U")

    for label, seq in parse_fasta(fasta_f):
        fasta_labels.append(label.split()[0])

    return fasta_labels


def get_dup_labels_perc(fasta_labels):
    """ Calculates percentage of sequences with duplicate labels

    fasta_labels: list of fasta labels
    """
    fasta_labels_count = float(len(fasta_labels))
    fasta_labels_derep = float(len(set(fasta_labels)))

    perc_dup = "%1.3f" %\
        ((fasta_labels_count - fasta_labels_derep) / fasta_labels_count)

    label_counts = defaultdict(int)
    for curr_label in fasta_labels:
        label_counts[curr_label] += 1

    labels_from_dups = []
    for label in label_counts:
        if label_counts[label] > 1:
            labels_from_dups.append(label)

    return perc_dup, labels_from_dups


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

    total_seq_count = float(total_seq_count)
    valid_id_count = float(valid_id_count)
    matches_sampleid_count = float(matches_sampleid_count)

    perc_not_valid = "%1.3f" %\
        ((total_seq_count - valid_id_count) / total_seq_count)
    perc_nosampleid_match = "%1.3f" %\
        ((total_seq_count - matches_sampleid_count) / total_seq_count)

    return perc_not_valid, perc_nosampleid_match


def check_fasta_seqs(input_fasta_fp,
                     barcodes,
                     linkerprimerseqs,
                     total_seq_count,
                     valid_chars=frozenset(['A', 'T', 'C', 'G', 'N', 'a',
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
    barcodes_at_start = 0

    # Get max barcode length to checking the beginning of seq for barcode
    if barcodes:
        max_bc_len = max([len(bc_len) for bc_len in barcodes])
    else:
        max_bc_len = 0

    for label, seq in parse_fasta(input_fasta_f):

        # Only count one offending problem
        for curr_nt in seq:
            if curr_nt not in valid_chars:
                invalid_chars_count += 1
                break

        sliced_seq = seq[0:max_bc_len]

        for curr_bc in barcodes:
            if curr_bc in sliced_seq:
                barcodes_at_start += 1
                break

        for curr_bc in barcodes:
            if curr_bc in seq:
                barcodes_count += 1
                break

        for curr_primer in linkerprimerseqs:
            if curr_primer in seq:
                linkerprimers_count += 1
                break

    invalid_chars_count = float(invalid_chars_count)
    barcodes_count = float(barcodes_count)
    linkerprimers_count = float(linkerprimers_count)
    total_seq_count = float(total_seq_count)
    barcodes_at_start_count = float(barcodes_at_start)

    perc_invalid_chars = "%1.3f" %\
        (invalid_chars_count / total_seq_count)
    perc_barcodes_detected = "%1.3f" %\
        (barcodes_count / total_seq_count)
    perc_primers_detected = "%1.3f" %\
        (linkerprimers_count / total_seq_count)
    perc_barcodes_at_start_detected = "%1.3f" %\
        (barcodes_at_start_count / total_seq_count)

    return perc_invalid_chars, perc_barcodes_detected, perc_primers_detected,\
        perc_barcodes_at_start_detected


def check_fasta_seqs_lens(input_fasta_fp):
    """ Creates bins of sequence lens

    Useful for checking for valid aligned sequences.

    input_fasta_fp:  input fasta filepath
    """

    seq_lens = defaultdict(int)

    input_fasta_f = open(input_fasta_fp, "U")

    for label, seq in parse_fasta(input_fasta_f):
        seq_lens[len(seq)] += 1

    input_fasta_f.close()

    formatted_seq_lens = []

    for curr_key in seq_lens:
        formatted_seq_lens.append((seq_lens[curr_key], curr_key))

    formatted_seq_lens.sort(reverse=True)

    return formatted_seq_lens


def check_all_ids(fasta_labels,
                  sample_ids):
    """ Tests that all sample IDs from mapping file are found in seq labels

    fasta_labels: list of fasta labels
    sample_ids: set of sample ids from mapping file
    """

    # Need to get modified fasta labels with underscore stripped

    raw_fasta_labels = set([label.split('_')[0] for label in fasta_labels])

    sample_ids_not_found = []

    for curr_id in sample_ids:
        if curr_id not in raw_fasta_labels:
            sample_ids_not_found.append(curr_id)

    # Return True if all were found, otherwise list of sampleIDs
    if len(sample_ids_not_found) == 0:
        sample_ids_not_found = True

    return sample_ids_not_found


def check_tree_subset(fasta_labels,
                      tree_fp):
    """ Returns a list of all fasta labels that are not a subset of the tree

    fasta_labels:  list of fasta labels
    tree_fp: tree filepath
    """

    # Need to get modified fasta labels with underscore stripped

    raw_fasta_labels = set([label.split('_')[0] for label in fasta_labels])

    tree_f = open(tree_fp, "U")

    tree = DndParser(tree_f)

    # Get a set of tree tip names
    tree_tips = set(tree.getTipNames())

    labels_not_in_tips = []

    for curr_label in raw_fasta_labels:
        if curr_label not in tree_tips:
            labels_not_in_tips.append(curr_label)

    # Return True if all found in tree tips
    if len(labels_not_in_tips) == 0:
        labels_not_in_tips = True

    return labels_not_in_tips


def check_tree_exact_match(fasta_labels,
                           tree_fp):
    """Checks fasta labels to exact match to tree tips

    Returns a list of two lists, the fasta labels not in tips, and tips not
     in fasta labels.
    fasta_labels: list of fasta labels
    tree_fp: tree filepath
    """

    # Need to get modified fasta labels with underscore stripped

    raw_fasta_labels = set([label.split('_')[0] for label in fasta_labels])

    tree_f = open(tree_fp, "U")

    tree = DndParser(tree_f)

    # Get a set of tree tip names
    tree_tips = set(tree.getTipNames())

    labels_not_in_tips = []

    for curr_label in raw_fasta_labels:
        if curr_label not in tree_tips:
            labels_not_in_tips.append(curr_label)

    # Return True if all found in tree tips
    if len(labels_not_in_tips) == 0:
        labels_not_in_tips = True

    tips_not_in_labels = []

    for curr_tip in tree_tips:
        if curr_tip not in raw_fasta_labels:
            tips_not_in_labels.append(curr_tip)

    if len(tips_not_in_labels) == 0:
        tips_not_in_labels = True

    return [labels_not_in_tips, tips_not_in_labels]


def run_fasta_checks(input_fasta_fp,
                     mapping_fp,
                     tree_fp=None,
                     tree_subset=False,
                     tree_exact_match=False,
                     same_seq_lens=False,
                     all_ids_found=False,
                     suppress_barcode_checks=False,
                     suppress_primer_checks=False):
    """ Returns dictionary of records for different fasta checks

    input_fasta_fp: fasta filepath
    mapping_fp: mapping filepath
    tree_fp: newick tree filepath
    tree_subset: If True, will test that SampleIDs are a subset of the tree tips
    tree_exact_match: If True, will test that SampleIDs are an exact match to
     the tree
    same_seq_lens: If True, will determine if sequences are of different lens.
    all_ids_found: If True, will determine if all SampleIDs are represented
     in the sequence labels.
    suppress_barcode_checks=If True, will skip getting barcodes from mapping
     file and searching for these in sequences.
    suppress_primer_checks=If True, will skip getting primers from mapping
     file and searching for these in sequences"""

    # Stores details of various checks
    fasta_report = {}

    # get sets of data for testing fasta labels/seqs
    sample_ids, barcodes, linkerprimerseqs = get_mapping_details(mapping_fp,
                                                                 suppress_barcode_checks, suppress_primer_checks)

    fasta_labels = get_fasta_labels(input_fasta_fp)

    total_seq_count = len(fasta_labels)

    fasta_report['duplicate_labels'], fasta_report['duplicate_ids'] =\
        get_dup_labels_perc(fasta_labels)

    fasta_report['invalid_labels'], fasta_report['nosample_ids_map'] =\
        check_labels_sampleids(fasta_labels, sample_ids, total_seq_count)

    fasta_report['invalid_seq_chars'], fasta_report['barcodes_detected'],\
        fasta_report['linkerprimers_detected'],\
        fasta_report['barcodes_at_start'] = check_fasta_seqs(input_fasta_fp,
                                                             barcodes, linkerprimerseqs, total_seq_count)

    if same_seq_lens:
        fasta_report['same_seq_lens'] = check_fasta_seqs_lens(input_fasta_fp)
    else:
        fasta_report['same_seq_lens'] = False

    if all_ids_found:
        fasta_report['all_ids_found'] = check_all_ids(fasta_labels, sample_ids)
    else:
        fasta_report['all_ids_found'] = False

    if tree_subset:
        fasta_report['tree_subset'] = check_tree_subset(fasta_labels, tree_fp)
    else:
        fasta_report['tree_subset'] = False

    if tree_exact_match:
        fasta_report['tree_exact_match'] =\
            check_tree_exact_match(fasta_labels, tree_fp)
    else:
        fasta_report['tree_exact_match'] = False

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

    output_f.write("Percent duplicate labels: %s\n" %
                   fasta_report['duplicate_labels'])
    output_f.write("Percent QIIME-incompatible fasta labels: %s\n" %
                   fasta_report['invalid_labels'])
    output_f.write("Percent of labels that fail to map to SampleIDs: %s\n" %
                   fasta_report['nosample_ids_map'])
    output_f.write("Percent of sequences with invalid characters: %s\n" %
                   fasta_report['invalid_seq_chars'])
    output_f.write("Percent of sequences with barcodes detected: %s\n" %
                   fasta_report['barcodes_detected'])
    output_f.write("Percent of sequences with barcodes detected at the " +
                   "beginning of the sequence: %s\n" % fasta_report['barcodes_at_start'])
    output_f.write("Percent of sequences with primers detected: %s\n" %
                   fasta_report['linkerprimers_detected'])

    # Optional tests

    if fasta_report['same_seq_lens']:
        output_f.write("Sequence lengths report\n")
        output_f.write("Counts of sequences, followed by their sequence " +
                       "lengths:\n")
        for len_data in fasta_report['same_seq_lens']:
            output_f.write("%s\t%s\n" % (len_data[0], len_data[1]))

    # need to make an explicit check for true as there can be other non-boolean

    # values stored under this key in the dictionary; this needs to be fixed

    if fasta_report['all_ids_found']:
        output_f.write("Sample ID in fasta sequences report\n")
        if fasta_report['all_ids_found'] == True:

            output_f.write("All SampleIDs found in sequence labels.\n")
        else:
            output_f.write("The following SampleIDs were not found:\n")
            for curr_id in fasta_report['all_ids_found']:
                output_f.write("%s\n" % curr_id)

    if fasta_report['tree_subset']:
        output_f.write("Fasta label subset in tree tips report\n")
        if fasta_report['tree_subset'] == True:

            output_f.write("All fasta labels were a subset of tree tips.\n")
        else:
            output_f.write("The following labels were not in tree tips:\n")
            for curr_id in fasta_report['tree_subset']:
                output_f.write("%s\n" % curr_id)

    if fasta_report['tree_exact_match']:
        output_f.write("Fasta label/tree tip exact match report\n")
        if fasta_report['tree_exact_match'][0] == True:

            output_f.write("All fasta labels found in tree tips.\n")
        else:
            output_f.write("The following labels were not in tree tips:\n")
            for curr_label in fasta_report['tree_exact_match'][0]:
                output_f.write("%s\n" % curr_label)
        if fasta_report['tree_exact_match'][1] == True:

            output_f.write("All tree tips found in fasta labels.\n")
        else:
            output_f.write("The following tips were not in fasta labels:\n")
            for curr_tip in fasta_report['tree_exact_match'][1]:
                output_f.write("%s\n" % curr_tip)

    if fasta_report['duplicate_ids']:
        output_f.write("Duplicate labels found:\n")
        for id in fasta_report['duplicate_ids']:
            output_f.write("%s\n" % id)


def validate_fasta(input_fasta_fp,
                   mapping_fp,
                   output_dir,
                   tree_fp=None,
                   tree_subset=False,
                   tree_exact_match=False,
                   same_seq_lens=False,
                   all_ids_found=False,
                   suppress_barcode_checks=False,
                   suppress_primer_checks=False):
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
    suppress_barcode_checks=If True, will skip getting barcodes from mapping
     file and searching for these in sequences.
    suppress_primer_checks=If True, will skip getting primers from mapping
     file and searching for these in sequences
    """

    # First test is valid fasta format, can't do other tests if file can't be
    # parsed so will bail out here if invalid

    verify_valid_fasta_format(input_fasta_fp)

    fasta_report = run_fasta_checks(input_fasta_fp, mapping_fp, tree_fp,
                                    tree_subset, tree_exact_match, same_seq_lens, all_ids_found,
                                    suppress_barcode_checks, suppress_primer_checks)

    write_log_file(output_dir, input_fasta_fp, fasta_report)
