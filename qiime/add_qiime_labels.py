#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Emily TerAvest"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from os.path import join, basename
from string import letters, digits

from skbio.parse.sequences import parse_fasta

from qiime.util import duplicates_indices
from qiime.check_id_map import process_id_map


def add_qiime_labels(mapping_f,
                     fasta_dir,
                     filename_column,
                     output_dir=".",
                     count_start=0):
    """ Main function for combining fasta files, writing valid QIIME labels

    mapping_f:  open file object of the metadata mapping file
    fasta_dir:  Directory of fasta files to combine into a single file
    filename_column:  Column of metadata mapping file containing fasta filenames
    output_dir:  Directory to write output combined file to
    count_start:  Number to start enumeration of fasta labels with

    """

    headers, mapping_data, run_description, errors, warnings = \
        process_id_map(mapping_f, has_barcodes=False,
                       disable_primer_check=True, added_demultiplex_field=None,
                       variable_len_barcodes=False)

    fasta_name_to_sample_id = check_mapping_data(mapping_data, headers,
                                                 filename_column)

    fasta_files = get_fasta_fps(fasta_dir, fasta_name_to_sample_id.keys())

    write_combined_fasta(fasta_name_to_sample_id, fasta_files, output_dir,
                         counter=count_start)


def check_mapping_data(mapping_data, headers, filename_column):
    """ Checks mapping data for MIMARKS SampleIDs, unique IDs, fasta file names

    Also returns a dict of fasta file name: SampleID

    mapping_data:  list of lines of data from mapping file
    headers: list of header strings
    filename_column:  Column of metadata mapping file containing fasta filenames
    """

    # First make sure there is a SampleID and filename_column present
    try:
        sample_id_ix = headers.index("SampleID")
    except ValueError:
        raise ValueError("SampleID column not found in mapping file, please " +
                         "check mapping file with validate_mapping_file.py")

    try:
        filename_col_ix = headers.index(filename_column)
    except ValueError:
        raise ValueError("Specified column %s not found in mapping file." %
                         filename_column)

    valid_mimarks = letters + digits + "."

    fasta_name_to_sample_id = {}

    fasta_names = []
    sample_ids = []
    for line in mapping_data:

        try:
            fasta_name_to_sample_id[basename(line[filename_col_ix].strip())] =\
                line[sample_id_ix]
        except IndexError:
            raise IndexError("Missing filename column data in line %s " %
                             line)

        for curr_char in line[sample_id_ix]:
            if curr_char not in valid_mimarks:
                raise ValueError("Found invalid character in line: %s\n" %
                                 line + "SampleIDs must be alphanumeric and . characters " +
                                 "only")
        sample_ids.append(line[sample_id_ix].strip())
        fasta_names.append(line[filename_col_ix].strip())

    fasta_name_dups = duplicates_indices(fasta_names)
    if fasta_name_dups:
        raise ValueError("Found duplicate fasta names: %s" %
                         "\t".join([fasta_name for fasta_name in fasta_name_dups.keys()]))

    sample_id_dups = duplicates_indices(sample_ids)
    if sample_id_dups:
        raise ValueError("Found duplicate SampleID names: %s" %
                         "\t".join([sample_id for sample_id in sample_id_dups.keys()]))

    return fasta_name_to_sample_id


def get_fasta_fps(fasta_dir, fasta_files):
    """ Returns list of fasta filepaths (only .fna, .fasta, and .fa files)

    fasta_dir:  Directory of fasta files to check
    fasta_files:  list of fasta filenames to open
    """

    fasta_filepaths = []
    for curr_file in fasta_files:
        curr_fp = join(fasta_dir, curr_file)
        try:
            file_test = open(curr_fp, "U")
            file_test.close()
        except IOError:
            raise IOError("Unable to open %s" % curr_fp)
        fasta_filepaths.append(curr_fp)

    return fasta_filepaths


def write_combined_fasta(fasta_name_to_sample_id,
                         fasta_files,
                         output_dir=".",
                         counter=0):
    """ Writes combined, enumerated fasta file

    fasta_name_to_sample_id:  dict of fasta file name to SampleID
    fasta_files: list of filepaths to iterate through
    output_dir:  output directory to write combined file to
    counter:  Starting number to enumerate sequences with
    """

    combined_file_out = open(join(output_dir + "/", "combined_seqs.fna"), "w")

    for curr_fasta in fasta_files:
        for label, seq in parse_fasta(open(curr_fasta, "U")):
            combined_file_out.write(">%s_%d %s\n" %
                                    (fasta_name_to_sample_id[basename(curr_fasta)], counter, label))
            combined_file_out.write("%s\n" % seq)
            counter += 1
