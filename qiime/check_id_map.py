#!/usr/bin/env python
from __future__ import division

"""Parse mapping file, checking for a number of undesirable characteristics.

Specifically, we check that:
    - The BarcodeSequence, LinkerPrimerSequences, and ReversePrimer fields
       have valid IUPAC DNA characters, and BarcodeSequence characters
       are non-degenerate (error)
    - The SampleID, BarcodeSequence, LinkerPrimerSequence, and Description
       headers are present. (error)
    - There are not duplicate header fields (error)
    - There are not duplicate barcodes (error)
    - Barcodes are of the same length.  Suppressed when
       variable_len_barcode flag is passed (warning)
    - The headers do not contain invalid characters (alphanumeric and
       underscore only) (warning)
    - The data fields do not contain invalid characters (alphanumeric,
       underscore, space, and +-%./:,; characters) (warning)
    - SampleID fields are MIENS compliant (only alphanumeric
      and . characters) (warning)
    - There are no duplicates when the primer and variable length
       barcodes are appended (error)
    - There are no duplicates when barcodes and added demultiplex
       fields (-j option) are combined (error)
    - Data fields are not found beyond the Description column (warning)

Errors and warnings will be appended to an initially empty list and returned.

A log file will be created containing a list of all errors and warnings
detected.  Header errors can preclude the generation of errors or warnings
in data cells.

A _corrected_mapping.txt file will be created in the output directory with
these characters replaced with the char_replace parameter character.  If *any*
other warnings or errors are detected, the _corrected_mapping.txt file will
contain a message that the _corrected file to this effect, and will
reference the log and html file for the user to correct the problems manually.

"""

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "William Walters"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"

from collections import defaultdict
from string import letters, digits, upper
from os.path import basename, join
from operator import itemgetter
from copy import deepcopy
from shutil import copyfile

from skbio.sequence import DNASequence

from qiime.util import get_qiime_project_dir, duplicates_indices
from qiime.parse import parse_mapping_file
from qiime.format import format_mapping_html_data


def check_mapping_file(mapping_fp,
                       output_dir=".",
                       has_barcodes=True,
                       char_replace="_",
                       verbose=True,
                       variable_len_barcodes=False,
                       disable_primer_check=False,
                       added_demultiplex_field=None,
                       suppress_html=False):
    """ Main program function for checking mapping file

    Checks mapping file for errors, warnings, writes log file, html file,
    and corrected mapping file.

    mapping_fp:  path to metadata mapping file
    output_dir:  output directory for log, html, corrected mapping file.
    has_barcodes:  If True, will test for perform barcodes test (presence,
     uniqueness, valid IUPAC DNA chars).
    char_replace:  Character used to replace invalid characters in data
     fields.  SampleIDs always use periods to be MIENS compliant.
    verbose: If True, a message about warnings and/or errors will be printed
     to stdout.
    variable_len_barcodes:  If True, suppresses warnings about barcodes of
     varying length.
    disable_primer_check:  If True, disables tests for valid primer sequences.
    added_demultiplex_field:  If specified, references a field in the mapping
     file to use for demultiplexing.  These are to be read from fasta labels
     during the actual demultiplexing step.  All combinations of barcodes,
     primers, and the added_demultiplex_field must be unique."""

    header, mapping_data, run_description, errors, warnings =\
        process_id_map(open(mapping_fp, 'U'), disable_primer_check,
                       has_barcodes, char_replace, variable_len_barcodes,
                       added_demultiplex_field, strip_quotes=False, suppress_stripping=True)

    if not suppress_html:
        formatted_html = format_mapping_html_data(header, mapping_data,
                                                  errors, warnings)

        output_html = join(output_dir +
                           basename(mapping_fp).replace('.txt', '') + ".html")

        html_f = open(output_html, "w")
        html_f.write(formatted_html)

        # get QIIME directory
        qiime_dir = get_qiime_project_dir()

        # Write javascript file necessary for mouseover tooltips.
        # move javascript file to javascript output directory
        copyfile(join(qiime_dir, 'qiime', 'support_files',
                      'js/overlib.js'), join(output_dir, 'overlib.js'))

    corrected_mapping_data = correct_mapping_data(mapping_data,
                                                  header, char_replace)

    output_corrected_fp = join(output_dir +
                               basename(mapping_fp).replace('.txt', '') + "_corrected.txt")

    write_corrected_mapping(output_corrected_fp, header, run_description,
                            corrected_mapping_data)

    output_log_fp = join(output_dir +
                         basename(mapping_fp).replace('.txt', '') + ".log")

    write_log_file(output_log_fp, errors, warnings)

    if verbose:
        if errors or warnings:
            print "Errors and/or warnings detected in mapping file.  Please " +\
                "check the log and html file for details."
        else:
            print "No errors or warnings were found in mapping file."


def process_id_map(mapping_f,
                   disable_primer_check=False,
                   has_barcodes=True,
                   char_replace="_",
                   variable_len_barcodes=False,
                   added_demultiplex_field=None,
                   strip_quotes=True,
                   suppress_stripping=False):
    """ Reads mapping file, returns data, warnings, and errors

    mapping_f:  list of lines (open metadata mapping file object)
    has_barcodes:  If True, will test for perform barcodes test (presence,
     uniqueness, valid IUPAC DNA chars).
    char_replace:  Character used to replace invalid characters in data
     fields.  SampleIDs always use periods to be MIENS compliant.
    variable_len_barcodes:  If True, suppresses warnings about barcodes of
     varying length.
    disable_primer_check:  If True, disables tests for valid primer sequences.
    added_demultiplex_field:  If specified, references a field in the mapping
     file to use for demultiplexing.  These are to be read from fasta labels
     during the actual demultiplexing step.  All combinations of barcodes,
     primers, and the added_demultiplex_field must be unique.
    strip_quotes: Sets stripping of quote characters from the mapping file.
    suppress_stripping: suppresses stripping of white space from mapping file.

    """

    errors = []
    warnings = []

    mapping_data, header, comments = parse_mapping_file(
        mapping_f, strip_quotes,
        suppress_stripping)

    sample_id_ix = 0
    # Get index of last field of header
    desc_ix = len(header) - 1
    bc_ix = 1
    linker_primer_ix = 2

    # Find errors/warnings in header fields
    errors, warnings = check_header(header, errors,
                                    warnings, sample_id_ix, desc_ix, bc_ix, linker_primer_ix,
                                    added_demultiplex_field)

    # Find errors/warnings in data fields, get corrected form with invalid
    # characters replaced
    errors, warnings = check_data_fields(header, mapping_data,
                                         errors, warnings, disable_primer_check, has_barcodes, char_replace,
                                         variable_len_barcodes, added_demultiplex_field)

    return header, mapping_data, comments, errors, warnings

# Being data field checking functions


def check_data_fields(header,
                      mapping_data,
                      errors,
                      warnings,
                      disable_primer_check=False,
                      has_barcodes=True,
                      char_replace="_",
                      variable_len_barcodes=False,
                      added_demultiplex_field=None):
    """ Handles all functions for valid data fields

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    errors:  list of errors
    warnings:  list of warnings
    has_barcodes:  If True, will test for perform barcodes test (presence,
     uniqueness, valid IUPAC DNA chars).
    char_replace:  Character used to replace invalid characters in data
     fields.  SampleIDs always use periods to be MIENS compliant.
    variable_len_barcodes:  If True, suppresses warnings about barcodes of
     varying length.
    disable_primer_check:  If True, disables tests for valid primer sequences.
    added_demultiplex_field:  If specified, references a field in the mapping
     file to use for demultiplexing.  These are to be read from fasta labels
     during the actual demultiplexing step.  All combinations of barcodes,
     primers, and the added_demultiplex_field must be unique.
    """

    # Check for valid IUPAC DNA characters in barcode, primer, and reverse
    # primer fields.  Separate check for barcodes and primers, because primer
    # pools separated by commas are allowed, only single barcode allowed.
    # Even if primer check disabled, may still have situation where reverse
    # primers are used, so still do check
    errors = check_dna_chars_primers(header, mapping_data, errors,
                                     disable_primer_check)

    # Skip barcode presence/valid DNA char checks if not barcoded
    if has_barcodes:
        errors = check_dna_chars_bcs(header, mapping_data, errors,
                                     has_barcodes)

    # Check that barcodes all have the same length if not variable length
    if not variable_len_barcodes and has_barcodes:
        warnings = check_bcs_lengths(header, mapping_data, warnings)

    # Check for duplicate barcodes/added demultiplexing fields
    errors = check_bc_duplicates(header, mapping_data, errors, has_barcodes,
                                 variable_len_barcodes, added_demultiplex_field)

    # Check for duplicate SampleIDs
    errors = check_sampleid_duplicates(header, mapping_data, errors)

    # Check for invalid characters
    warnings = check_chars_data_fields(header, mapping_data, warnings)

    # Check for data fields after Description column
    warnings = check_fields_past_bounds(header, mapping_data, warnings)

    # Check for empty data fields before Description column
    warnings = check_empty_fields_before_bounds(header, mapping_data, warnings)

    return errors, warnings


def check_empty_fields_before_bounds(header,
                                     mapping_data,
                                     warnings):
    """ Checks for empty fields before Description header, adds to warnings

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    warnings:  list of warnings
    """

    desc_field = "Description"
    correction = 1
    primer_field = "LinkerPrimerSequence"

    try:
        desc_field_ix = header.index(desc_field) + correction
        primer_field_ix = header.index(primer_field) + correction
    except ValueError:
        # Skip if Description field not present, already get header error
        return warnings

    for curr_row in range(len(mapping_data)):
        for curr_col in range(primer_field_ix, desc_field_ix):
            curr_field = mapping_data[curr_row][curr_col].replace('\n', '')
            if not curr_field:
                warnings.append('Empty data field ' +
                                '%s found\t%d,%d' %
                                (mapping_data[
                                    curr_row][curr_col].replace('\n', ''),
                                 curr_row + correction, curr_col))

    return warnings


def check_fields_past_bounds(header,
                             mapping_data,
                             warnings):
    """ Checks for fields after Description header, adds to warnings

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    warnings:  list of warnings
    """

    desc_field = "Description"
    correction = 1

    try:
        desc_field_ix = header.index(desc_field)
    except ValueError:
        # Skip if Description field not present, already get header error
        return warnings

    for curr_row in range(len(mapping_data)):
        for curr_col in range(len(mapping_data[curr_row])):
            if curr_col > desc_field_ix:
                warnings.append('Data field ' +
                                '%s found after Description column\t%d,%d' %
                                (mapping_data[
                                    curr_row][curr_col].replace('\n', ''),
                                 curr_row + correction, curr_col))

    return warnings


def check_chars_data_fields(header,
                            mapping_data,
                            warnings):
    """ Checks for valid SampleID (MIENS) and other data field characters

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    warnings:  list of warnings
    """

    allowed_data_field_chars = "+-%./ :,;_" + digits + letters
    allowed_sampleid_chars = "." + digits + letters
    correction = 1

    sample_id_field = "SampleID"
    fields_to_skip = ["BarcodeSequence", "LinkerPrimerSequence",
                      "ReversePrimer"]

    for curr_field in range(len(header)):
        if header[curr_field] in fields_to_skip:
            continue
        if header[curr_field] == sample_id_field:
            valid_chars = allowed_sampleid_chars
        else:
            valid_chars = allowed_data_field_chars
        for curr_data in range(len(mapping_data)):
            # Need to skip newline characters
            curr_cell = mapping_data[curr_data][curr_field].replace('\n', '')
            for curr_char in curr_cell:
                if curr_char not in valid_chars:
                    warnings.append("Invalid characters found in %s\t%d,%d" %
                                    (mapping_data[
                                        curr_data][curr_field].replace(
                                        '\n', ''),
                                     curr_data + correction, curr_field))
                    break

    return warnings


def check_dna_chars_primers(header,
                            mapping_data,
                            errors,
                            disable_primer_check=False
                            ):
    """ Checks for valid DNA characters in primer fields

    Also flags empty fields as errors unless flags are passed to suppress
    barcode or primer checks.

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    errors:  list of errors
    disable_primer_check:  If True, disables tests for valid primer sequences.
    """

    valid_dna_chars = DNASequence.iupac_characters()
    valid_dna_chars.add(',')

    # Detect fields directly, in case user does not have fields in proper
    # order in the mapping file (this will generate error separately)
    header_fields_to_check = ["ReversePrimer"]
    if not disable_primer_check:
        header_fields_to_check.append("LinkerPrimerSequence")

    check_indices = []

    for curr_field in range(len(header)):
        if header[curr_field] in header_fields_to_check:
            check_indices.append(curr_field)

    # Correction factor for header being the first line
    correction_ix = 1
    # Check for missing data
    for curr_data in range(len(mapping_data)):
        for curr_ix in check_indices:
            if len(mapping_data[curr_data][curr_ix]) == 0:
                errors.append("Missing expected DNA sequence\t%d,%d" %
                              (curr_data + correction_ix, curr_ix))

    # Check for non-DNA characters
    for curr_data in range(len(mapping_data)):
        for curr_ix in check_indices:
            for curr_nt in mapping_data[curr_data][curr_ix]:
                if curr_nt not in valid_dna_chars:
                    errors.append("Invalid DNA sequence detected: %s\t%d,%d" %
                                  (mapping_data[curr_data][curr_ix],
                                   curr_data + correction_ix, curr_ix))
                    continue

    return errors


def check_dna_chars_bcs(header,
                        mapping_data,
                        errors,
                        has_barcodes=True):
    """ Checks for valid DNA characters in barcode field

    Also flags empty fields as errors unless flags are passed to suppress
    barcode or primer checks.

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    errors:  list of errors
    has_barcodes:  If True, will test for perform barcodes test (presence,
     uniqueness, valid IUPAC DNA chars).
    """

    valid_dna_chars = DNASequence.iupac_standard_characters()
    # Detect fields directly, in case user does not have fields in proper
    # order in the mapping file (this will generate error separately)
    header_fields_to_check = []
    if has_barcodes:
        header_fields_to_check.append("BarcodeSequence")

    check_indices = []

    for curr_field in range(len(header)):
        if header[curr_field] in header_fields_to_check:
            check_indices.append(curr_field)

    # Correction factor for header being the first line
    correction_ix = 1
    # Check for missing data
    for curr_data in range(len(mapping_data)):
        for curr_ix in check_indices:
            if len(mapping_data[curr_data][curr_ix]) == 0:
                errors.append("Missing expected DNA sequence\t%d,%d" %
                              (curr_data + correction_ix, curr_ix))
                continue
            for curr_nt in mapping_data[curr_data][curr_ix]:
                if curr_nt not in valid_dna_chars:
                    errors.append("Invalid DNA sequence detected: %s\t%d,%d" %
                                  (mapping_data[curr_data][curr_ix],
                                   curr_data + correction_ix, curr_ix))
                    continue

    return errors


def check_bcs_lengths(header,
                      mapping_data,
                      warnings):
    """ Adds warnings if barcodes have different lengths

    As this is mostly intended to find typos in barcodes, this will find the
    mode of the barcode lengths, and flag barcodes that are different from
    this.

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    warnings:  list of warnings
    """

    len_counts = defaultdict(int)

    header_field_to_check = "BarcodeSequence"

    # Skip if not field BarcodeSequence
    try:
        check_ix = header.index(header_field_to_check)
    except ValueError:
        return warnings

    for curr_data in range(len(mapping_data)):
        len_counts[len(mapping_data[curr_data][check_ix])] += 1

    # length of the mode
    expected_bc_len = max(len_counts.iteritems(), key=itemgetter(1))[0]

    correction_ix = 1
    for curr_data in range(len(mapping_data)):
        if len(mapping_data[curr_data][check_ix]) != expected_bc_len:
            warnings.append('Barcode %s differs than length %d\t%d,%d' %
                            (mapping_data[
                                curr_data][check_ix], expected_bc_len,
                             curr_data + correction_ix, check_ix))

    return warnings


def check_bc_duplicates(header,
                        mapping_data,
                        errors,
                        has_barcodes=True,
                        variable_len_barcodes=False,
                        added_demultiplex_field=None):
    """ Checks for barcode and other demultiplexing duplicates

    Default check is for unique barcodes.  A potential tricky situation to
    handle is variable length barcodes, than when combined with the 5' end of
    a primer sequence, are indistinguishable, and these are tested for as well.
    Finally, combinations of barcode sequences and added_demultiplex_field
    values are tested to ensure that combinations of these values are unique.

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    errors:  list of errors
    has_barcodes:  If True, will test for perform barcodes test (presence,
     uniqueness, valid IUPAC DNA chars).
    variable_len_barcodes:  If True, suppresses warnings about barcodes of
     varying length.
    added_demultiplex_field:  If specified, references a field in the mapping
     file to use for demultiplexing.  These are to be read from fasta labels
     during the actual demultiplexing step.  All combinations of barcodes,
     primers, and the added_demultiplex_field must be unique.
    """

    if (has_barcodes and not variable_len_barcodes
            and not added_demultiplex_field):
        errors = check_fixed_len_bcs_dups(header, mapping_data, errors)
    if (has_barcodes and variable_len_barcodes
            and not added_demultiplex_field):
        errors = check_variable_len_bcs_dups(header, mapping_data, errors)
    if added_demultiplex_field:
        errors = check_added_demultiplex_dups(header, mapping_data, errors,
                                              has_barcodes, added_demultiplex_field)
    # Special case of has_barcodes = False and no added_demultiplex_field,
    # need to check that only a single SampleID is passed in this case so
    # we have "unique" demultiplexing.
    if (not has_barcodes and not added_demultiplex_field):
        # only one line of mapping data for one sample
        if len(mapping_data) != 1:
            errors.append("If no barcodes are present, and the " +
                          "added_demultiplex_field option isn't used, only a single " +
                          "SampleID can be present.\t-1,-1")
    return errors


def check_fixed_len_bcs_dups(header,
                             mapping_data,
                             errors):
    """ Checks barcodes of same length for duplicates, adds to errors if found

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    errors:  list of errors
    """

    header_field_to_check = "BarcodeSequence"

    # Skip if no field BarcodeSequence
    try:
        check_ix = header.index(header_field_to_check)
    except ValueError:
        return errors

    barcodes = []

    correction = 1

    for curr_data in mapping_data:
        barcodes.append(upper(curr_data[check_ix]))

    dups = duplicates_indices(barcodes)

    for curr_dup in dups:
        for curr_loc in dups[curr_dup]:
            errors.append('Duplicate barcode %s found.\t%d,%d' %
                          (curr_dup, curr_loc + correction, check_ix))

    return errors


def check_variable_len_bcs_dups(header,
                                mapping_data,
                                errors):
    """ Checks variable length barcodes plus sections of primers for dups

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    errors:  list of errors
    """

    header_field_to_check = "BarcodeSequence"

    # Skip if no field BarcodeSequence
    try:
        check_ix = header.index(header_field_to_check)
    except ValueError:
        return errors

    linker_primer_field = "LinkerPrimerSequence"

    try:
        linker_primer_ix = header.index(linker_primer_field)
        no_primers = False
    except ValueError:
        no_primers = True

    barcodes = []
    bc_lens = []

    correction = 1

    for curr_data in mapping_data:
        barcodes.append(upper(curr_data[check_ix]))
        bc_lens.append(len(curr_data[check_ix]))

    # Get max length of barcodes to determine how many primer bases to slice
    barcode_max_len = max(bc_lens)

    # Have to do second pass to append correct number of nucleotides to
    # check for duplicates between barcodes and primer sequences

    bcs_added_nts = []
    for curr_data in mapping_data:
        if no_primers:
            bcs_added_nts.append(upper(curr_data[check_ix]))
        else:
            adjusted_len = barcode_max_len - len(curr_data[check_ix])
            bcs_added_nts.append(upper(curr_data[check_ix] +
                                       curr_data[linker_primer_ix][0:adjusted_len]))

    dups = duplicates_indices(bcs_added_nts)

    for curr_dup in dups:
        for curr_loc in dups[curr_dup]:
            if no_primers:
                errors.append('Duplicate barcode %s found.\t%d,%d' %
                              (curr_dup, curr_loc + correction, check_ix))
            else:
                errors.append('Duplicate barcode and primer fragment sequence ' +
                              '%s found.\t%d,%d' % (curr_dup, curr_loc + correction, check_ix))

    return errors


def check_added_demultiplex_dups(header,
                                 mapping_data,
                                 errors,
                                 has_barcodes=True,
                                 added_demultiplex_field=None):
    """ Checks that all barcodes and added demultiplex fields are unique

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    errors:  list of errors
    has_barcodes:  True if barcode fields are to be used.
    added_demultiplex_field:  If specified, references a field in the mapping
     file to use for demultiplexing.  These are to be read from fasta labels
     during the actual demultiplexing step.  All combinations of barcodes,
     primers, and the added_demultiplex_field must be unique.
    """

    # Treat as variable length to test combinations of barcodes and the
    # added demultiplex field (should return the same result for the barcode
    # component)
    correction = 1

    header_field_to_check = "BarcodeSequence"
    bc_found = False

    # Skip if no field BarcodeSequence
    if has_barcodes:
        try:
            bc_ix = header.index(header_field_to_check)
            bc_found = True
        except ValueError:
            pass

    linker_primer_field = "LinkerPrimerSequence"

    try:
        linker_primer_ix = header.index(linker_primer_field)
        no_primers = False
    except ValueError:
        no_primers = True

    try:
        added_demultiplex_ix = header.index(added_demultiplex_field)
    except ValueError:
        # Skip out at this point, header check will have error for missing
        # field
        return errors

    barcodes = []
    bc_lens = []
    bcs_added_field = []

    if has_barcodes and bc_found:
        for curr_data in mapping_data:
            barcodes.append(upper(curr_data[bc_ix]))
            bc_lens.append(len(curr_data[bc_ix]))

        # Get max length of barcodes to determine how many primer bases to
        # slice
        barcode_max_len = max(bc_lens)

        # Have to do second pass to append correct number of nucleotides to
        # check for duplicates between barcodes and primer sequences

        for curr_data in mapping_data:
            if no_primers:
                bcs_added_field.append(curr_data[bc_ix] +
                                       curr_data[added_demultiplex_ix])
            else:
                adjusted_len = barcode_max_len - len(curr_data[bc_ix])
                bcs_added_field.append(curr_data[bc_ix] +
                                       curr_data[linker_primer_ix][0:adjusted_len] +
                                       curr_data[added_demultiplex_ix])
    else:
        for curr_data in mapping_data:
            bcs_added_field.append(curr_data[added_demultiplex_ix])

    dups = duplicates_indices(bcs_added_field)

    for curr_dup in dups:
        if has_barcodes and bc_found:
            for curr_loc in dups[curr_dup]:
                errors.append('Duplicate barcode and added demultiplex field ' +
                              '%s found.\t%d,%d' % (curr_dup, curr_loc + correction, bc_ix))
        else:
            for curr_loc in dups[curr_dup]:
                errors.append('Duplicate added demultiplex field ' +
                              '%s found.\t%d,%d' % (curr_dup, curr_loc + correction,
                                                    added_demultiplex_ix))

    return errors


def check_sampleid_duplicates(header,
                              mapping_data,
                              errors):
    """ Flags duplicate, missing SampleIDs as errors

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    errors:  list of errors
    """

    sample_id_field = "SampleID"
    correction = 1

    try:
        sample_id_ix = header.index(sample_id_field)
    except ValueError:
        # Skip out at this point, header check will have error for missing
        # field
        return errors

    sample_ids = []

    # Need to save locations of missing IDs so they aren't flagged twice
    missing_sample_ids = []

    for curr_data in range(len(mapping_data)):
        if len(mapping_data[curr_data][sample_id_ix]) == 0:
            errors.append('Missing SampleID.\t%d,%d' %
                          (curr_data + correction, sample_id_ix))
            missing_sample_ids.append(curr_data + correction)
        sample_ids.append(mapping_data[curr_data][sample_id_ix])

    dups = duplicates_indices(sample_ids)

    for curr_dup in dups:
        for curr_loc in dups[curr_dup]:
            if (curr_loc + correction) not in missing_sample_ids:
                errors.append('Duplicate SampleID %s found.\t%d,%d' %
                              (curr_dup, curr_loc + correction, sample_id_ix))

    return errors


# End data field checking functions
# Begin header field checking functions
def check_header(header,
                 errors,
                 warnings,
                 sample_id_ix,
                 desc_ix,
                 bc_ix,
                 linker_primer_ix,
                 added_demultiplex_field=None):
    """ Checks header for valid characters, unique and required fields

    header:  list of header strings
    errors:  list of errors
    warnings:  list of warnings
    sample_id_ix:  index of SampleID in header
    desc_ix: index of Description in header
    bc_ix:  index of BarcodeSequence in header
    linker_primer_ix:  index of LinkerPrimerSequence in header
    added_demultiplex_field:  If specified, references a field in the mapping
     file to use for demultiplexing.  These are to be read from fasta labels
     during the actual demultiplexing step.  All combinations of barcodes,
     primers, and the added_demultiplex_field must be unique.
    """

    # Check for duplicates, append to errors if found
    errors = check_header_dups(header, errors)

    # Check for valid characters
    warnings = check_header_chars(header, warnings)

    # Check for required header fields
    errors = check_header_required_fields(header, errors, sample_id_ix,
                                          desc_ix, bc_ix, linker_primer_ix, added_demultiplex_field)

    return errors, warnings


def check_header_dups(header,
                      errors):
    """ Checks for duplicates in headers, appends to errors if found

    header:  list of header strings
    errors:  list of errors
    """

    for curr_elem in range(len(header)):
        if header.count(header[curr_elem]) != 1:
            errors.append('%s found in header %d times.  ' %
                          (header[curr_elem], header.count(header[curr_elem])) +
                          'Header fields must be unique.\t%d,%d' % (0, curr_elem))

    return errors


def check_header_chars(header,
                       warnings,
                       allowed_chars_header='_' + digits + letters):
    """ Checks for valid characters in headers, appends to warnings

    header:  list of header strings
    warnings:  list of warnings
    """

    for curr_elem in range(len(header)):
        for curr_char in header[curr_elem]:
            if curr_char not in allowed_chars_header:
                warnings.append('Found invalid character in %s ' %
                                header[curr_elem] + 'header field.\t%d,%d' % (0, curr_elem))
                break

    return warnings


def check_header_required_fields(header,
                                 errors,
                                 sample_id_ix,
                                 desc_ix,
                                 bc_ix,
                                 linker_primer_ix,
                                 added_demultiplex_field=None):
    """ Checks for required header fields, appends to errors if not found

    header:  list of header strings
    errors:  list of errors
    sample_id_ix:  index of SampleID in header
    desc_ix: index of Description in header
    bc_ix:  index of BarcodeSequence in header
    linker_primer_ix:  index of LinkerPrimerSequence in header
    """

    header_checks = {
        sample_id_ix: "SampleID",
        desc_ix: "Description",
        bc_ix: "BarcodeSequence",
        linker_primer_ix: "LinkerPrimerSequence"
    }

    for curr_check in header_checks:
        if (header[curr_check] != header_checks[curr_check] and
                header_checks[curr_check] == "Description"):
            errors.append('Found header field %s, last field should be %s' %
                          (header[curr_check], header_checks[curr_check]) +
                          '\t%d,%d' % (0, curr_check))
        elif (header[curr_check] != header_checks[curr_check] and
              header_checks[curr_check] != "Description"):
            errors.append('Found header field %s, expected field %s' %
                          (header[curr_check], header_checks[curr_check]) +
                          '\t%d,%d' % (0, curr_check))

    if added_demultiplex_field:
        if added_demultiplex_field not in header:
            errors.append('Missing added demultiplex field %s\t%d,%d' %
                          (added_demultiplex_field, -1, -1))

    return errors

# End header field checking functions

# Misc functions


def correct_mapping_data(mapping_data,
                         header,
                         char_replace="_"):
    """ Replaces invalid characters in mapping data

    mapping_data:  list of lists of raw metadata mapping file data
    header:  list of header strings
    char_replace:  Character used to replace invalid characters in data
     fields.  SampleIDs always use periods to be MIENS compliant.
    """

    corrected_data = deepcopy(mapping_data)

    valid_sample_id_chars = letters + digits + "."
    valid_data_field_chars = letters + digits + "+-%./ :,;_"

    sample_id_char_replace = "."

    sample_id_field = "SampleID"
    fields_to_skip = ["BarcodeSequence", "LinkerPrimerSequence",
                      "ReversePrimer"]

    try:
        sample_id_ix = header.index(sample_id_field)
    except ValueError:
        sample_id_ix = -1

    fields_to_skip_ixs = []
    for curr_field in fields_to_skip:
        try:
            fields_to_skip_ixs.append(header.index(curr_field))
        except ValueError:
            continue

    for curr_row in range(len(mapping_data)):
        for curr_col in range(len(mapping_data[curr_row])):
            if curr_col in fields_to_skip_ixs:
                continue
            elif (sample_id_ix != -1) and (curr_col == sample_id_ix):
                curr_replacement = sample_id_char_replace
                curr_valid_chars = valid_sample_id_chars
            else:
                curr_replacement = char_replace
                curr_valid_chars = valid_data_field_chars

            curr_corrected_field = ""
            for curr_char in mapping_data[curr_row][curr_col].replace('\n', ''):
                if curr_char not in curr_valid_chars:
                    curr_corrected_field += curr_replacement
                else:
                    curr_corrected_field += curr_char

            corrected_data[curr_row][curr_col] = curr_corrected_field

    return corrected_data


def write_corrected_mapping(output_corrected_fp,
                            header,
                            run_description,
                            corrected_mapping_data):
    """ Writes corrected mapping file with invalid characters replaced

    output_corrected_fp:  Filepath to write corrected mapping file to.
    header:  list of strings of header data
    run_description:  Comment lines, written after header.
    corrected_mapping_data:  list of lists of corrected mapping data.
    """

    out_f = open(output_corrected_fp, "w")

    out_f.write("#" + "\t".join(header).replace('\n', '') + "\n")

    for curr_comment in run_description:
        out_f.write("#" + curr_comment.replace('\n', '') + "\n")

    for curr_data in corrected_mapping_data:
        out_f.write("\t".join(curr_data).replace('\n', '') + "\n")


def write_log_file(output_log_fp,
                   errors,
                   warnings):
    """ Writes log file with details of errors, warnings in mapping file.

    output_log_fp:  output filepath for log file.
    errors:  list of errors
    warnings:  list of warnings
    """

    out_f = open(output_log_fp, "w")

    if not errors and not warnings:
        out_f.write("No errors or warnings found in mapping file.")
        return

    out_f.write("# Errors and warnings are written as a tab separated " +
                "columns, with the first column showing the error or warning, and the " +
                "second column contains the location of the error or warning, written " +
                "as row,column, where 0,0 is the top left header item (SampleID).  " +
                "Problems not specific to a particular data cell will be listed as " +
                "having 'no location'.\n")
    out_f.write("Errors -----------------------------\n")
    for error in errors:
        if error.split('\t')[1] == "-1,-1":
            curr_err = error.split('\t')[0] + "\tno location"
        else:
            curr_err = error
        out_f.write(curr_err + "\n")
    out_f.write("Warnings ---------------------------\n")
    for warning in warnings:
        if warning.split('\t')[1] == "-1,-1":
            curr_warning = warning.split('\t')[0] + "\tno location"
        else:
            curr_warning = warning
        out_f.write(curr_warning + "\n")
