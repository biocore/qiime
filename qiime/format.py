#!/usr/bin/env python
from __future__ import division
__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski", "Jeremy Widmann",
               "Antonio Gonzalez Pena", "Daniel McDonald", "Jai Ram Rideout"]
# remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import numpy
from numpy import asarray, isnan, log10, median
from StringIO import StringIO
from re import compile, sub
from os import walk
from os.path import join, splitext, exists, isfile, abspath

from skbio.sequence import BiologicalSequence
from biom.table import Table

from qiime.util import get_qiime_library_version, load_qiime_config
from qiime.colors import data_color_hsv

"""Contains formatters for the files we expect to encounter in 454 workflow.

A lot of this might migrate into cogent at some point.
"""


def format_mapping_file(headers, mapping_data, comments=None):
    """ returns a large formatted string representing the entire mapping file

    each input is a list, and all data should be strings, not e.g. ints
    * headers defines column labels, and SampleID should not include a '#'
    * mapping_data is a list of lists, each sublist is a row in the mapping file
    each mapping_data sublist must be the same length as headers - use ''
    for absent data
    * if included, commments will be inserted above the header line
    comments should not include a # - that will be appended in this formatter
    """

    result = []  # each elem is a string representing a line

    result.append('#' + '\t'.join(headers))

    if comments is not None:
        for comment in comments:
            result.append('#' + comment)

    for mapping_line in mapping_data:
        if not (len(mapping_line) == len(headers)):
            raise RuntimeError('error formatting mapping file, does each ' +
                               'sample have the same length of data as the headers?')
        result.append('\t'.join(mapping_line))

    str_result = '\n'.join(result)
    return str_result


def format_p_value_for_num_iters(p, num_iters):
    """adjust p to str w correct num of decimals for num monte carlo iters
    """
    if num_iters < 10:
        # this can be the last step of a long process, so we don't
        # want to fail
        return "Too few iters to compute p-value (num_iters=%d)" % num_iters
    decimal_places = int(log10(num_iters + 1))
    result = ('%1.' + '%df' % decimal_places) % p
    return result


def format_qiime_parameters(params, header="#QIIME parameters"):
    """Formats lines for qiime_parameters.txt"""
    qiime_params = [header]
    for script, options in sorted(params.items()):
        for option, value in sorted(options.items()):
            specific_option = ':'.join([script, option])

            # Based on how qiime_parameters is parsed
            if value is None:
                value = "True"

            # cast value to string just in case
            full_line = '\t'.join([specific_option, str(value)])

            qiime_params.append(full_line)
    return qiime_params


def format_add_taxa_summary_mapping(summary, tax_order, mapping, header,
                                    delimiter=';'):
    """Formats a summarized taxonomy with mapping information"""
    tax_order = [delimiter.join(tax) for tax in tax_order]
    header.extend(tax_order)
    yield "#%s\n" % '\t'.join(header)

    for row in mapping:
        sample_id = row[0]

        # only save samples we have summaries for
        if sample_id not in summary:
            continue

        # grab otu counts for each taxon
        row.extend(map(str, summary[sample_id]))
        yield "%s\n" % '\t'.join(row)


def write_add_taxa_summary_mapping(summary, tax_order, mapping, header,
                                   output_fp, delimiter=';'):
    """ """
    of = open(output_fp, 'w')
    for line in format_add_taxa_summary_mapping(summary, tax_order, mapping,
                                                header, delimiter):
        of.write(line)
    of.close()


def format_taxa_summary(taxa_summary):
    """Formats a taxa summary to be suitable for writing to a file.

    Returns a string where the first line contains 'Taxon' and then a list of
    tab-separated sample IDs. The next lines contain each taxon and its
    abundance in each sample separated by a tab.

    Arguments:
        taxa_summary - the taxa summary tuple (i.e. the output of
            qiime.parse.parse_taxa_summary_table) to format
    """
    result = 'Taxon\t' + '\t'.join(taxa_summary[0]) + '\n'
    for taxon, row in zip(taxa_summary[1], taxa_summary[2]):
        row = map(str, row)
        result += '%s\t' % taxon + '\t'.join(row) + '\n'
    return result


def format_correlation_vector(corr_vector, num_permutations, header=''):
    """Formats a correlation vector to be suitable for writing to a file.

    Returns a string where each line contains nine tab-separated fields: the
    two sample IDs that were compared, the correlation coefficient, the
    parametric p-value, the Bonferroni-corrected parametric p-value, the
    nonparametric p-value, the Bonferroni-corrected nonparametric p-value, and
    the confidence interval (lower and upper endpoints as separate fields).

    If the confidence intervals are not valid for this dataset (i.e. the
    input CI is (None, None)), the confidence intervals will be formatted as
    'N/A' for both lower and upper endpoints.

    Arguments:
        corr_vector - a list of 8-element tuples, where the first
            element is a sample ID, the second element is a sample ID, the
            third element is the correlation coefficient computed between the
            two samples (a double), the fourth element is the parametric
            p-value, the fifth value is the Bonferroni-corrected parametric
            p-value, the sixth value is the nonparametric p-value, the seventh
            value is the Bonferroni-corrected nonparametric p-value, and the
            eighth element is a tuple containing the low and high ends of the
            confidence interval
        num_permutations - the number of permutations that were used to
            calculate the nonparametric p-values. Will be used to format the
            correct number of digits for these p-values. If less than 1, the
            p-values will be 'N/A'
        header - if provided, this string will be inserted at the beginning of
            the returned string. For example, might be useful to add a comment
            describing what correlation coefficient was used. This string does
            not need to contain a newline at the end
    """
    result = ''
    if header != '':
        result += header + '\n'
    result += 'Sample ID\tSample ID\tCorrelation coefficient\t' + \
              'Parametric p-value\tParametric p-value ' + \
              '(Bonferroni-corrected)\tNonparametric p-value\t' + \
              'Nonparametric p-value (Bonferroni-corrected)\t' + \
              'CI (lower)\tCI (upper)\n'
    for samp_id1, samp_id2, corr_coeff, param_p_val, param_p_val_corr, \
            nonparam_p_val, nonparam_p_val_corr, conf_interval in corr_vector:
        if num_permutations > 0:
            nonparam_p_val_str = format_p_value_for_num_iters(nonparam_p_val,
                                                              num_permutations)
            nonparam_p_val_corr_str = format_p_value_for_num_iters(
                nonparam_p_val_corr, num_permutations)
        else:
            nonparam_p_val_str = 'N/A'
            nonparam_p_val_corr_str = 'N/A'
        if conf_interval == (None, None):
            conf_interval_str = 'N/A\tN/A'
        else:
            conf_interval_str = '%.4f\t%.4f' % conf_interval

        result += '%s\t%s\t%.4f\t%.4f\t%.4f\t%s\t%s\t%s\n' % (
            samp_id1, samp_id2, corr_coeff, param_p_val, param_p_val_corr,
            nonparam_p_val_str, nonparam_p_val_corr_str, conf_interval_str)
    return result


def format_correlation_info(corr_coeff, param_p_val, nonparam_p_val,
                            conf_interval, num_permutations, header=''):
    """Formats correlation information to be suitable for writing to a file.

    Returns a string containing a header and a single line (with a newline at
    the end) that has the input correlation information in tab-separated
    format, with nonparametric p-value formatted according to the number of
    permutations.

    If the confidence interval is not valid for this dataset (i.e. the
    input CI is (None, None)), the confidence interval will be formatted as
    'N/A' for both lower and upper endpoints.

    Arguments:
        corr_coeff - the correlation coefficient (a float)
        param_p_val - the parametric p-value (a float)
        nonparam_p_val - the nonparametric p-value (a float)
        conf_interval - a tuple containing the lower and upper bounds of the
            confidence interval
        num_permutations - the number of permutations that were used to
            calculate the nonparametric p-value. Will be used to format the
            correct number of digits for this p-value. If less than 1, the
            p-value will be 'N/A'
        header - if provided, this string will be inserted at the beginning of
            the returned string. For example, might be useful to add a comment
            describing what correlation coefficient was used. This string does
            not need to contain a newline at the end
    """
    result = ''
    if header != '':
        result += header + '\n'
    result += 'Correlation coefficient\tParametric p-value\t' + \
              'Nonparametric p-value\tCI (lower)\tCI (upper)\n'

    if num_permutations > 0:
        nonparam_p_val_str = format_p_value_for_num_iters(nonparam_p_val,
                                                          num_permutations)
    else:
        nonparam_p_val_str = 'N/A'

    if conf_interval == (None, None):
        conf_interval_str = 'N/A\tN/A'
    else:
        conf_interval_str = '%.4f\t%.4f' % conf_interval

    result += '%.4f\t%.4f\t%s\t%s\n' % (corr_coeff, param_p_val,
                                        nonparam_p_val_str, conf_interval_str)
    return result


def write_otu_map(otu_map, output_fp, otu_id_prefix=''):
    """
    """
    of = open(output_fp, 'w')
    for line in format_otu_map(otu_map, otu_id_prefix):
        of.write(line)
    of.close()


def format_otu_map(otu_map, otu_id_prefix):
    """ Takes list of format [(otu_id,[seq_ids]), ... ]
    """
    # raise error if prefix contains chars other than . or alnums
    # this functionality needs to be centralized
    for c in otu_id_prefix:
        if not c.isalnum() and not c == '.':
            raise ValueError("%s char is not allowed in OTU IDs" % c)

    for otu_id, seq_ids in otu_map:
        yield '%s%s\t%s\n' % (otu_id_prefix,
                              otu_id,
                              '\t'.join(seq_ids))
    return
format_observation_map = format_otu_map


def format_distance_matrix(labels, data):
    """Writes distance matrix as tab-delimited text, uses format_matrix"""
    return format_matrix(data, labels, labels,
                         convert_matching_names_to_zero=True)


def format_matrix(data, row_names, col_names,
                  convert_matching_names_to_zero=False):
    """Writes matrix as tab-delimited text.

    format is rows: samples,
    cols: metrics/ confidence levels, etc

    data: array or 2d list.

    convert_matching_names_to_zero : bool
        If ``True``, each element in `data` with a matching row name and column
        name will be converted to 0.0 if the original value is close to zero
        (as determined by ``numpy.isclose(value, 0.0)``). This is a necessary
        evil to handle distance matrices with diagonals that are close to zero.
        We have to perform this conversion here instead of simply checking the
        distance matrix's diagonal because parallelized beta diversity
        calculations will compute subsets of the entire distance matrix and we
        don't know where the true diagonal will fall. See
        https://github.com/biocore/qiime/issues/1933 for details.

    """
    len_col = len(col_names)
    try:
        if data.shape != (len(row_names), len_col):
            raise ValueError(
                "Data shape of %s doesn't match header sizes %s %s" %
                (data.shape, len(row_names), len(col_names)))
    except AttributeError:
        # must be list of list
        try:
            if not numpy.all([len_col == len(row) for row in data]) or\
                    len(row_names) != len(data):
                raise ValueError(
                    "Data shape doesn't match header sizes %s %s" %
                    (len(row_names), len(col_names)))
        except:
            raise ValueError("Unsupported data type for format_matrix")

    lines = []
    row_names = map(str, row_names)
    col_names = map(str, col_names)
    # just in case they weren't strings initially
    lines.append('\t'.join([''] + col_names))
    for row_name, values in zip(row_names, data):
        line = [row_name]
        for col_name, value in zip(col_names, values):
            if convert_matching_names_to_zero and col_name == row_name and \
                    numpy.isclose(value, 0.0):
                value = 0.0
            line.append(str(value))
        lines.append('\t'.join(line))
    return '\n'.join(lines)


def format_nmds_coords(samples, points, stress):
    """ samples is list, points is samples by axis coord (typ many by 2 mtx)
    """
    result = []
    col_headers = ["NMDS" + str(aa) for aa in range(1, points.shape[1] + 1)]
    result.append('samples\t' + '\t'.join(col_headers))
    for name, row in zip(samples, points):
        result.append('\t'.join([name] + map(str, row)))
    result.append('')
    result.append('stress\t' + str(stress) + '\t' +
                  '0\t' * (len(col_headers) - 2) + '0')
    result.append('% variation explained\t' + '0\t' * (len(col_headers)))
    return '\n'.join(result)


def build_prefs_string(mapping_headers_to_use, background_color,
                       monte_carlo_dist, headers, otu_ids, ball_scale,
                       arrow_line_color, arrow_head_color):
    """Create a preferences file, which can be used for some of the \
    visualization scripts."""

    # Open up the prefs dictionary
    pref_lines = ["{\n"]

    # Define and add the background_color dictionary to prefs dictionary
    bk_color = "'background_color':'%s',\n" % (background_color)
    pref_lines.append(bk_color)

    # create a unique field dictionary for use with the FIELDS dictionary
    unique_dist_fields = {}

    # Iterate through the user-supplied fields or all the fields in the mapping
    # file, then validate that the fields exist
    if mapping_headers_to_use == 'ALL':
        fields = headers
        for key in fields:
            unique_dist_fields[key] = monte_carlo_dist
    else:
        # If '&&' is used, split into multiple fields and valid them against
        # the mapping file
        fields = mapping_headers_to_use.split(',')
        for field_id in fields:
            f_str = field_id.split('&&')
            for id_ in f_str:
                validity = False
                for head_id in headers:
                    if id_ == head_id:
                        unique_dist_fields[id_] = monte_carlo_dist
                        validity = True
                if not validity:
                    raise ValueError(
                        "%s is not a header in your mapping file" %
                        field_id)

    # Syntax for sample_coloring dictionary
    sample_coloring = ["\n'sample_coloring':\n\t{"]
    sample_colors = \
        "\t\t'%s':" + \
        "\n\t\t{" + \
        "\n\t\t\t'column':'%s'," + \
        "\n\t\t\t'colors':(('red',(0,100,100)),('blue',(240,100,100)))" + \
        "\n\t\t}"

    # Syntax for monte_carlo dictionary
    monte_carlo_main = ["\n'MONTE_CARLO_GROUP_DISTANCES':\n\t{\n"]
    monte_carlo = []
    monte_carlo_distances = "\t\t'%s': %s"

    # Syntax for fields dictionary
    field_dict_main = ["\n'FIELDS':\n\t[\n"]
    field_dict = []
    dist_fields = "\t\t'%s'"

    # This iterates through the fields and creates a sample_color dictionary
    # values for each field
    first = True
    for field in fields:
        if first:
            first = False
            sample_coloring.append('\n')
        else:
            sample_coloring.append(',\n')
        sample_coloring.append(sample_colors % (field, field))
        monte_carlo.append(monte_carlo_distances % (field, monte_carlo_dist))

    # Syntax for taxonomy_coloring dictionary
    taxon_start = \
        "\n'taxonomy_coloring':\n\t{\n" + \
        "\t\t'Level_%s':" + \
        "\n\t\t{" + \
        "\n\t\t\t'column':'%s'," + \
        "\n\t\t\t'colors':\n\t\t\t{"
    taxon_colors = "\n\t\t\t\t'%s':('red%s',(%d,100,100))"
    taxon_coloring = []

    if otu_ids:
        level = max([len(t.split(';')) - 1 for t in otu_ids])
        taxon_coloring.append(taxon_start % (str(level), str(level)))
        taxons = []
        otu_id_iter = (240.0 / (len(otu_ids) + 1))
        counter = 0
        for i in otu_ids:
            taxons.append(taxon_colors % (i, str(counter), counter))
            counter = counter + otu_id_iter
        taxon_coloring.append(','.join(taxons))
    else:
        taxon_coloring.append(taxon_start % (str(1), str(1)))
        taxon_coloring.append(taxon_colors % ('Root;Bacteria', str(0), 0))

    taxon_coloring.append("\n\t\t\t}\n\t\t}\n\t},\n")
    taxonomy_coloring_str = ''.join(taxon_coloring)

    # Close and convert the sample_coloring dictionary to a string
    sample_coloring.append('\n\t},')
    sample_coloring_str = ''.join(sample_coloring)

    # Close and convert the monte_carlo dictionary to a string
    monte_carlo1 = ',\n'.join(monte_carlo)
    monte_carlo_main.append(monte_carlo1)
    monte_carlo_main.append('\n\t},')
    monte_carlo_str = ''.join(monte_carlo_main)

    # This iterates through the fields and creates the monte_carlo and fields
    # dictionary values
    for field in unique_dist_fields:
        field_dict.append(dist_fields % (field))

    # Close and convert the fields dictionary to a string
    field_dict1 = ',\n'.join(field_dict)
    field_dict_main.append(field_dict1)
    field_dict_main.append('\n\t],')
    field_dict_str = ''.join(field_dict_main)

    # Add all the dictionary values to the prefs dictionary
    pref_lines.append(sample_coloring_str)
    pref_lines.append(monte_carlo_str)
    pref_lines.append(field_dict_str)
    pref_lines.append(taxonomy_coloring_str)

    # Add ball_scale
    b_scale = "'ball_scale':'%f',\n" % (ball_scale)
    pref_lines.append(b_scale)

    # Add arrow_line_color
    alc = "'arrow_line_color':'%s',\n" % (arrow_line_color)
    pref_lines.append(alc)

    # Add arrow_head_color
    ahc = "'arrow_head_color':'%s'" % (arrow_head_color)
    pref_lines.append(ahc)

    # Closing tabs
    pref_lines.append('\n}')

    return ''.join(pref_lines)


def format_map_file(headers, id_map, desc_key, sample_id_key,
                    description_map=None, run_description=None):
    """Generates string for formatted map file.

    Input:
        headers: list of strings corresponding to col headers
        id_map: dict of {id:{header:val}}
        description_map: dict of {id:description}
        run_description: either string, or list of strings
    """
    result = []
    if desc_key in headers:
        headers.remove(desc_key)
    if sample_id_key in headers:
        headers.remove(sample_id_key)
    header_line = '\t'.join([sample_id_key] + headers + [desc_key])
    if not header_line.startswith('#'):
        header_line = '#' + header_line
    result.append(header_line)
    if run_description:
        if not isinstance(run_description, str):
            run_description = '\n#'.join(run_description)
        if not run_description.startswith('#'):
            run_description = '#' + run_description
        result.append(run_description)
    for id_, fields in sorted(id_map.items()):
        curr_line = [id_]
        curr_line.extend([fields.get(h, '') for h in headers])
        curr_line.append(description_map.get(id_, ''))
        result.append('\t'.join(map(str, curr_line)))
    return '\n'.join(result)


def format_histograms_two_bins(pre_hist, post_hist, bin_edges):
    """Returns text-formatted histogram."""
    lines = []
    lines.append('Length\tBefore\tAfter')
    for edge, pre, post in zip(bin_edges, pre_hist, post_hist):
        lines.append('\t'.join(map(str, [edge, pre, post])))
    return '\n'.join(lines)


def format_histograms(raw_hist, pre_hist, post_hist, bin_edges):
    """Returns text-formatted histogram.  Needs to take 3 bins of data"""
    lines = []
    lines.append('# bins raw sequence lengths, length of sequences that '
                 'pass quality filters before processing, and lengths of sequences that '
                 'pass quality filters post processing.')
    lines.append('Length\tRaw\tBefore\tAfter')
    for edge, raw, pre, post in zip(bin_edges, raw_hist, pre_hist, post_hist):
        lines.append('\t'.join(map(str, [edge, raw, pre, post])))
    return '\n'.join(lines)


def format_histogram_one_count(counts, bin_edges):
    """Returns text-formatted histogram with only one count."""
    lines = []
    lines.append('Length\tCount')
    for edge, count in zip(bin_edges, counts):
        lines.append('\t'.join(map(str, [edge, count])))
    return '\n'.join(lines)


def format_split_libraries_fastq_log(count_barcode_not_in_map,
                                     count_too_short,
                                     count_too_many_N,
                                     count_bad_illumina_qual_digit,
                                     count_barcode_errors_exceed_max,
                                     input_sequence_count,
                                     sequence_lengths,
                                     seqs_per_sample_counts):
    """ Format the split libraries log """
    log_out = ["Quality filter results"]
    log_out.append(
        "Total number of input sequences: %d" %
        input_sequence_count)
    log_out.append(
        "Barcode not in mapping file: %d" %
        count_barcode_not_in_map)
    log_out.append(
        "Read too short after quality truncation: %d" %
        count_too_short)
    log_out.append(
        "Count of N characters exceeds limit: %d" %
        count_too_many_N)
    log_out.append(
        "Illumina quality digit = 0: %d" %
        count_bad_illumina_qual_digit)
    log_out.append(
        "Barcode errors exceed max: %d" %
        count_barcode_errors_exceed_max)

    log_out.append("")

    log_out.append("Result summary (after quality filtering)")
    log_out.append("Median sequence length: %1.2f" % median(sequence_lengths))
    counts = sorted([(v, k) for k, v in seqs_per_sample_counts.items()])
    counts.reverse()
    for sequence_count, sample_id in counts:
        log_out.append('%s\t%d' % (sample_id, sequence_count))

    total_seqs_written = 0
    for curr_count in counts:
        total_seqs_written += curr_count[0]

    log_out.append('\nTotal number seqs written\t%d' % total_seqs_written)
    return '\n'.join(log_out)


def format_unifrac_sample_mapping(sample_ids, otu_ids, otu_table_array):
    """Returns a unifrac sample mapping file from output of parse_otu_table
    """
    out = []
    for i, row in enumerate(otu_table_array):
        for j, val in enumerate(row):
            if val > 0:
                line = [otu_ids[i], sample_ids[j], str(val)]
                out.append('\t'.join(line))
    return out


def write_Fasta_from_name_seq_pairs(name_seqs, fh):
    """writes a list of (name,seqs) to filehandle.

    name_seqs: (name,seqs) pair such as from parse_fasta
    fh: an open filehandle
    """
    if fh is None:
        raise ValueError("Need open file handle to write to.")

    for (name, seq) in name_seqs:
        fh.write("%s\n" % BiologicalSequence(seq, id=name).to_fasta())


def illumina_data_to_fastq(record_data, number_of_bases=None):
    """ given data from an Illumina qseq file, write to fastq

        read data: generator of tuples with the following data
         (machine ID,
          ... [GREG: NEED TO FILL IN DETAILS!]
          illumina quality digit (0:failed screen; 1: passed screen)
          read number,
          sequence,
          quality string)

        number_of_bases: number of bases to keep, starting from
         beginnng of the read - useful when additional cycles are
         applied (e.g., sometimes happens when sequencing barcodes)

    """
    seq_index = 8
    qual_index = 9
    pass_filter_index = 10

    try:
        pass_filter = int(record_data[pass_filter_index])
    except IndexError:
        pass_filter = 2

    if number_of_bases is None:
        seq = record_data[seq_index].replace('.', 'N')
        qual = record_data[qual_index]
    else:
        seq = record_data[seq_index][:number_of_bases].replace('.', 'N')
        qual = record_data[qual_index][:number_of_bases]

    header = '%s_%s:%s:%s:%s:%s#%s/%s' % (
        record_data[0],
        record_data[1],
        record_data[2],
        record_data[3],
        record_data[4],
        record_data[5],
        record_data[6],
        record_data[7])

    return '@%s\n%s\n+\n%s' % (header, seq, qual), pass_filter


def format_mapping_html_data(header,
                             mapping_data,
                             errors,
                             warnings):
    """ Generates list of lines of html for displaying problems with metadata

    This function assumes that errors and warnings each are tab separated with
    an index in the form of row,column where 0,0 is the SampleID position.  If
    the position -1,-1 is given, this indicates that the error is something
    missing entirely from the mapping file (eg., an added demultiplex field
    that does not exist).

    header:  list of header strings
    mapping_data:  list of lists of raw metadata mapping file data
    errors:  list of errors, should contain indices
    warnings:  list of warnings, should contain indices
    """

    html_lines = HTML_LINES_INIT

    if not errors and not warnings:
        html_lines += "<h1>No errors or warnings detected.<br></h1>"
    # Find errors/warnings that are not in particular cells
    general_errors = ""
    for curr_err in errors:
        loc = curr_err.split('\t')[1].strip()
        if loc == "-1,-1":
            general_errors += '<td bgcolor="red"><font color="white">' +\
                curr_err.split('\t')[0] + '<font color="black"></td>'
    general_warnings = ""
    for curr_warning in warnings:
        loc = curr_warning.split('\t')[1].strip()
        if loc == "-1,-1":
            general_warnings += '<td bgcolor="yellow">' +\
                curr_err.split('\t')[0] + "</td>"

    html_lines += HTML_LINES_MSG % ("+-%./ :,;_",
                                    general_errors, general_warnings)

    # Check header fields, color and add popup messages if errors/warnings
    # are present
    formatted_header = ""
    for curr_field in range(len(header)):
        all_errs_warnings = ""
        curr_pos = "%s,%s" % (0, curr_field)
        for curr_warning in warnings:
            loc = curr_warning.split('\t')[1].strip()
            if loc == curr_pos:
                bg_color = "yellow"
                font_color = "black"
                all_errs_warnings += curr_warning.split('\t')[0] + "<br>"

        for curr_err in errors:
            loc = curr_err.split('\t')[1].strip()
            if loc == curr_pos:
                bg_color = "red"
                font_color = "white"
                all_errs_warnings += curr_err.split('\t')[0] + "<br>"
        if not all_errs_warnings:
            formatted_header += "<th>%s</th>" % header[curr_field]
        elif not header[curr_field]:
            formatted_header += """<th bgcolor=%s><a href="javascript:void(0);" onmouseover="return overlib('%s');" onmouseout="return nd();"><font color=%s>%s</a></th>""" % (
                bg_color,
                all_errs_warnings.replace(
                    '"',
                    ''),
                font_color,
                "missing data")

        else:
            formatted_header += """<th bgcolor=%s><a href="javascript:void(0);" onmouseover="return overlib('%s');" onmouseout="return nd();"><font color=%s>%s</a></th>""" % (
                bg_color,
                all_errs_warnings.replace('"',
                                          ''),
                font_color,
                header[curr_field])

    html_lines += HTML_LINES_HEADER % formatted_header

    formatted_data = ""
    correction_ix = 1

    for curr_row in range(len(mapping_data)):
        formatted_data += "<tr>"
        for curr_cell in range(len(mapping_data[curr_row])):
            all_errs_warnings = ""
            append_location = False
            curr_pos = "%s,%s" % (curr_row + correction_ix, curr_cell)
            for curr_warning in warnings:
                loc = curr_warning.split('\t')[1].strip()
                if loc == curr_pos:
                    append_location = True
                    bg_color = "yellow"
                    font_color = "black"
                    all_errs_warnings += curr_warning.split('\t')[0] + "<br>"

            for curr_err in errors:
                loc = curr_err.split('\t')[1].strip()
                if loc == curr_pos:
                    append_location = True
                    bg_color = "red"
                    font_color = "white"
                    all_errs_warnings += curr_err.split('\t')[0] + "<br>"
            if append_location:
                if len(mapping_data[curr_row][0]) == 0:
                    sample_id_name = "missing sample id"
                else:
                    sample_id_name = mapping_data[curr_row][0]
                try:
                    header_location = header[curr_cell]
                except IndexError:
                    header_location = "no header"
                location_desc = "Location (SampleID,Header Field)<br>%s,%s" %\
                    (sample_id_name, header_location)
            if not all_errs_warnings:
                formatted_data += "<th><tt>%s</tt></th>" %\
                    mapping_data[curr_row][curr_cell]
            elif len(mapping_data[curr_row][curr_cell].replace('\n', '')) == 0:
                formatted_data += """<th bgcolor=%s><a href="javascript:void(0);" onmouseover="return overlib('%s');" onmouseout="return nd();"><font color=%s><tt>%s</tt></a></th>""" % (
                    bg_color,
                    all_errs_warnings.replace(
                        '"',
                        '').replace(
                        "'",
                        "") + location_desc,
                    font_color,
                    "missing data")
            else:
                formatted_data += """<th bgcolor=%s><a href="javascript:void(0);" onmouseover="return overlib('%s');" onmouseout="return nd();"><font color=%s><tt>%s</tt></a></th>""" % (
                    bg_color,
                    all_errs_warnings.replace('"',
                                              '').replace("'",
                                                          "") + location_desc,
                    font_color,
                    mapping_data[curr_row][curr_cell])

        formatted_data += "</tr>"
    html_lines += HTML_LINES_DATA % formatted_data

    return html_lines


def format_te_prefs(prefs_dict):
    """ Format a preferences file for use with TopiaryExplorer tep file """

    # get the sample colors
    sample_coloring = prefs_dict['sample_coloring']
    lines = []

    # for each sample get the color as hsv and write preference lines
    for k in sample_coloring:
        for t in sample_coloring[k]['colors']:
            if(isinstance(t, tuple)):
                lines.append(''.join([str(i) + ',' for i in t[1]]) + '\n')
            if(isinstance(t, str)):
                lines.append(t + ':' +
                             ''.join([str(i) + ',' for i in data_color_hsv[sample_coloring[k]['colors'][t]]]) + '\n')
        lines.append('>default' + k + ':' + k + '\n')

    return lines


def format_tep_file_lines(otu_table_data, mapping_lines, tree_lines,
                          prefs_dict):
    """ Format the tep file for TopiaryExplorer """

    # write tree file lines
    lines = ['>>tre\n']
    lines += [tree_lines.read()]
    lines += '\n'

    # get otu table data
    if(otu_table_data.metadata(axis='observation')):
        lines += ['>>otm\n#OTU ID\tOTU Metadata\n']
        for i in range(len(otu_table_data.ids(axis='observation'))):
            new_string = otu_table_data.ids(axis='observation')[i] + '\t'
            for m in otu_table_data.metadata(axis='observation')[i]['taxonomy']:
                new_string += m + ';'
            lines += [new_string]
            lines += '\n'

    # format and write otu table and taxonomy lines
    lines += ['>>osm\n']
    if otu_table_data.metadata(axis='observation') is None:
        lines += [str(otu_table_data.delimited_self())]
    elif "taxonomy" in otu_table_data.metadata(axis='observation')[0]:
        lines += [str(otu_table_data.delimited_self(header_key="taxonomy",
                                                   header_value="Consensus Lineage",
                                                   metadata_formatter=lambda x: ';'.join(x)))]

    # write mapping file lines
    lines += ['\n>>sam\n']
    lines += mapping_lines.readlines()

    # if prefs file supplied, write pref lines
    if prefs_dict:
        te_prefs = format_te_prefs(prefs_dict)
        lines += ['\n>>pre\n']
        lines += te_prefs

    return lines


def format_jnlp_file_lines(web_flag, url, tep_fp):
    """ Format the jnlp file for TopiaryExplorer """

    # write the jnlp header
    lines = [jnlp_top_block]

    # write the location of TopiaryExplorer location
    if(web_flag):
        lines += ['http://topiaryexplorer.sourceforge.net/app/']
    else:
        topiaryexplorer_project_dir =\
            load_qiime_config()['topiaryexplorer_project_dir']
        if topiaryexplorer_project_dir:
            lines += ['file:' + topiaryexplorer_project_dir]
        else:
            print "WARNING: Couldn't create jnlp file - topiaryexplorer_project_dir is not defined in your qiime_config. tep file was created sucessfully."

    # write the jnlp body text
    lines += [jnlp_middle_block]
    if(url):
        lines += [url]
    else:
        lines += [abspath(tep_fp)]

    # write the jnlp footer
    lines += [jnlp_bottom_block]

    return lines


def format_fastq_record(label,
                        seq,
                        qual):
    """ Formats a line of fastq data to be written

    label: fastq label
    seq: nucleotide sequence
    qual: quality scores
    """

    return "@%s\n%s\n+\n%s\n" % (label, seq, qual)


HTML_LINES_INIT = """<html>
<head>

<script type="text/javascript" src="./overlib.js"></script>
</head>
<body bgcolor="white"> """

HTML_LINES_MSG = """<h1>Mapping file error and warning details.</h1>
Notes for interpreting this report:
<ul>
    <li>Errors will be listed in red, warnings in yellow.
    <li>Mouse over an error or warning in a cell for more details.
    <li>Errors in the header row may mask other errors, so these should be corrected first.
    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>validate_mapping_file.py</tt> until no errors (nor warnings, ideally) are found.
</ul>
<p>
Some general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):
<ul>
    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.
    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>
    <li>Other fields allow alphanumeric and <tt>%s</tt> characters.
</ul>
General issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr>%s%s</tr></table><br>"""


HTML_LINES_HEADER = """
<table border="2" cellspacing="0" cellpadding="5">

<tr></tr>
<tr>
%s
</tr>
"""

HTML_LINES_DATA = """
<tr>
%s
</tr>
</table>

</body>
</html>"""


jnlp_top_block = """
<?xml version="1.0" encoding="utf-8"?>

<jnlp codebase=\""""

jnlp_middle_block = """\">

    <information>
    <title>TopiaryExplorer</title>
    <vendor>University of Colorado</vendor>
    <description>TopiaryExplorer</description>

    <offline-allowed/>

    </information>

    <security>
        <all-permissions/>
    </security>

    <resources>
        <j2se version="1.6+" initial-heap-size="500M" max-heap-size="2000m" />

        <jar href="topiaryexplorer1.0.jar" />
        <jar href="lib/core.jar" />
        <jar href="lib/itext.jar" />
        <jar href="lib/pdf.jar" />
        <jar href="lib/ojdbc14.jar" />
        <jar href="lib/opengl.jar" />
        <jar href="lib/mysql-connector-java-5.1.10-bin.jar" />
        <jar href="lib/javaws.jar" />
        <jar href="lib/classes12.jar" />
        <jar href="lib/jogl.jar" />
        <jar href="lib/guava-r09.jar" />
    </resources>

    <application-desc main-class="topiaryexplorer.TopiaryExplorer">
    <argument>"""

jnlp_bottom_block = """</argument>
    </application-desc>
</jnlp>
"""
