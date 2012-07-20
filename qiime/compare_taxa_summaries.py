#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Contains functions used in the compare_taxa_summaries.py script."""

from os.path import basename, splitext
from numpy import array, sqrt
from cogent.maths.stats.test import pearson

# Define valid choices for comparison mode and correlation type here so that it
# can be used by the library code and the script.
comparison_modes = ['paired', 'expected']
correlation_types = ['pearson', 'spearman']

def compare_taxa_summaries(taxa_summary1, taxa_summary2, comparison_mode,
                           correlation_type='pearson', sample_id_map=None):
    """Compares two taxa summaries using the specified comparison mode.

    Taxa summaries are compared by computing the correlation coefficient
    between samples in the two taxa summaries based on the abundance of various
    taxa. The two taxa summaries are sorted and filled such that taxa that are
    missing in one summary but are present in the other are represented with
    zero abundance.

    Returns the sorted and filled taxa summaries (in a format ready to be
    written to a file) and a correlation vector (also in a format ready to be
    written to a file) where each line shows the samples that were compared and
    the associated correlation coefficient.

    Arguments:
        taxa_summary1 - the first taxa summary to compare. This should be a
            tuple containing the sample IDs, taxa, and taxonomic data (i.e. the
            output of qiime.parse.parse_taxa_summary_table)
        taxa_summary2 - the second taxa summary to compare
        comparison_mode - the type of comparison to perform on the two taxa
            summaries. Can be either 'paired' or 'expected'. If 'paired', the
            samples that match between the two taxa summaries will be compared,
            unless a sample_id_map is specified. If 'expected', each sample in
            the first taxa summary will be compared to the single 'expected'
            sample in the second taxa summary. If 'expected', the second taxa
            summary must only contain a single sample
        correlation_type - the type of correlation coefficient to calculate
            when comparing samples in the taxa summaries. Can be either
            'pearson' or 'spearman'
        sample_id_map - a dictionary mapping sample IDs in the first taxa
            summary to be compared to sample IDs in the second taxa summary.
            Multiple keys (i.e. sample IDs in the first taxa summary) may map
            to the same value (i.e. sample IDs in the second taxa summary).
            Comparisons will only be made between mappings that exist in this
            dictionary. This argument is only used if the comparison mode is
            'paired'. If not provided, only matching sample IDs between the two
            taxa summaries will be compared
    """
    # Sort and fill the taxa summaries so that we can compare them.
    correlation_fn = _get_correlation_function(correlation_type)
    filled_ts1, filled_ts2 = _sort_and_fill_taxa_summaries([taxa_summary1,
                                                            taxa_summary2])
    if comparison_mode == 'paired':
        # Make sure that each sample is paired up to the sample it needs to be
        # compared against according to the sample ID map.
        compatible_ts1, compatible_ts2 = _make_compatible_taxa_summaries(
                filled_ts1, filled_ts2, sample_id_map)
        correlations = _compute_paired_sample_correlations(compatible_ts1,
                compatible_ts2, correlation_fn)
    elif comparison_mode == 'expected':
        correlations = _compute_all_to_expected_correlations(filled_ts1,
                filled_ts2, correlation_fn)
    else:
        raise ValueError("Invalid comparison mode '%s'. Must be one of %r." %
                         (comparison_mode, comparison_modes))

    # Format the correlation vector with a header describing the correlation
    # coefficient that was used.
    header = "# Correlation coefficient: %s" % correlation_type
    correlation_vector = _format_correlation_vector(correlations, header)
    return (_format_taxa_summary(filled_ts1), _format_taxa_summary(filled_ts2),
           correlation_vector)

def parse_sample_id_map(sample_id_map_f):
    """Parses the lines of a sample ID map file into a dictionary.

    Returns a dictionary with sample IDs as the keys and values.

    Arguments:
        sample_id_map_f - the lines of a sample ID map file to parse. Each line
            should contain two sample IDs separated by a tab. Each value in the
            first column must be unique, since the returned data structure is a
            dictionary using those values as keys
    """
    result = {}
    for line in sample_id_map_f:
        line = line.strip()
        if line:
            # Only try to parse lines that aren't just whitespace.
            samp_id, mapped_id = line.split('\t')
            if samp_id in result:
                raise ValueError("The first column of the sample ID map must "
                                 "contain unique sample IDs ('%s' is "
                                 "repeated). The second column, however, may "
                                 "contain repeats." % samp_id)
            else:
                result[samp_id] = mapped_id
    return result

def add_filename_suffix(filepath, suffix):
    """Adds a suffix to the filepath, inserted before the file extension.

    Returns the new filepath string. For example, if filepath is 'foo.txt' and
    suffix is '_bar', 'foo_bar.txt' will be returned.

    Arguments:
        filepath - any filepath to append the suffix to (before the file
            extension, if it exists). Most useful if the filepath points to a
            file instead of a directory
    """
    root, extension = splitext(basename(filepath))
    return root + suffix + extension

def _format_correlation_vector(correlations, header=''):
    """Formats a correlation vector to be suitable for writing to a file.

    Returns a string where each line contains three tab-separated fields: the
    two sample IDs that were compared, and the computed correlation
    coefficient.

    Arguments:
        correlations - a list of 3-element tuples, where the first element is a
            sample ID, the second element is a sample ID, and the third element
            is the correlation coefficient computed between the two samples (a
            double)
        header - if provided, this string will be inserted at the beginning of
            the returned string. For example, might be useful to add a comment
            describing what correlation coefficient was used. This string does
            not need to contain a newline at the end
    """
    result = ''
    if header != '':
        result += header + '\n'
    for samp_id1, samp_id2, correlation in correlations:
        result += '%s\t%s\t%.4f\n' % (samp_id1, samp_id2, correlation)
    return result

def _format_taxa_summary(taxa_summary):
    """Formats a taxa summary to be suitable for writing to a file.

    Returns a string where the first line contains 'Taxon' and then a list of
    tab-separated sample IDs. The next lines contain each taxon and its
    abundance separated by a tab.

    Arguments:
        taxa_summary - the taxa summary tuple (i.e. the output of
            qiime.parse.parse_taxa_summary_table) to format
    """
    result = 'Taxon\t' + '\t'.join(taxa_summary[0]) + '\n'
    for taxon, row in zip(taxa_summary[1], taxa_summary[2]):
        row = map(str, row)
        result += '%s\t' % taxon + '\t'.join(row) + '\n'
    return result

def _get_correlation_function(correlation_type):
    """Returns the correlation function to use based on the correlation type.

    The correlation function that is returned will always accept two lists of
    numbers to compute the correlation coefficient on.

    Arguments:
        correlation_type - a string specifying the correlation type (may be
            either 'pearson' or 'spearman')
    """
    if correlation_type == 'pearson':
        correlation_fn = _pearson_correlation
    elif correlation_type == 'spearman':
        correlation_fn = _spearman_correlation
    else:
        raise ValueError("Invalid correlation type '%s'. Must be one of %r." %
                         (correlation_type, correlation_types))
    return correlation_fn

def _make_compatible_taxa_summaries(ts1, ts2, sample_id_map=None):
    """Returns two taxa summaries that are ready for direct comparison.

    The returned taxa summaries will have their samples ordered such that
    direct comparisons may be made between each of the corresponding samples in
    the two summaries. For example, if ts1 has samples 'S1' and 'S2' and ts2
    has samples 'S2' and 'S3' (and no sample ID map is provided), the resulting
    taxa summaries will only have the sample 'S2'.

    As another example, assume ts1 has samples 'S1' and 'S2' and ts2 has
    samples 'T1' and 'T2'. A sample ID map may be provided that maps 'S1' to
    'T1' and 'S2' also to 'T1'. The first resulting compatible taxa summary
    will have the samples 'S1' and 'S2' and the second resulting compatible
    taxa summary will have samples 'T1' and 'T1' (repeated due to the
    many-to-one mapping in the sample ID map). Thus, these resulting taxa
    summaries can be directly compared because each sample lines up to the one
    it needs to be compared to.

    The input taxa summaries do not need to be already sorted and filled, but
    it is okay if they are. The taxonomic information is not altered by this
    function, only the order and presence/absence of samples (i.e. columns in
    the table). It is a good idea to sort and filter the taxa summaries (either
    before or after) before computing the correlation coefficients between
    samples.

    Arguments:
        ts1 - the first taxa summary
        ts2 - the second taxa summary
        sample_id_map - a dictionary describing which samples in the first taxa
            summary should be compared to which samples in the second taxa
            summary. If not provided, only samples whose sample IDs directly
            match will be compared
    """
    # Check to make sure the sample ID map looks sane.
    if sample_id_map:
        bad_samp_ids = set(sample_id_map.keys()).difference(set(ts1[0]))
        if len(bad_samp_ids) != 0:
            raise ValueError("The sample IDs %r in the sample ID map do not "
                             "match any of the sample IDs in the taxa summary "
                             "file." % bad_samp_ids)

    # For each sample ID in the first taxa summary file, try to find a matching
    # sample ID (using the sample ID map if one was provided) in the second
    # file. Create two new taxa summary files containing the matching samples
    # in the same order.
    new_samp_ids1, new_samp_ids2, new_data1, new_data2 = [], [], [], []
    for samp_idx, samp_id in enumerate(ts1[0]):
        matching_samp_id = None
        if sample_id_map:
            if samp_id in sample_id_map:
                matching_samp_id = sample_id_map[samp_id]
        else:
            if samp_id in ts2[0]:
                matching_samp_id = samp_id
        if matching_samp_id:
            new_samp_ids1.append(samp_id)
            new_samp_ids2.append(matching_samp_id)

            # Transpose the data matrix so that we can index by sample index.
            new_data1.append(ts1[2].T[samp_idx])
            try:
                new_data2.append(ts2[2].T[ts2[0].index(matching_samp_id)])
            except ValueError:
                raise ValueError("The sample ID '%s' was not in the second "
                                 "taxa summary file. Please check your sample "
                                 "ID map." % matching_samp_id)
    if len(new_samp_ids1) == 0:
        raise ValueError("No sample IDs matched between the taxa summaries. "
                         "The taxa summaries are incompatible.")
    return (new_samp_ids1, ts1[1], array(new_data1).T), \
           (new_samp_ids2, ts2[1], array(new_data2).T)

def _sort_and_fill_taxa_summaries(taxa_summaries):
    """Sorts and fills the taxonomic information of taxa summaries.

    Returns a list of sorted and filled taxa summaries. The taxonomic
    information (i.e. rows of the taxa summary tables) are sorted and filled
    such that any missing taxa are included with abundances of zero. This makes
    it possible to compute correlation coefficients between samples because the
    taxonomic information will be the same length and in the same order. This
    is also useful for obtaining taxa summaries that can be compared in QIIME's
    plot_taxa_summary.py script.

    Arguments:
        taxa_summaries - a list of any number of taxa summaries to be sorted
            and filled
    """
    # Build up a sorted list of all taxa found in the taxa summaries
    # (no repeats).
    master_taxa = []
    for ts in taxa_summaries:
        master_taxa += ts[1]
    master_taxa = list(set(master_taxa))
    master_taxa.sort()

    # Put each taxon and data in the order it appears in the master taxa list
    # for each taxa summary. If it doesn't exist, record an abundance of zero
    # for that row.
    result = []
    for ts in taxa_summaries:
        samples = ts[0]
        orig_taxa = ts[1]
        orig_data = ts[2]
        data = []
        for taxa in master_taxa:
            try:
                taxa_index = orig_taxa.index(taxa)
                data.append(orig_data[taxa_index])
            except ValueError:
                data.append([0.] * len(samples))
        result.append((samples, master_taxa, array(data)))
    return result

def _compute_all_to_expected_correlations(observed_taxa_summary,
                                          expected_taxa_summary,
                                          correlation_fn):
    """Computes the correlation between all samples and an expected sample.

    The input taxa summaries MUST already be sorted and filled (see
    _sort_and_fill_taxa_summaries) so that the various taxa line up and contain
    the same number of taxa (e.g. the first taxon in both files is 'Bacteria'
    and not mismatched).

    Returns a correlation vector, which is a list of 3-element tuples, where
    the first element is a sample ID from observed_taxa_summary, the second
    element is the single sample ID from expected_taxa_summary, and the third
    element is the correlation coefficient computed between the two samples (a
    double).

    Arguments:
        observed_taxa_summary - the taxa summary containing samples that will
            be compared to the expected sample
        expected_taxa_summary - the taxa summary containing a single sample
            that all samples in observed_taxa_summary will be compared to
        correlation_fn - the correlation function to use when comparing
            samples. This function must accept two lists of numbers and return
            a number (the correlation coefficient)
    """
    # Make sure that the second taxa summary has only one sample.
    if len(expected_taxa_summary[0]) != 1:
        raise ValueError("The second taxa summary file must contain a single "
                "sample (column) to compare all samples in the first taxa "
                "summary file against when the comparison mode is 'expected'. "
                "You provided %d samples." % len(expected_taxa_summary[0]))

    # Make sure the taxa information is the same (i.e. the summaries have been
    # sorted and filled).
    if observed_taxa_summary[1] != expected_taxa_summary[1]:
        raise ValueError("The taxa do not match exactly between the two taxa "
                         "summary files. The taxa must be sorted and filled "
                         "before attempting to compare them.")

    # Compute the correlation between each sample and the expected sample.
    result = []
    for sample_idx, sample_id in enumerate(observed_taxa_summary[0]):
        result.append((sample_id, expected_taxa_summary[0][0],
            correlation_fn(observed_taxa_summary[2].T[sample_idx],
                           expected_taxa_summary[2].T[0])))
    return result

def _compute_paired_sample_correlations(ts1, ts2, correlation_fn):
    """Computes the correlation between each sample pair in the taxa summaries.

    The input taxa summaries MUST already be sorted and filled (see
    _sort_and_fill_taxa_summaries) so that the various taxa line up and contain
    the same number of taxa (e.g. the first taxon in both files is 'Bacteria'
    and not mismatched). Additionally, the input taxa summaries must already
    have been made compatible, meaning the number of samples match between the
    two taxa summaries. This is very important as the first sample in ts1 will
    be compared to the first sample in ts2, the second sample in ts1 will be
    compared to the second sample in ts2, and so on. The sample IDs are not
    checked by this function to ensure they are in the correct order or mapping
    (this is the job of compare_taxa_summaries).

    Returns a correlation vector, which is a list of 3-element tuples, where
    the first element is a sample ID from ts1, the second element is a sample
    ID from ts2, and the third element is the correlation coefficient computed
    between the two samples (a double).

    Arguments:
        ts1 - the first taxa summary to be compared
        ts2 - the second taxa summary to be compared
        correlation_fn - the correlation function to use when comparing
            samples. This function must accept two lists of numbers and return
            a number (the correlation coefficient)
    """
    # Make sure the number of samples match between the two files (the IDs do
    # not have to match because of the sample ID map).
    if len(ts1[0]) != len(ts2[0]):
        raise ValueError("The two taxa summaries are incompatible because "
                         "they do not have the same number of sample IDs. The "
                         "taxa summaries must be made compatible before "
                         "attempting to perform pairwise-comparisons between "
                         "samples.")

    # Make sure the taxa information is the same (i.e. the summaries have been
    # sorted and filled).
    if ts1[1] != ts2[1]:
        raise ValueError("The taxa do not match exactly between the two taxa "
                         "summary files. The taxa must be sorted and filled "
                         "before attempting to compare them.")

    # Compute the correlation between each paired sample.
    result = []
    for samp_idx, samp_id in enumerate(ts1[0]):
        corr_coeff = correlation_fn(ts1[2].T[samp_idx], ts2[2].T[samp_idx])
        result.append((samp_id, ts2[0][samp_idx], corr_coeff))
    return result

def _pearson_correlation(vec1, vec2):
    """Wraps PyCogent's pearson function to provide better error handling."""
    if len(vec1) != len(vec2):
        raise ValueError("The length of the two vectors must be the same in "
                         "order to calculate the Pearson correlation "
                         "coefficient.")
    if len(vec1) < 2:
        raise ValueError("The two vectors must both contain at least 2 "
                "elements. The vectors are of length %d." % len(vec1))
    return pearson(vec1, vec2)

# These next two functions are taken from the qiime.stats.BioEnv class, written
# by Michael Dwan. We don't use the class itself because we don't have the
# required data to instantiate one, so the methods are copied here for
# convenience. TODO: this should be moved into pycogent at some point.
def _spearman_correlation(vec1, vec2, ranked=False):
    """Calculates the the Spearman distance of two vectors."""
    try:
        temp = len(vec1)
    except ValueError:
        raise ValueError('First input vector is not a list.')
    try:
        temp = len(vec2)
    except ValueError:
        raise ValueError('Second input vector is not a list.')
    if len(vec1) == 0 or len(vec2) == 0:
        raise ValueError('One or both input vectors has/have zero elements')
    if len(vec1) != len(vec2):
        raise ValueError('Vector lengths must be equal')

    if not ranked:
        rank1, ties1 = _get_rank(vec1)
        rank2, ties2 = _get_rank(vec2)
    else:
        rank1, ties1 = vec1
        rank2, ties2 = vec2

    if ties1 == 0 and ties2 == 0:
        n = len(rank1)
        sum_sqr = sum([(x-y)**2 for x,y in zip(rank1,rank2)])
        rho = 1 - (6*sum_sqr/(n*(n**2 - 1)))
        return rho

    avg = lambda x: sum(x)/len(x)

    x_bar = avg(rank1)
    y_bar = avg(rank2)

    numerator = sum([(x-x_bar)*(y-y_bar) for x,y in zip(rank1, rank2)])
    denominator = sqrt(sum([(x-x_bar)**2 for x in rank1])*
                       sum([(y-y_bar)**2 for y in rank2]))
    # Calculate rho.
    return numerator/denominator

def _get_rank(data):
    """Ranks the elements of a list. Used in Spearman correlation."""
    indices = range(len(data))
    ranks = range(1,len(data)+1)
    indices.sort(key=lambda index:data[index])
    ranks.sort(key=lambda index:indices[index-1])
    data_len = len(data)
    i = 0
    ties = 0
    while i < data_len:
        j = i + 1
        val = data[indices[i]]
        try:
            val += 0
        except TypeError:
            raise(TypeError)

        while j < data_len and data[indices[j]] == val:
            j += 1
        dup_ranks = j - i
        val = float(ranks[indices[i]]) + (dup_ranks-1)/2.0
        for k in range(i, i+dup_ranks):
            ranks[indices[k]] = val
        i += dup_ranks
        ties += dup_ranks-1
    return ranks, ties
