#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Contains functions used in the compare_taxa_summaries.py script."""

from numpy import array, sqrt
from qiime.stats import correlation_t
from qiime.format import (format_correlation_info, format_correlation_vector,
                          format_taxa_summary)

# Define valid choices for various options here so that they can be used by the
# library code and the script.
comparison_modes = ['paired', 'expected']
correlation_types = ['pearson', 'spearman']
tail_types = ['low', 'high', 'two-sided']


def compare_taxa_summaries(taxa_summary1, taxa_summary2, comparison_mode,
                           correlation_type='pearson', tail_type='two-sided',
                           num_permutations=999, confidence_level=0.95,
                           perform_detailed_comparisons=False,
                           sample_id_map=None, expected_sample_id=None):
    """Compares two taxa summaries using the specified comparison mode.

    Taxa summaries are compared by computing the correlation coefficient
    between samples in the two taxa summaries based on the abundance of various
    taxa. The two taxa summaries are sorted and filled such that taxa that are
    missing in one summary but are present in the other are represented with
    zero abundance.

    Returns a four-element tuple containing the following values: the sorted
    and filled taxa summaries (in a format ready to be written to a file) and a
    report detailing the correlation between the taxa summaries (also in a
    format ready to be written to a file), including the correlation
    coefficient, parametric and nonparametric p-values, and confidence interval
    for the overall comparison.

    If perform_detailed_comparisons is True, a correlation vector is returned
    (also in a format ready to be written to a file) where each line shows the
    samples that were compared and the associated correlation coefficient,
    parametric and nonparametric p-values (uncorrected and
    Bonferroni-corrected), and the confidence interval. If
    perform_detailed_comparisons is False, None will be returned for this
    value.

    Arguments:
        taxa_summary1 - the first taxa summary to compare. This should be a
            tuple containing the sample IDs, taxa, and taxonomic data (i.e. the
            output of qiime.parse.parse_taxa_summary_table)
        taxa_summary2 - the second taxa summary to compare
        comparison_mode - the type of comparison to perform on the two taxa
            summaries. Can be either 'paired' or 'expected'. If 'paired', the
            samples that match between the two taxa summaries will be compared,
            unless a sample_id_map is specified. If 'expected', each sample in
            the first taxa summary will be compared to an 'expected' sample in
            the second taxa summary. If 'expected', the second taxa summary
            must only contain a single sample unless expected_sample_id is
            provided
        correlation_type - the type of correlation coefficient to calculate
            when comparing samples in the taxa summaries. Can be either
            'pearson' or 'spearman'
        tail_type - if 'two-sided', a two-sided test is performed. 'high'
            for a one-tailed test for positive association, or 'low' for a
            one-tailed test for negative association. This parameter affects
            both the parametric and nonparametric tests, but the confidence
            interval will always be two-sided
        num_permutations - the number of permutations to use in the
            nonparametric test. Must be a number greater than or equal to 0. If
            0, the nonparametric test will not be performed. In this case, the
            nonparametric p-value will be 'N/A' in the formatted results string
        confidence_level - the confidence level to use when constructing the
            confidence interval. Must be between 0 and 1 (exclusive)
        perform_detailed_comparisons - if True, computes the correlation
            between pairs of samples in addition to computing the overall
            correlation between taxa summaries
        sample_id_map - a dictionary mapping original sample IDs to new sample
            IDs. New sample IDs that match will be compared. All original
            sample IDs must be mapped. This argument is only used if the
            comparison mode is 'paired'. If not provided, only matching sample
            IDs between the two taxa summaries will be compared
        expected_sample_id - the sample ID in taxa_summary2 to compare all
            samples in taxa_summary1 to. This argument is only used if the
            comparison mode is 'expected'. If not provided, taxa_summary2 must
            only contain a single sample, and all samples in taxa_summary1 will
            be compared to it
    """
    # Perform some initial error checking before getting into the heavy
    # processing.
    if correlation_type not in correlation_types:
        raise ValueError("Invalid correlation type '%s'. Must be one of %r." %
                         (correlation_type, correlation_types))
    if tail_type not in tail_types:
        raise ValueError("Invalid tail type '%s'. Must be one of %r." %
                         (tail_type, tail_types))
    if num_permutations < 0:
        raise ValueError("Invalid number of permutations: %d. Must be greater "
                         "than or equal to zero." % num_permutations)
    if confidence_level <= 0 or confidence_level >= 1:
        raise ValueError("Invalid confidence level: %.4f. Must be between "
                         "zero and one (exclusive)." % confidence_level)

    # Define some comments to be put in the result strings.
    header = "# Correlation coefficient: %s.\n" % correlation_type
    header += "# The parametric p-value(s) were calculated using a "

    if tail_type == 'two-sided':
        tail_type_desc = tail_type
    elif tail_type == 'high':
        tail_type_desc = "one-sided (positive association)"
    elif tail_type == 'low':
        tail_type_desc = "one-sided (negative association)"

    header += tail_type_desc + " test of significance using a " + \
        "t-distribution.\n"

    if num_permutations > 0:
        header += "# The nonparametric p-value(s) were calculated using " + \
                  "a " + tail_type_desc + " permutation test with " + \
                  str(num_permutations) + " permutations.\n"

    header += "# The confidence interval(s) were constructed at a " + \
              "confidence level of " + str(confidence_level * 100) + \
              "% using Fisher's z-transformation (see Sokal and Rohlf " + \
              "3rd edition pg. 575). The confidence interval(s) are two-sided."

    spearman_overall_warning = "# Since there were 10 or fewer " + \
        "observations when calculating Spearman's rank correlation " + \
        "coefficient, the parametric p-value is "
    spearman_detailed_warning = "# Since there were 10 or fewer taxa in " + \
        "the sorted and filled taxa summary files, the parametric " + \
        "p-values and Bonferroni-corrected parametric p-values are "
    spearman_warning_suffix = "not accurate when using the " + \
        "t-distribution. Please see Biometry (Sokal and Rohlf, " + \
        "3rd edition) page 600 for more details."
    spearman_overall_warning += spearman_warning_suffix
    spearman_detailed_warning += spearman_warning_suffix

    # Sort and fill the taxa summaries so that we can compare them.
    filled_ts1, filled_ts2 = _sort_and_fill_taxa_summaries([taxa_summary1,
                                                            taxa_summary2])
    if comparison_mode == 'paired':
        # Make sure that each sample is paired up to the sample it needs to be
        # compared against according to the sample ID map.
        compatible_ts1, compatible_ts2 = _make_compatible_taxa_summaries(
            filled_ts1, filled_ts2, sample_id_map)
        overall_corr, corr_vec = _compute_correlation(compatible_ts1,
                                                      compatible_ts2, comparison_mode, correlation_type, tail_type,
                                                      num_permutations, confidence_level,
                                                      perform_detailed_comparisons)

        # Calculate the length of the vectors that were used to compute
        # correlation of.
        num_overall_observations = len(compatible_ts1[0]) * \
            len(compatible_ts1[1])

        # Report the number of samples that matched.
        header += "\n# Number of samples that matched between the taxa " + \
                  "summary files: %d" % len(compatible_ts1[0])
    elif comparison_mode == 'expected':
        overall_corr, corr_vec = _compute_correlation(filled_ts1, filled_ts2,
                                                      comparison_mode, correlation_type, tail_type, num_permutations,
                                                      confidence_level, perform_detailed_comparisons,
                                                      expected_sample_id)
        num_overall_observations = len(filled_ts1[0]) * len(filled_ts1[1])
    else:
        raise ValueError("Invalid comparison mode '%s'. Must be one of %r." %
                         (comparison_mode, comparison_modes))

    # Format the overall correlation into a string that is writable to a file.
    # Include a warning in the header if the correlation coefficient was
    # spearman and the number of observations was <= 10.
    overall_corr_str_header = header
    if correlation_type == 'spearman' and num_overall_observations <= 10:
        overall_corr_str_header += '\n' + spearman_overall_warning
    overall_corr_str = format_correlation_info(overall_corr[0],
                                               overall_corr[1], overall_corr[
                                                   2], overall_corr[3],
                                               num_permutations, overall_corr_str_header)

    # Format the correlation vector.
    corr_vec_str = None
    if perform_detailed_comparisons:
        detailed_header = header
        if correlation_type == 'spearman' and len(filled_ts1[1]) <= 10:
            detailed_header += '\n' + spearman_detailed_warning
        corr_vec_str = format_correlation_vector(corr_vec, num_permutations,
                                                 detailed_header)

    return (format_taxa_summary(filled_ts1), format_taxa_summary(filled_ts2),
            overall_corr_str, corr_vec_str)


def _make_compatible_taxa_summaries(ts1, ts2, sample_id_map=None):
    """Returns two taxa summaries that are ready for direct comparison.

    The returned taxa summaries will have their samples ordered such that
    direct comparisons may be made between each of the corresponding samples in
    the two summaries. For example, if ts1 has samples 'S1' and 'S2' and ts2
    has samples 'S2' and 'S3' (and no sample ID map is provided), the resulting
    taxa summaries will only have the sample 'S2'.

    As another example, assume ts1 has samples 'S1' and 'S2' and ts2 has
    samples 'T1' and 'T2'. A sample ID map may be provided that maps 'S1' to
    'T1' and 'S2' to 'T2'. The first resulting compatible taxa summary
    will have the samples 'S1' and 'S2' and the second resulting compatible
    taxa summary will have samples 'T1' and 'T2'. Thus, these resulting taxa
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
    # Check to make sure the sample ID map looks sane. All samples IDs in both
    # taxa summaries should be mapped, though extra mappings (i.e. if the taxa
    # summaries have been filtered) are allowed.
    if sample_id_map:
        for samp_id in ts1[0] + ts2[0]:
            if samp_id not in sample_id_map:
                raise ValueError("The original sample ID '%s' does not have a "
                                 "mapping in the sample ID map. All sample IDs must "
                                 "have a mapping." % samp_id)

    # For each sample ID in the first taxa summary file, try to find a matching
    # sample ID (using the sample ID map if one was provided) in the second
    # file. Create two new taxa summary files containing the matching samples
    # in the same order.
    new_samp_ids1, new_samp_ids2, new_data1, new_data2 = [], [], [], []
    for samp_idx, samp_id in enumerate(ts1[0]):
        matching_samp_id = None
        if sample_id_map:
            # Find the matching sample ID that the current sample ID maps to.
            # Only one-to-one mappings are allowed when the sample ID map is
            # first parsed, though 'broken' mappings are allowed (i.e. a sample
            # maps to a new sample ID, but no other sample maps to that ID). In
            # this case, we simply ignore the broken mapping.
            new_samp_id = sample_id_map[samp_id]
            for orig_samp_id in sample_id_map:
                if (orig_samp_id != samp_id and
                        sample_id_map[orig_samp_id] == new_samp_id):
                    matching_samp_id = orig_samp_id
        else:
            if samp_id in ts2[0]:
                matching_samp_id = samp_id
        if matching_samp_id:
            new_samp_ids1.append(samp_id)
            new_samp_ids2.append(matching_samp_id)

            # Transpose the data matrix so that we can index by sample index.
            new_data1.append(ts1[2].T[samp_idx])
            new_data2.append(ts2[2].T[ts2[0].index(matching_samp_id)])
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
    master_taxa = sorted(set(master_taxa))

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


def _compute_correlation(ts1, ts2, comparison_mode, correlation_type,
                         tail_type, num_permutations, confidence_level,
                         perform_detailed_comparisons=False,
                         expected_sample_id=None):
    """Computes the correlation between two taxa summary files.

    The input taxa summaries MUST already be sorted and filled (see
    _sort_and_fill_taxa_summaries) so that the various taxa line up and contain
    the same number of taxa (e.g. the first taxon in both files is 'Bacteria'
    and not mismatched).

    Additionally, if comparison_mode is 'paired', the input taxa summaries must
    already have been made compatible, meaning the number of samples match
    between the two taxa summaries. This is very important as the first sample
    in ts1 will be compared to the first sample in ts2, the second sample in
    ts1 will be compared to the second sample in ts2, and so on. The sample IDs
    are not checked by this function to ensure they are in the correct order or
    mapping (this is the job of _make_compatible_taxa_summaries).

    Returns a two-element tuple: the first element is a four-element tuple
    containing the correlation coefficient, parametric p-value, nonparametric
    p-value, and a tuple for the confidence interval for the overall
    comparison.

    If perform_detailed_comparisons is True, the second element is a
    correlation vector, which is a list of 8-element tuples, where the first
    element is a sample ID from ts1, the second element is a sample ID from
    ts2, the third element is the correlation coefficient computed between the
    two samples (a double), the fourth element is the parametric p-value, the
    fifth element is the Bonferroni-corrected parametric p-value, the sixth
    element is the nonparametric p-value, the seventh element is the
    Bonferroni-corrected nonparametric p-value, and the eighth element is a
    tuple containing the low and high ends of the confidence interval. If
    perform_detailed_comparisons is False, None will be returned for the second
    element.

    Arguments:
        ts1 - the first taxa summary to be compared
        ts2 - the second taxa summary to be compared
        comparison_mode - the type of comparison to perform on the two taxa
            summaries. Can be either 'paired' or 'expected'. If 'paired', each
            positional pair of samples between the two taxa summaries will be
            compared. If 'expected', each sample in ts1 will be compared to an
            'expected' sample in ts2. If 'expected', ts2 must only contain a
            single sample unless expected_sample_id is provided
        correlation_type - the type of correlation coefficient to calculate
            when comparing samples in the taxa summaries. Can be either
            'pearson' or 'spearman'
        tail_type - if 'two-sided', a two-sided test is performed. 'high'
            for a one-tailed test for positive association, or 'low' for a
            one-tailed test for negative association. This parameter affects
            both the parametric and nonparametric tests, but the confidence
            interval will always be two-sided
        num_permutations - the number of permutations to use in the
            nonparametric test. Must be a number greater than or equal to 0. If
            0, the nonparametric test will not be performed. In this case, the
            nonparametric p-values will be None
        confidence_level - the confidence level to use when constructing the
            confidence interval. Must be between 0 and 1 (exclusive)
        perform_detailed_comparisons - if True, computes the correlation
            between each pair of samples in addition to computing the overall
            correlation between taxa summaries
        expected_sample_id - the sample ID in ts2 to compare all samples in ts1
            to. If not provided, ts2 must only contain a single sample, and all
            samples in ts1 will be compared to it
    """
    if comparison_mode != 'paired' and comparison_mode != 'expected':
        raise ValueError("Invalid comparison mode '%s'. Must be one of %r." %
                         (comparison_mode, comparison_modes))

    # Make sure that the second taxa summary has only one sample if we weren't
    # provided an expected sample ID to compare against.
    if (comparison_mode == 'expected' and expected_sample_id is None and
            len(ts2[0]) != 1):
        raise ValueError("The second taxa summary file must contain a single "
                         "sample (column) to compare all samples in the first taxa "
                         "summary file against when the comparison mode is 'expected' "
                         "and an expected sample ID is not provided. You provided a "
                         "file with %d samples."
                         % len(ts2[0]))

    if comparison_mode == 'paired':
        # Make sure the number of samples match between the two files (the IDs
        # do not have to match because of the sample ID map).
        if len(ts1[0]) != len(ts2[0]):
            raise ValueError("The two taxa summaries are incompatible because "
                             "they do not have the same number of sample IDs. "
                             "The taxa summaries must be made compatible "
                             "before attempting to perform "
                             "pairwise-comparisons between samples.")

    # Make sure the taxa information is the same (i.e. the summaries have been
    # sorted and filled).
    if ts1[1] != ts2[1]:
        raise ValueError("The taxa do not match exactly between the two taxa "
                         "summary files. The taxa must be sorted and filled "
                         "before attempting to compare them.")

    # Find the index of the expected sample ID.
    if comparison_mode == 'expected':
        if expected_sample_id:
            try:
                expected_idx = ts2[0].index(expected_sample_id)
            except ValueError:
                raise ValueError("The expected sample ID '%s' is not in the "
                                 "taxa summary file." % expected_sample_id)
        else:
            # We know the 'expected' taxa summary has a single sample in it, so
            # this is the only possible index.
            expected_idx = 0

    # Compute the overall correlation between each sample and the expected
    # sample, or each of the paired samples, and optionally the correlation
    # between each pair of samples individually.
    corr_vec = None
    if perform_detailed_comparisons:
        corr_vec = []
        num_comparisons = len(ts1[0])

    all_ts1_data = []
    all_ts2_data = []
    for samp_idx, samp_id in enumerate(ts1[0]):
        if comparison_mode == 'paired':
            paired_idx = samp_idx
        elif comparison_mode == 'expected':
            paired_idx = expected_idx
        else:
            # Redundant check, but here for safety in case the one above is
            # changed or removed.
            raise ValueError("Invalid comparison mode '%s'. Must be one of "
                             "%r." % (comparison_mode, comparison_modes))

        # Grab the columns of data for the current sample and its pair.
        ts1_data = ts1[2].T[samp_idx]
        ts2_data = ts2[2].T[paired_idx]
        all_ts1_data.extend(ts1_data)
        all_ts2_data.extend(ts2_data)

        if perform_detailed_comparisons:
            # Compare the current sample and its pair.
            corr_coeff, param_p_val, unused, nonparam_p_val, conf_interval = \
                correlation_t(ts1_data, ts2_data,
                                 method=correlation_type,
                                 tails=tail_type,
                                 permutations=num_permutations,
                                 confidence_level=confidence_level)

            # Compute the Bonferroni-corrected p-values.
            param_p_val_corr = min(param_p_val * num_comparisons, 1)
            nonparam_p_val_corr = None if nonparam_p_val is None else \
                min(nonparam_p_val * num_comparisons, 1)

            corr_vec.append((samp_id, ts2[0][paired_idx], corr_coeff,
                             param_p_val, param_p_val_corr, nonparam_p_val,
                             nonparam_p_val_corr, conf_interval))

    # Compare all paired samples at once.
    results = correlation_t(all_ts1_data, all_ts2_data,
                               method=correlation_type, tails=tail_type,
                               permutations=num_permutations,
                               confidence_level=confidence_level)
    # We don't need to return all of the permuted correlation coefficients.
    overall_corr = (results[0], results[1], results[3], results[4])
    return overall_corr, corr_vec
