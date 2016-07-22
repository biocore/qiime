#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os import path

from skbio.stats import p_value_to_str
from skbio.stats.distance import DistanceMatrix, mantel

from qiime.util import make_compatible_distance_matrices
from qiime.stats import MantelCorrelogram, PartialMantel


def run_mantel_test(method, fps, distmats, num_perms, tail_type, comment,
                    control_dm_fp=None, control_dm=None,
                    sample_id_map=None):
    """Runs a Mantel test on all pairs of distance matrices.

    Returns a string suitable for writing out to a file containing the results
    of the test.

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.

    Arguments:
        method - which Mantel test to run (either 'mantel' or 'partial_mantel')
        fps - list of filepaths of the distance matrices
        distmats - list of tuples containing dm labels and dm data (i.e. the
            output of parse_distmat)
        num_perms - the number of permutations to use to calculate the
            p-value(s)
        tail_type - the type of tail test to use when calculating the
            p-value(s). Can be 'two-sided', 'greater', or 'less'. Only applies
            when method is mantel
        comment - comment string to add to the beginning of the results string
        control_dm_fp - filepath of the control distance matrix. Only applies
            when method is partial_mantel (it is required then)
        control_dm - tuple containing control distance matrix labels and matrix
            data. Only applies when method is partial_mantel (it is required
            then)
        sample_id_map - dict mapping sample IDs (i.e. what is expected by
            make_compatible_distance_matrices)
    """
    if len(fps) != len(distmats):
        raise ValueError("Must provide the same number of filepaths as there "
                         "are distance matrices.")
    if comment is None:
        comment = ''
    result = comment

    if method == 'mantel':
        result += 'DM1\tDM2\tNumber of entries\tMantel r statistic\t' + \
                  'p-value\tNumber of permutations\tTail type\n'
    elif method == 'partial_mantel':
        if not control_dm_fp or not control_dm:
            raise ValueError("You must provide a control matrix filepath and "
                             "control matrix when running the partial Mantel "
                             "test.")
        result += 'DM1\tDM2\tCDM\tNumber of entries\t' + \
            'Mantel r statistic\tp-value\tNumber of permutations\t' +\
            'Tail type\n'
    else:
        raise ValueError("Invalid method '%s'. Must be either 'mantel' or "
                         "'partial_mantel'." % method)

    # Loop over all pairs of dms.
    for i, (fp1, (dm1_labels, dm1_data)) in enumerate(zip(fps, distmats)):
        for fp2, (dm2_labels, dm2_data) in zip(fps, distmats)[i + 1:]:
            # Make the current pair of distance matrices compatible by only
            # keeping samples that match between them, and ordering them by
            # the same sample IDs.
            (dm1_labels, dm1_data), (dm2_labels, dm2_data) = \
                make_compatible_distance_matrices((dm1_labels, dm1_data),
                                                  (dm2_labels, dm2_data), lookup=sample_id_map)
            if method == 'partial_mantel':
                # We need to intersect three sets (three matrices).
                (dm1_labels, dm1_data), (cdm_labels, cdm_data) = \
                    make_compatible_distance_matrices(
                        (dm1_labels, dm1_data), control_dm,
                        lookup=sample_id_map)
                (dm1_labels, dm1_data), (dm2_labels, dm2_data) = \
                    make_compatible_distance_matrices(
                        (dm1_labels, dm1_data), (dm2_labels, dm2_data),
                        lookup=sample_id_map)
                if len(dm1_labels) < 3:
                    result += '%s\t%s\t%s\t%d\tToo few samples\n' % (fp1,
                                                                     fp2, control_dm_fp, len(dm1_labels))
                    continue
            elif len(dm1_labels) < 3:
                result += '%s\t%s\t%d\tToo few samples\n' % (fp1, fp2,
                                                             len(dm1_labels))
                continue

            dm1 = DistanceMatrix(dm1_data, dm1_labels)
            dm2 = DistanceMatrix(dm2_data, dm2_labels)

            if method == 'mantel':
                corr_coeff, p_value, n = mantel(dm1, dm2, method='pearson',
                                 permutations=num_perms, alternative=tail_type,
                                 strict=True)
                p_str = p_value_to_str(p_value, num_perms)
                result += "%s\t%s\t%d\t%.5f\t%s\t%d\t%s\n" % (
                    fp1, fp2, n, corr_coeff, p_str, num_perms, tail_type)
            elif method == 'partial_mantel':
                cdm = DistanceMatrix(cdm_data, cdm_labels)
                results = PartialMantel(dm1, dm2, cdm)(num_perms)
                p_str = p_value_to_str(results['mantel_p'], num_perms)
                result += "%s\t%s\t%s\t%d\t%.5f\t%s\t%d\t%s\n" % (
                    fp1, fp2, control_dm_fp, len(dm1_labels),
                    results['mantel_r'], p_str, num_perms, 'greater')
    return result


def run_mantel_correlogram(fps, distmats, num_perms, comment, alpha,
                           sample_id_map=None,
                           variable_size_distance_classes=False):
    """Runs a Mantel correlogram analysis on all pairs of distance matrices.

    Returns a string suitable for writing out to a file containing the results
    of the test, a list of correlogram filepath names, and a list of matplotlib
    Figure objects representing each correlogram.

    The correlogram filepaths can have an extension string appended to the end
    of them and then be used to save each of the correlogram Figures to a file.
    Each correlogram filepath will be a combination of the two distance matrix
    filepaths that were used to create it.

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.

    Arguments:
        fps - list of filepaths of the distance matrices
        distmats - list of tuples containing dm labels and dm data (i.e. the
            output of parse_distmat)
        num_perms - the number of permutations to use to calculate the
            p-value(s)
        comment - comment string to add to the beginning of the results string
        alpha - the alpha value to use to determine significance in the
            correlogram plots
        sample_id_map - dict mapping sample IDs (i.e. what is expected by
            make_compatible_distance_matrices)
        variable_size_distance_classes - create distance classes that vary in
            size (i.e. width) but have the same number of distances in each
            class
    """
    if len(fps) != len(distmats):
        raise ValueError("Must provide the same number of filepaths as there "
                         "are distance matrices.")
    if comment is None:
        comment = ''
    result = comment + 'DM1\tDM2\tNumber of entries\t' + \
                       'Number of permutations\tClass index\t' + \
                       'Number of distances\tMantel r statistic\t' + \
                       'p-value\tp-value (Bonferroni corrected)\tTail type\n'
    correlogram_fps = []
    correlograms = []

    # Loop over all pairs of dms.
    for i, (fp1, (dm1_labels, dm1_data)) in enumerate(zip(fps, distmats)):
        for fp2, (dm2_labels, dm2_data) in zip(fps, distmats)[i + 1:]:
            # Make the current pair of distance matrices compatible by only
            # keeping samples that match between them, and ordering them by
            # the same sample IDs.
            (dm1_labels, dm1_data), (dm2_labels, dm2_data) = \
                make_compatible_distance_matrices((dm1_labels, dm1_data),
                                                  (dm2_labels, dm2_data), lookup=sample_id_map)
            if len(dm1_labels) < 3:
                result += '%s\t%s\t%d\tToo few samples\n' % (fp1, fp2,
                                                             len(dm1_labels))
                continue

            dm1 = DistanceMatrix(dm1_data, dm1_labels)
            dm2 = DistanceMatrix(dm2_data, dm2_labels)

            # Create an instance of our Mantel correlogram test and run it with
            # the specified number of permutations.
            mc = MantelCorrelogram(dm1, dm2, alpha=alpha,
                                   variable_size_distance_classes=variable_size_distance_classes)
            results = mc(num_perms)

            # Generate a name for the current correlogram and save it and the
            # correlogram itself.
            dm1_name = path.basename(fp1)
            dm2_name = path.basename(fp2)
            correlogram_fps.append('_'.join((dm1_name, 'AND', dm2_name,
                                             'mantel_correlogram')) + '.')
            correlograms.append(results['correlogram_plot'])

            # Iterate over the results and write them to the text file.
            first_time = True
            for class_idx, num_dist, r, p, p_corr in zip(
                    results['class_index'], results['num_dist'],
                    results['mantel_r'], results['mantel_p'],
                    results['mantel_p_corr']):
                # Format p-values and figure out which tail type we have based
                # on the sign of r.
                p_str = None
                if p is not None:
                    p_str = p_value_to_str(p, num_perms)
                p_corr_str = None
                if p_corr is not None:
                    p_corr_str = p_value_to_str(p_corr, num_perms)
                if r is None:
                    tail_type = None
                elif r < 0:
                    tail_type = 'less'
                else:
                    tail_type = 'greater'

                if first_time:
                    result += '%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\n' % (
                        fp1, fp2, len(dm1_labels), num_perms, class_idx,
                        num_dist, r, p_str, p_corr_str, tail_type)
                    first_time = False
                else:
                    result += '\t\t\t\t%s\t%d\t%s\t%s\t%s\t%s\n' % (class_idx,
                                                                    num_dist, r, p_str, p_corr_str, tail_type)
    return result, correlogram_fps, correlograms
