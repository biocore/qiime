#!/usr/bin/env python
# file make_3d_plots.py

__author__ = "Jesse Stombaugh, Rob Knight, and Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stombaugh", "Rob Knight", "Micah Hamady", "Dan Knights",
               "Antonio Gonzalez Pena", "Yoshiki Vazquez Baeza", "Will Van Treuren"]
    # remember to add yourself
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

import os
import re
import numpy as np
from copy import deepcopy
from random import choice
from time import strftime
from numpy.linalg import norm
from biplots import make_mage_taxa
from numpy import abs as numpy_abs
from cogent.util.misc import flatten
from qiime.format import format_coords
from cogent.maths.stats.util import Numbers
from qiime.sort import natsort, signed_natsort
from qiime.pycogent_backports.test import ANOVA_one_way
from qiime.parse import parse_coords, group_by_field, parse_mapping_file
from qiime.util import load_pcoa_files, summarize_pcoas, MissingFileError
from numpy import (array, shape, apply_along_axis, dot, delete, vstack, sqrt,
                   average, isnan, nan, diff, mean, std, concatenate, ones,
                   append, zeros)
from qiime.colors import (get_group_colors, color_groups, make_color_dict,
                          combine_map_label_cols, process_colorby,
                          linear_gradient, iter_color_groups, get_map,
                          kinemage_colors)


'''
xdata_colors = {
        'aqua':     (180, 100, 100),
        'blue':     (240,100,100),
        'fuchsia':  (300,100,100),
        'gray':     (300,0,50.2),
        'green':    (120,100,50.2),
        'lime':     (120,100,100),
        'maroon':   (0,100,50.2),
        'olive':    (60,100,50.2),
        'purple':   (300,100,50.2),
        'red':      (0,100,100),
        'silver':   (0, 0, 75.3),
        'teal':     (180,100,50.2),
        'yellow':   (60,100,100),
}
'''

data_colors = {'blue': '#0000FF', 'lime': '#00FF00', 'red': '#FF0000',
               'aqua': '#00FFFF', 'fuchsia': '#FF00FF', 'yellow': '#FFFF00',
               'green': '#008000', 'maroon': '#800000', 'teal': '#008080',
               'purple': '#800080', 'olive': '#808000',
               'silver': '#C0C0C0', 'gray': '#808080'}


def create_dir(dir_path, plot_type):
    """Creates directory where data is stored.  If directory is not supplied in\
       the command line, a random folder is generated"""

    alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
    alphabet += alphabet.lower()
    alphabet += "01234567890"

    if dir_path is None or dir_path == '':
        dir_path = ''
        random_dir_name = ''.join([choice(alphabet) for i in range(10)])
        dir_path = './' + plot_type + \
            strftime("%Y_%m_%d_%H_%M_%S") + random_dir_name + '/'

    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    return dir_path


def make_3d_plots(coord_header, coords, pct_var, mapping, prefs,
                  background_color, label_color, taxa=None, custom_axes=None,
                  edges=None, coords_low=None, coords_high=None,
                  ellipsoid_prefs=None, user_supplied_edges=False,
                  ball_scale=1.0,
                  arrow_colors={'line_color': 'white', 'head_color': 'red'},
                  add_vectors=None, plot_scaled=False, plot_unscaled=True):
    """Makes 3d plots given coords, mapping file, and prefs.

    Added quick-and-dirty hack for gradient coloring of columns, should
    replace with more general method. Current solution is to pass in any
    of the following:

    'colors':(('white', (0,100,100)),('red',(100,100,100)))

    makes gradient between white and red, applies to all samples

    'colors':{'RK':(('white',(0,0,100)),('red',(0,100,100))),
              'NF':(('white',(120,0,100)),('green',(120,100,100)))
             }
    pulls the combination samples starting with RK, colors with
    first gradient, then pulls the combination samples starting
    with NF, colors with the next gradient.

    """
    sorting_column = {}
    result = []
    # Iterate through prefs and color by given mapping labels
    # Sort by the column name first
    groups_and_colors = iter_color_groups(mapping, prefs)
    groups_and_colors = list(groups_and_colors)

    # processing the vectors
    if add_vectors:

        # getting the index of the columns in the mapping file
        ind_group = mapping[0].index(add_vectors['vectors'][0])
        if len(add_vectors['vectors']) >= 2:
            ind_sort = mapping[0].index(add_vectors['vectors'][1])
        else:
            ind_sort = 0

        # creating groups
        groups = {}
        for l in mapping[1:]:
            if l[0] not in coord_header:
                continue
            if l[ind_group] not in groups:
                groups[l[ind_group]] = [(l[ind_sort], l[0])]
            else:
                groups[l[ind_group]].append((l[ind_sort], l[0]))
        for g in groups:
            groups[g] = signed_natsort(groups[g])
        add_vectors['vectors'] = groups

        # managing the options to weight the vectors by a given metadata column
        if add_vectors['weight_by_vector'] and ind_sort == 0:
            raise ValueError("The weighting_vector can't be the SampleID, " +
                             "please specify a numeric column in the --add_vectors option.")
        elif add_vectors['weight_by_vector']:
            # if we can sort the column, try to create a dictionary with the
            # value of this sorting column as values and the SampleID as key
            for mapping_row in mapping[1:]:
                try:
                    sorting_column[mapping_row[0]] = float(
                        mapping_row[ind_sort])
                except ValueError:
                    raise ValueError("The sorting column for --add_vectors " +
                                     "must be numeric, i. e. DOB or days_since_epoch")

    for i in range(len(groups_and_colors)):
        # Write to kinemage file using the groups, colors and coords
        labelname = groups_and_colors[i][0]
        groups = groups_and_colors[i][1]
        colors = groups_and_colors[i][2]
        data_colors = groups_and_colors[i][3]
        data_color_order = groups_and_colors[i][4]

        if plot_unscaled:
            result.extend(
                make_mage_output(groups, colors, coord_header, coords,
                                 pct_var, background_color, label_color, data_colors,
                                 taxa, custom_axes, name=labelname,
                                 scaled=False, edges=edges,
                                 coords_low=coords_low, coords_high=coords_high,
                                 ellipsoid_prefs=ellipsoid_prefs,
                                 user_supplied_edges=user_supplied_edges,
                                 ball_scale=ball_scale, arrow_colors=arrow_colors,
                                 add_vectors=add_vectors, sorting_column=sorting_column))
        if plot_scaled:
            result.extend(
                make_mage_output(groups, colors, coord_header, coords,
                                 pct_var, background_color, label_color, data_colors,
                                 taxa, custom_axes, name=labelname,
                                 scaled=True, edges=edges,
                                 coords_low=coords_low, coords_high=coords_high,
                                 ellipsoid_prefs=ellipsoid_prefs,
                                 user_supplied_edges=user_supplied_edges,
                                 ball_scale=ball_scale, arrow_colors=arrow_colors,
                                 add_vectors=add_vectors, sorting_column=sorting_column))

    return result


def scale_pc_data_matrix(coords, pct_var):
    """Scales pc data matrix by percent variation"""
    return coords * (pct_var / pct_var.max())


def auto_radius(coords, ratio=0.01):
    """Determine radius from coords"""
    dim1 = coords[:, 0]
    range = dim1.max() - dim1.min()

    return ratio * range


def make_mage_output(groups, colors, coord_header, coords, pct_var,
                     background_color, label_color, data_colors,
                     taxa=None, custom_axes=None, name='',
                     radius=None, alpha=.75, num_coords=10, scaled=False,
                     coord_scale=1.05, edges=None, coords_low=None,
                     coords_high=None, ellipsoid_prefs=None,
                     user_supplied_edges=False, ball_scale=1.0,
                     arrow_colors={'line_color': 'white', 'head_color': 'red'},
                     add_vectors=None, sorting_column=None):
    """Convert groups, colors, coords and percent var into mage format
        weighting_colum: is a dictionary with SampleIDs as keys and the
        corresponding values for the column selected to sort the vectors.
    """
    result = []

    # Scale the coords and generate header labels
    if scaled:
        scalars = pct_var
        if custom_axes:
            # create a dummy vector of ones to avoid scaling custom axes
            custom_scalars = scalars[0] * np.ones(len(custom_axes))
            scalars = np.append(custom_scalars, scalars)
        coords = scale_pc_data_matrix(coords, scalars)
        if not coords_low is None:
            coords_low = scale_pc_data_matrix(coords_low, scalars)
        if not coords_high is None:
            coords_high = scale_pc_data_matrix(coords_high, scalars)
        header_suffix = '_scaled'
    else:
        header_suffix = '_unscaled'

    if radius is None:
        radius = float(auto_radius(coords)) * float(ball_scale)
    else:
        radius = float(radius) * float(ball_scale)

    maxes = coords.max(0)[:num_coords]
    mins = coords.min(0)[:num_coords]
    pct_var = pct_var[:num_coords]  # scale from fraction

    # check that we didn't get fewer dimensions than we wanted
    if len(mins) < num_coords:
        num_coords = len(mins)
    min_maxes = flatten(zip(mins, maxes))

    if custom_axes:
        axis_names = ['PC%s' % (i + 1)
                      for i in xrange(num_coords - len(custom_axes))]
        axis_names = custom_axes + axis_names
    else:
        axis_names = ['PC%s' % (i + 1) for i in xrange(num_coords)]

    # Write the header information
    result.append('@kinemage {%s}' % (name + header_suffix))
    result.append(
        '@dimension ' + ' '.join(['{%s}' %
                                (aname) for aname in axis_names]))
    result.append('@dimminmax ' + ' '.join(map(str, min_maxes)))
    result.append('@master {points}')
    result.append('@master {labels}')
    if edges:
        result.append('@master {edges}')

    if not taxa is None:
        result.append('@master {taxa_points}')
        result.append('@master {taxa_labels}')

    for cname, color in sorted(data_colors.items()):
        result.append(color.toMage())

    if background_color == 'white':
        result.append('@whitebackground')
        result.append('@hsvcolor {black} 0.0 0.0 0.0')
    else:
        result.append('@hsvcolor {white} 180.0 0.0 100.0')

    # Write the groups, colors and coords
    coord_dict = dict(zip(coord_header, coords))
    if not coords_low is None:
        coord_low_dict = dict(zip(coord_header, coords_low))
    if not coords_high is None:
        coord_high_dict = dict(zip(coord_header, coords_high))
    for group_name in natsort(groups):
        ids = groups[group_name]
        result.append(
            '@group {%s (n=%s)} collapsible' %
            (group_name, len(ids)))

        color = colors[group_name]
        coord_lines = []
        for id_ in sorted(ids):
            if id_ in coord_dict:
                coord_lines.append('{%s} %s' %
                                   (id_, ' '.join(map(str, coord_dict[id_][:num_coords]))))

        # create list of balls, one for each sample
        result.append('@balllist color=%s radius=%s alpha=%s dimension=%s \
master={points} nobutton' % (color, radius, alpha, num_coords))
        result.append('\n'.join(coord_lines))
        # make ellipsoids if low and high coord bounds were received
        if (not coords_low is None) and (not coords_high is None):
            # create one trianglelist for each sample to define ellipsoids
            result += make_mage_ellipsoids(ids, coord_dict, coord_low_dict,
                                           coord_high_dict, color, ellipsoid_prefs)

        # create list of labels
        result.append('@labellist color=%s radius=%s alpha=%s dimension=%s \
master={labels} nobutton' % (color, radius, alpha, num_coords))
        result.append('\n'.join(coord_lines))

        # Write vectors if requested
        if add_vectors:
            result += make_vectors_output(coord_dict,
                                          add_vectors, num_coords, color, ids)

            # Calculate the average traces for this group
            avg_coord_dict, avg_add_vectors = avg_vector_for_group(group_name,
                                                                   ids, coord_dict, custom_axes, deepcopy(add_vectors))

            # Add the average traces, these are by default non-visible
            result += make_vectors_output(avg_coord_dict, avg_add_vectors,
                                          num_coords, color, avg_coord_dict.keys(
                                          ),
                                          display_name='average', state='off')

            avg_coord_dict = {}
            avg_add_vectors = {}
            avg_ids = {}

            if not scaled and add_vectors['vectors_path'] and add_vectors['vectors_algorithm']:

                # the sorting_column dict is optional, so get it if possible
                try:
                    # for the current group of ids, get the needed values, only
                    # if they are in coord_dict
                    weighting_vector = array(
                        [sorting_column[single_id] for single_id in ids if single_id in coord_dict])
                except KeyError:
                    weighting_vector = None

                try:
                    vector_result = make_subgroup_vectors(
                        coord_dict, add_vectors['eigvals'],
                        ids, add_vectors[
                            'vectors_algorithm'],
                        add_vectors[
                            'vectors_axes'], custom_axes,
                        weight=add_vectors[
                            'weight_by_vector'],
                        weighting_vector=weighting_vector,
                        window_size=add_vectors['window_size'])
                except TypeError:
                    continue

                if isinstance(vector_result['calc'], dict) or not isnan(vector_result['calc']):
                    if name not in add_vectors['vectors_output']:
                        add_vectors['vectors_output'][name] = {}
                    if group_name not in add_vectors['vectors_output'][name]:
                        add_vectors['vectors_output'][name][group_name] = {}

                    add_vectors[
                        'vectors_output'][
                        name][
                        group_name][
                        'vectors_vector'] = vector_result[
                        'vector']
                    add_vectors[
                        'vectors_output'][
                        name][
                        group_name][
                        'vectors_result'] = vector_result[
                        'calc']
                    add_vectors[
                        'vectors_output'][
                        name][
                        group_name][
                        'vectors_message'] = vector_result[
                        'message']

    if not taxa is None:
        result += make_mage_taxa(taxa, num_coords, pct_var,
                                 scaled=scaled, scalars=None, radius=radius)

    # Write the axes on the bottom of the graph
    result.append('@group {axes} collapsible')
    state = 'on'
    axis_mins = mins * coord_scale
    axis_maxes = maxes * coord_scale

    if not custom_axes:
        custom_axes = []
    # draw each axis
    for i in xrange(num_coords):
        if i == 3:
            state = 'off'
        result.append('@vectorlist {%s line} dimension=%s %s' %
                      (axis_names[i], num_coords, state))

        result.append(' '.join(map(str, axis_mins)) + ' ' + label_color)
        end = axis_mins.copy()
        end[i] = axis_maxes[i]
        result.append(' '.join(map(str, end)) + ' ' + label_color)
        end[i] *= coord_scale  # add scale factor to offset labels a little

        # custom axes come first, no "percent variance" shown
        if i < len(custom_axes):
            result.append('@labellist {%s} dimension=%s %s' %
                          (axis_names[i], num_coords, state))
            result.append(('{%s}' % (axis_names[i])) +
                          ' '.join(map(str, end)) + ' ' + label_color)
        # if all custom axes have been drawn, draw normal PC axes
        else:
            pct = pct_var[i - len(custom_axes)]
            result.append('@labellist {%s (%0.2g%%)} dimension=%s %s' %
                          (axis_names[i], pct, num_coords, state))
            result.append(('{%s (%0.2g%%)}' % (axis_names[i], pct)) +
                          ' '.join(map(str, end)) + ' ' + label_color)

    # Write edges if requested
    if edges:
        result += make_edges_output(coord_dict, edges, num_coords, label_color,
                                    arrow_colors=arrow_colors,
                                    user_supplied_edges=user_supplied_edges)

    return result


def make_subgroup_vectors(
        coord_dict, eigvals, ids, method='avg', vectors_axes=3,
        custom_axes=None, weight=False, weighting_vector=None,
        window_size=None):
    """Creates either the first difference or the rms value, vector, of a subgroup (ids) of coord_dict

       Params:
        coord_dict: a dict of (sampleID, coords), where coords is a numpy array (row)
        eigvals: the eigvals of the coords
        ids: the list of ids that should be used to caclulate the first difference or the rms
        from coord_dict
        method: 'avg', 'trajectory', 'diff' or 'wdiff'
        custom_axes: the list of custom axis used in other methods
        weight: whether or not the vector should be weighted using the
        qiime.make_3d_plots.weight_by_vector function
        window_size: required only if 'wdiff' is selected as a method.

       Returns:
        a dict with 2 members: 'vector', the vector representing from which the first
        difference or RMS value was calculated; and 'calc' the resulting calculation of the first
        difference or the RMS value with the selected algorithm.
    """
    result = {}
    message_buffer = ''

    # We multiply the coord values with the value of the eigvals represented
    if custom_axes:
        vectors = [coord_dict[id][1:][:vectors_axes] * eigvals[:vectors_axes]
                   for id in ids if id in coord_dict]
    else:
        vectors = [coord_dict[id][:vectors_axes] * eigvals[:vectors_axes]
                   for id in ids if id in coord_dict]

    if len(vectors) <= 0:
        raise TypeError(
            "No samples to process, an empty list cannot be processed")

    # the weighting can only be done over vectors with a lenght greater than 1
    if weight and weighting_vector is not None and len(ids) > 1:
        try:
            vectors_copy = deepcopy(vectors)
            vectors = weight_by_vector(vectors_copy, weighting_vector)
        except (FloatingPointError, ValueError):
            message_buffer = "Could not weight group, no gradient in the " +\
                "the weighting vector.\n"
            vectors = vectors_copy

    elif weight and weighting_vector is None:
        raise ValueError("You must pass a weighting vector if you want to" +
                         "weight your data")

    # to check for the length, look at vectors, ids might be misleading
    if method == 'avg':
        center = average(vectors, axis=0)
        if len(vectors) == 1:
            result['vector'] = [norm(center)]
            result['calc'] = {'avg': result['vector']}
        else:
            result['vector'] = [norm(i) for i in vectors - center]
            result['calc'] = {'avg': average(result['vector'])}

    elif method == 'trajectory':
        if len(vectors) == 1:
            result['vector'] = [norm(vectors)]
            result['calc'] = {'trajectory': result['vector']}
        else:
            result['vector'] = [norm(vectors[i - 1] - vectors[i])
                                for i in range(len(vectors) - 1)]
            result['calc'] = {'trajectory': norm(result['vector'])}

    elif method == 'diff':
        if len(vectors) == 1:
            result['vector'] = [norm(vectors)]
            result['calc'] = {'mean': result['vector'], 'std': 0}
        elif len(vectors) == 2:
            result['vector'] = [norm(vectors[1] - vectors[0])]
            result['calc'] = {'mean': result['vector'], 'std': 0}
        else:
            vec_norm = [norm(vectors[i - 1] - vectors[i])
                        for i in range(len(vectors) - 1)]
            vec_diff = diff(vec_norm)
            result['calc'] = {'mean': mean(vec_diff), 'std': std(vec_diff)}
            result['vector'] = vec_diff
    elif method == 'wdiff':
        if len(vectors) == 1:
            result['vector'] = [norm(vectors)]
            result['calc'] = {'mean': result['vector'], 'std': 0}
        elif len(vectors) == 2:
            result['vector'] = [norm(vectors[1] - vectors[0])]
            result['calc'] = {'mean': result['vector'], 'std': 0}
        else:
            vec_norm = [norm(vectors[i - 1] - vectors[i])
                        for i in range(len(vectors) - 1)]

            # windowed first differences won't be able on every group, specially
            # given the variation of size that a vector tends to have
            try:
                vec_diff = windowed_diff(vec_norm, window_size)
            except ValueError:
                vec_diff = vec_norm
                message_buffer = message_buffer + 'Cannot calculate the ' +\
                    'first difference with a window of size (%d).\n' % window_size

            result['calc'] = {'mean': mean(vec_diff), 'std': std(vec_diff)}
            result['vector'] = vec_diff

    elif method is None:
        result['calc'] = nan
    else:
        raise ValueError('The method "%s" for Vectors does not exist' % method)

    # If there is no message to write, then just set it to None
    if message_buffer == '':
        message_buffer = None
    result['message'] = message_buffer

    return result


def run_ANOVA_trajetories(groups):
    """Run ANOVA on trajectories categories

       Added some lines to control for labels and make everything a list vs. list/sets

       Params:
        groups, a dict of {trajectory_id: values}, where trajectory_id is
        each of the categories to test, values is a numpy array with the edges of
        that trajectory

       Returns:
        look at otu_category_significance.py -> run_single_ANOVA

        This was taken from qiime/otu_category_significance.py:run_single_ANOVA
    """
    values = []
    labels = []
    for trajectory in groups:
        labels.append(trajectory)
        values.append(groups[trajectory])

    # This is copied from otu_category_significance.py
    try:
        F, prob = ANOVA_one_way(values)
        group_means = [i.mean() for i in values]
    except ValueError:
        # set the p-value to 'diff' if the variances are 0.0 (within rounding
        # error) and the means are not all the same. If the means are all
        # the same and the variances are 0.0, set the p-value to 1
        group_means = []
        group_variances = []
        for i in values:
            group_means.append(i.Mean)
            group_variances.append(i.Variance)
        # Added list to always return the same type of output
        group_means = list(set(group_means))
        if sum(group_variances) < 1e-21 and len(group_means) > 1:
            prob = 0.0
        else:
            prob = 1.0

    return labels, group_means, prob


def make_vectors_output(coord_dict, add_vectors, num_coords,
                        color, ids, display_name='trace', state='on'):
    """Creates make output to display vectors (as a kinemage vectorlist).

       Params:
        coord_dict, a dict of (sampleID, coords), where coords is a numpy
        array edges, a dictionary of pairs of samples where
        ['label': [('order': sampleID) ...]
        num_coords, the number of included dimensions in the PCoA plot
        color, the color of this batch of vectors,
        ids the list of ids that are in this batch and where the first member
        of the groups id a should be a member of.
        state: whether the trace displayed will appear on screen by default,
        the permitted values are 'on' and 'off'; default is 'on'.

       Returns:
        result, a list of strings containing a kinemage vectorlist
    """
    # default parameters
    which_set = 0
    tip_fraction = 0.5
    arrow_colors = {'line_color': color, 'head_color': color}

    # header for section in kin file
    result = []
    result.append(
        '@vectorlist {%s} dimension=%s %s' %
        (display_name, num_coords, state))

    # selecting the edges to draw in this batch
    edges = []
    for k, v in add_vectors['vectors'].iteritems():
        edges = edges + [(v[i][1], v[i + 1][1])
                         for i in range(len(v[:-1])) if v[i][1] in ids]

    # creating "vectors"
    for edge in edges:
        id_fr, id_to = edge
        # extract the coords of each vertex
        pt_fr = coord_dict[id_fr][:num_coords]
        pt_to = coord_dict[id_to][:num_coords]

        # different tip color for each destination coords file
        #tip_color = kinemage_colors[which_set % len(kinemage_colors)]
        # plot a color 'tip' on the line (certain % of line length)
        # this gets the coords of the beginning of the 'tip'
        diffs = (pt_to - pt_fr) * (1 - tip_fraction)
        middles = pt_fr + diffs
        # add a default-color line segment

        # modified to use user defined
        tip_color = arrow_colors['head_color']
        label_color = arrow_colors['line_color']

        result.append('%s %s' %
                      (' '.join(map(str, pt_fr)), label_color))
        result.append('%s %s P' %
                      (' '.join(map(str, middles)), label_color))
        # add the tip-colored line segment
        result.append('%s %s' %
                      (' '.join(map(str, middles)), tip_color))
        result.append('%s %s P' %
                      (' '.join(map(str, pt_to)), tip_color))
    return result


def make_edges_output(coord_dict, edges, num_coords, label_color, arrow_colors,
                      tip_fraction=0.4, user_supplied_edges=False):
    """Creates make output to display edges (as a kinemage vectorlist).

       Params:
        coord_dict, a dict of (sampleID, coords), where coords is a numpy array
        edges, a list of pairs of sampleIDs (from, to)
        num_coords, the number of included dimensions in the PCoA plot
        label_color, the plain edge color.
        tip_fraction, the portion of each edge to be colored as the 'tip'

       Returns:
        result, a list of strings containing a kinemage vectorlist
    """
    result = []
    result.append(
        '@vectorlist {edges} dimension=%s on master={edges} nobutton' %
        (num_coords))
    for edge in edges:
        id_fr, id_to = edge
        # extract the coords of each vertex
        pt_fr = coord_dict[id_fr][:num_coords]
        pt_to = coord_dict[id_to][:num_coords]
        # get 'index' of the destination coords file from 'to' sampleID
        if user_supplied_edges:
            which_set = 0
        else:
            which_set = int(id_to[id_to.rindex('_') + 1:]) - 1

        # different tip color for each destination coords file
        #tip_color = kinemage_colors[which_set % len(kinemage_colors)]
        # plot a color 'tip' on the line (certain % of line length)
        # this gets the coords of the beginning of the 'tip'
        diffs = (pt_to - pt_fr) * (1 - tip_fraction)
        middles = pt_fr + diffs
        # add a default-color line segment

        # modified to use user defined
        tip_color = arrow_colors['head_color']
        label_color = arrow_colors['line_color']

        result.append('%s %s' %
                      (' '.join(map(str, pt_fr)), label_color))
        result.append('%s %s P' %
                      (' '.join(map(str, middles)), label_color))
        # add the tip-colored line segment
        result.append('%s %s' %
                      (' '.join(map(str, middles)), tip_color))
        result.append('%s %s P' %
                      (' '.join(map(str, pt_to)), tip_color))
    return result


def make_mage_ellipsoids(ids, coord_dict, coord_low_dict,
                         coord_high_dict, color, ellipsoid_prefs=
                         {"smoothness": 2, "alpha": .25}):
    """Makes ellipsoids with centers in coord_dict.
       coord_low_dict and coord_high_dict are used to scale
       each axis of the ellipsoid.
    """
    alpha = ellipsoid_prefs['alpha']
    nsubdivs = ellipsoid_prefs['smoothness']
    result = []
    coord_lines = []
    for id_ in sorted(ids):
        if id_ in coord_dict:
            center = coord_dict[id_][:3]
            dims = coord_high_dict[id_][:3] - coord_low_dict[id_][:3]

            faces = make_ellipsoid_faces(center, dims, nsubdivs=nsubdivs)
            for face in faces:
                result.append(
                    "@trianglelist color=%s alpha=%f master={points} nobutton" %
                    (color, alpha))
                for point in face:
                    result.append(' '.join(map(str, point)))
    return result


def make_ellipsoid_faces(center, dims, nsubdivs=2):
    """Returns a list of 3-tuples (triangles) of 3-tuples (points)
       defining an ellipsoid centered at center with axis
       dimensions given in dims.

       nsubdivs determines the number of recursive divisions of
       the faces that will be made.
       nsubdivs value     No. faces
       0                  20
       1

    """
    t = (1 + sqrt(5)) / 2.0
    s = sqrt(1 + t ** 2)

    vertices = [(t / s, 1 / s, 0), (-t / s, 1 / s, 0), (t / s, -1 / s, 0),
                (-t / s, -1 / s, 0), (1 / s, 0, t / s), (1 / s, 0,
                                                         -t / s), (
                    -1 / s, 0, t / s), (
                    -1 / s, 0, -t / s),
                (0, t / s, 1 / s), (0, -t / s, 1 / s), (0, t / s, -1 / s), (0, -t / s, -1 / s)]

    v = vertices
    faces = [(
        v[0], v[8], v[4]), (
        v[1], v[10], v[7]), (
        v[2], v[9], v[11]), (
        v[7], v[3], v[1]), (v[0], v[5], v[10]), (v[3], v[9], v[6]),
        (v[3], v[11], v[9]), (v[8], v[6], v[4]), (v[2], v[4],
                                                  v[9]), (
            v[3], v[7], v[11]), (
            v[4], v[2], v[0]),
        (v[9], v[4], v[6]), (v[2], v[11], v[5]), (v[0], v[10], v[8]
                                                  ), (v[
            5], v[0], v[2]), (v[10], v[5], v[7]), (v[1], v[6], v[8]),
        (v[1], v[8], v[10]), (v[6], v[1], v[3]), (v[11], v[7], v[5])]

    # subdivide each of the faces into 9 faces
    for i in xrange(nsubdivs):
        new_faces = []
        for face in faces:
            new_faces.extend(subdivide(face[0], face[1], face[2]))
        faces = new_faces
    faces = scale_faces(dims[0], dims[1], dims[2], faces)
    faces = translate_faces(center, faces)
    return faces


def subdivide(x, y, z):
    # look at x-y edge
    xy = [(
        x[0] * 2 / 3.0 + y[0] / 3.0,
        x[1] * 2 / 3.0 + y[1] / 3.0,
        x[2] * 2 / 3.0 + y[2] / 3.0),
        (x[0] / 3.0 + y[0] * 2 / 3.0,
         x[1] / 3.0 + y[1] * 2 / 3.0,
            x[2] / 3.0 + y[2] * 2 / 3.0)]
    # pull them to the surface of the sphere
    xy = [(
        i[0] / sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2),
        i[1] / sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2),
        i[2] / sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2)) for i in xy]

    # do the same for the other edges
    xz = [(
        x[0] * 2 / 3.0 + z[0] / 3.0,
        x[1] * 2 / 3.0 + z[1] / 3.0,
        x[2] * 2 / 3.0 + z[2] / 3.0),
        (x[0] / 3.0 + z[0] * 2 / 3.0,
         x[1] / 3.0 + z[1] * 2 / 3.0,
            x[2] / 3.0 + z[2] * 2 / 3.0)]
    xz = [(
        i[0] / sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2),
        i[1] / sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2),
        i[2] / sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2)) for i in xz]

    zy = [(
        z[0] * 2 / 3.0 + y[0] / 3.0,
        z[1] * 2 / 3.0 + y[1] / 3.0,
        z[2] * 2 / 3.0 + y[2] / 3.0),
        (z[0] / 3.0 + y[0] * 2 / 3.0,
         z[1] / 3.0 + y[1] * 2 / 3.0,
            z[2] / 3.0 + y[2] * 2 / 3.0)]
    zy = [(
        i[0] / sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2),
        i[1] / sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2),
        i[2] / sqrt(i[0] ** 2 + i[1] ** 2 + i[2] ** 2)) for i in zy]

    center = (
        (x[0] + y[0] + z[0]) / 3.0,
        (x[1] + y[1] + z[1]) / 3.0,
        (x[2] + y[2] + z[2]) / 3.0)
    center_len = sqrt(center[0] ** 2 + center[1] ** 2 + center[2] ** 2)
    center = (
        center[0] / center_len,
        center[1] / center_len,
        center[2] / center_len)

    # generate the new list of faces
    faces = [(
        x, xz[0], xy[0]), (
        xz[0], xy[0], center), (
        xz[0], xz[1], center), (
        xz[1], zy[0], center), (
        xz[1], z, zy[0]), (xy[0], center, xy[1]),
        (xy[1], y, zy[1]), (xy[1], zy[1], center), (zy[1], zy[0], center)]

    return faces


def scale_faces(a, b, c, faces):
    for i, face in enumerate(faces):
        faces[i] = list(faces[i])
        for j, point in enumerate(face):
            faces[i][j] = (point[0] * a, point[1] * b, point[2] * c)
    return faces


def translate_faces(center, faces):
    for i, face in enumerate(faces):
        faces[i] = list(faces[i])
        for j, point in enumerate(face):
            faces[i][j] = (
                point[0] + center[0],
                point[1] + center[1],
                point[2] + center[2])
    return faces


def process_custom_axes(axis_names):
    """Parses the custom_axes option from the command line"""
    return axis_names.strip().strip("'").strip('"').split(',')


def process_coord_filenames(coord_filenames):
    """Parses the custom_axes option from the command line"""
    return coord_filenames.strip().strip("'").strip('"').split(',')


def validate_coord_files(coord_files):
    """Returns True if there are the same number of column headers,
       coords, eigvals, and pct var's in every coord file"""
    # if coord_filenames is a string, and it's a directory, get list of files
    if isinstance(coord_files, str):
        if os.path.isdir(coord_files):
            dirname = coord_files
            coord_files = os.listdir(dirname)
            coord_files = [os.path.join(dirname, fn) for fn in coord_files]
        else:
            coord_files = [coord_files]

    # initialize n_columns using the first line of the first file
    n_columns = len(open(coord_files[0], 'U')
                    .readlines()[0].strip().split('\t'))
    for fn in coord_files:
        lines = open(fn, 'U').readlines()
        for line in lines:
            line = line.strip()
            if len(line) > 0 and len(line.split('\t')) != n_columns:
                return False
    return True


def get_custom_coords(axis_names, mapping, coords):
    """Gets custom axis coords from the mapping file.
       Appends custom as first column(s) of PCoA coords matrix.

       Params:
        axis_names, the names of headers of mapping file columns
        mapping, the mapping file object (with list of headers in element 0)
        coords, the PCoA coords object, with coords matrix in element 1
    """
    for i, axis in enumerate(reversed(axis_names)):
        if not axis in mapping[0]:
            raise ValueError('Warning: could not find custom axis %s in map headers: %s'
                             % (axis, mapping[0]))
        else:
            # get index of column in mapping file
            col_idx = mapping[0].index(axis)
            # extract column data
            col = zip(*mapping[1:])[col_idx]
            sample_IDs = zip(*mapping[1:])[0]
            new_coords = array([])
            # load custom coord for this axis for each sample ID
            for id in coords[0]:
                if id in sample_IDs:
                    row_idx = list(sample_IDs).index(id)
                    try:
                        as_float = float(col[row_idx])
                        new_coords = np.append(new_coords, as_float)
                    except ValueError:
                        new_coords = np.append(new_coords, np.nan)
            new_coords = np.transpose(np.column_stack(new_coords))
            # append new coords to beginning column of coords matrix
            coords[1] = np.hstack((new_coords, coords[1]))


def remove_nans(coords):
    """Deletes any samples with NANs in their coordinates"""
    s = np.apply_along_axis(sum, 1, np.isnan(coords[1])) == 0
    coords[0] = (np.asarray(coords[0])[s]).tolist()
    coords[1] = coords[1][s, :]


def scale_custom_coords(custom_axes, coords):
    """Scales custom coordinates to match min/max of PC1"""

    # the target min and max
    to_mn = min(coords[1][:, len(custom_axes)])
    to_mx = 2 * max(coords[1][:, len(custom_axes)])

    # affine transformation for each custom axis
    for i in xrange(len(custom_axes)):
        from_mn = min(coords[1][:, i])
        from_mx = max(coords[1][:, i])
        coords[1][:, i] = (coords[1][:, i] - from_mn) / (from_mx - from_mn)
        coords[1][:, i] = (coords[1][:, i]) * (to_mx - to_mn) + to_mn

# The following functions were not unit_tested, however the parts within
# the functions are unit_tested


def get_sample_ids(maptable):
    """Extracts list of sample IDs from mapping file."""
    return [line[0] for line in maptable[1:]]


def get_coord(coord_fname, method="IQR"):
    """Opens and returns coords location matrix and metadata.
       Also two spread matrices (+/-) if passed a dir of coord files.
       If only a single coord file, spread matrices are returned as None.
    """
    if not os.path.isdir(coord_fname):
        try:
            coord_f = open(coord_fname, 'U').readlines()
        except (TypeError, IOError):
            raise MissingFileError('Coord file required for this analysis')
        coord_header, coords, eigvals, pct_var = parse_coords(coord_f)
        return [coord_header, coords, eigvals, pct_var, None, None]
    else:
        master_pcoa, support_pcoas = load_pcoa_files(coord_fname)

        # get Summary statistics
        coords, coords_low, coords_high, eigval_average, coord_header = \
            summarize_pcoas(master_pcoa, support_pcoas, method=method)
        pct_var = master_pcoa[3]  # should be getting this from an average

        # make_3d_plots expects coord_header to be a python list
        coord_header = list(master_pcoa[0])
        return (
            [coord_header,
             coords,
             eigval_average,
             pct_var,
             coords_low,
             coords_high]
        )


def get_multiple_coords(coord_fnames, edges_file=None, serial=False):
    """Opens and returns coords data and edges from multiple coords files.

       Params:
        coord_fnames, the names of the coordinate files

       Returns:
        edges, a list of pairs of sample IDs, (from, to)
        coords
            a list of [coord_header, coords, eigvals, pct_var]
            all coords are put in a single data matrix.
            Sample IDs from ith file have _i appended to them.
            eigvals, pct_var are taken from first coords file

       If "serial" is True, connects points ending with _0 to those with _1,
       those with _1 to those with _2, etc. Otherwise all sets are connected
       back to those ending with _0.
    """
    # start with empty data structures
    coord_header = []
    coords = []
    edges = []

    # load predetermined edges if they were passed to us
    if not edges_file is None:
        edges = [ln.strip().split()
                 for ln in open(edges_file, 'U').readlines()]

    # load all coords files into same data matrix
    for i, f in enumerate(coord_fnames):
        try:
            coord_f = open(coord_fnames[i], 'U').readlines()
        except (TypeError, IOError):
            raise MissingFileError('Coord file required for this analysis')
        coord_header_i, coords_i, eigvals_i, pct_var_i = parse_coords(coord_f)
        sampleIDs = coord_header_i
        # append _i to this file's sampleIDs unless we have predetermined edges
        if edges_file is None:
            coord_header_i = ['%s_%d' % (h, i) for h in coord_header_i]

        # get eigvals, pct_var from first coords file
        if i == 0:
            eigvals = eigvals_i
            pct_var = pct_var_i
            coord_header = coord_header_i
            coords = coords_i
        # for second, third, etc coords files, just append to first file
        else:
            coord_header.extend(coord_header_i)
            coords = vstack((coords, coords_i))
    # add all edges unless we have predetermined edges
    if edges_file is None:
        for _id in sampleIDs:
            if serial:
                for i in xrange(len(coord_fnames) - 1):
                    # edges go from one set to the next
                    edges += [('%s_%d' % (_id, i), '%s_%d' % (_id, i + 1))]
            else:
                for i in xrange(1, len(coord_fnames)):
                    # edges go from first file's points to other files' points
                    edges += [('%s_%d' % (_id, 0), '%s_%d' % (_id, i))]

    return edges, [coord_header, coords, eigvals, pct_var, None, None]


def get_taxa(taxa_fname, sample_ids):
    """Opens and returns coords data"""
    try:
        lines = open(taxa_fname, 'U').readlines()
    except (TypeError, IOError):
        raise MissingFileError('Taxa summary file required for this analysis')
    map = parse_mapping_file(lines)
    return map


def remove_unmapped_samples(mapping, coords, edges=None):
    """Removes any samples not present in mapping file"""
    sample_IDs = zip(*mapping[1:])[0]

    # remove unmapped ids from headers and coords
    for i in xrange(len(coords[0]) - 1, -1, -1):
        if not coords[0][i] in sample_IDs:
            del(coords[0][i])
            coords[1] = np.delete(coords[1], i, 0)

    # remove unmapped ids from edges
    if edges:
        for i in xrange(len(edges) - 1, -1, -1):
            edge = edges[i]
            if not edge[0] in sample_IDs or not edge[1] in sample_IDs:
                del(edges[i])


def make_3d_plots_invue(data, groups_and_colors, intp_pts, polyh_pts, offset):
    """Makes 3d plots given the groups_and_colors output.
    """
    data3d = data['coord'][1][:, :3]
    centroid = np.average(data3d, axis=0)
    for i in range(polyh_pts):
        idx = np.sqrt(np.sum(np.square(data3d - centroid), axis=1)).argmax()
        if i == 0:
            polypts = [data3d[idx]]
            centroid = polypts
        elif i == 1:
            polypts = np.append(polypts, [data3d[idx]], axis=0)
            centroid = np.average(polypts[:i + 1], axis=0)
        elif i == 2:
            polypts = np.append(polypts, [data3d[idx]], axis=0)
            polypts = np.append(polypts, [polypts[0]], axis=0)
            centroid = np.average(polypts[:i + 1], axis=0)
            lastpts = polypts[:-1]
        elif i == 3:
            polypts = np.append(polypts, [data3d[idx]], axis=0)
            polypts = np.append(polypts, polypts[i - 2:i], axis=0)
            lastpts = np.append(lastpts, [data3d[idx]], axis=0)
            centroid = np.average(lastpts[:-3], axis=0)
        else:
            polypts = np.append(polypts, [data3d[idx]], axis=0)
            polypts = np.append(polypts, polypts[i - 2:i - 1], axis=0)
            polypts = np.append(polypts, [data3d[idx]], axis=0)
            lastpts = np.append(lastpts, [data3d[idx]], axis=0)
            centroid = np.average(lastpts[:-3], axis=0)

    smp_lbl = {}
    smp_lbl_grp = {}
    for i in range(len(groups_and_colors)):
        labelname = groups_and_colors[i][0]
        groups = groups_and_colors[i][1]
        colors = groups_and_colors[i][2]
        data_colors = groups_and_colors[i][3]

        # Binning data per metadata info
        smp_lbl[labelname] = {'coords': [], 'headrs': []}
        smp_lbl_grp[labelname] = {}
        for gr in groups:
            smp_lbl_grp[labelname][gr] = {'coords': [], 'headrs': []}

            for elm in groups[gr]:
                try:
                    idx = data['coord'][0].index(elm)
                except (ValueError):
                    continue
                    #raise ValueError, 'ValueError: list.index(x): %s not in list' % elm

                # Creating interpolation points
                if intp_pts == 0:
                    smp_lbl_grp[labelname][gr]['coords'].append(
                        data['coord'][1][idx][:3])
                    smp_lbl_grp[labelname][gr]['headrs'].append(elm)
                else:
                    if len(smp_lbl_grp[labelname][gr]['coords']) == 0:
                        smp_lbl_grp[labelname][gr]['coords'].append(
                            data['coord'][1][idx][:3])
                        smp_lbl_grp[labelname][gr]['headrs'].append(elm)
                        pass
                    else:
                        new_pts = linear_gradient(
                            prev_pts,
                            data['coord'][1][idx][:3],
                            intp_pts + 2)
                        for j, tmp in enumerate(new_pts[1:]):
                            smp_lbl_grp[labelname][gr]['headrs'].append(
                                "%s.%d" %
                                (elm, j))
                            smp_lbl_grp[labelname][gr][
                                'coords'].append(array(tmp))

                prev_pts = data['coord'][1][idx][:3]

                # Saving the coords
                smp_lbl[
                    labelname][
                    'coords'].append(np.append(data['coord'][1][idx][:3],
                                               [data_colors[colors[gr]].toInt()]))
                smp_lbl[labelname]['headrs'].append(elm)

    return smp_lbl, smp_lbl_grp, polypts * offset


def generate_3d_plots_invue(
        prefs, data, dir_path, filename, intp_pts, polyh_pts, offset):
    """ Make files to be imported to inVUE
        http://sourceforge.net/projects/invue/"""

    # Validating existance of all columns
    for col in prefs:
        if col not in data['map'][0]:
            raise ValueError('Column given "%s" does not exits in mapping \
                file' % col)

    # Split matrix by labelname, groups & give colors
    groups_and_colors = iter_color_groups(data['map'], prefs)
    groups_and_colors = list(groups_and_colors)

    smp_lbl, smp_lbl_grp, polypts = make_3d_plots_invue(
        data, groups_and_colors,
        intp_pts, polyh_pts, offset)

    # Looping to binning result to write full and binned files
    for lbl in smp_lbl:
        for grp in smp_lbl_grp[lbl]:
            # writting individual files
            ind_path = "%s/%s_%s_%s.txt" % (dir_path, filename, lbl, grp)
            smp = smp_lbl_grp[lbl][grp]
            outfile = open(ind_path, 'w')
            outfile.write(
                format_coords(smp['headrs'],
                              smp['coords'],
                              [],
                              [],
                              False))
            outfile.close()
        # writing full file
        full_path = "%s/%s_%s.txt" % (dir_path, filename, lbl)
        outfile = open(full_path, 'w')
        outfile.write(
            format_coords(smp_lbl[lbl]['headrs'], smp_lbl[lbl]['coords'],
                          [], [], False))
        outfile.close()

    # Writing tetraVertices.txt
    ind_path = "%s/tetraVertices.txt" % (dir_path)
    outfile = open(ind_path, 'w')
    outfile.write('\n'.join(['\t'.join(map(str, row)) for row in polypts]))
    outfile.write('\n')
    outfile.close()


def generate_3d_plots(prefs, data, custom_axes, background_color, label_color,
                      dir_path='', data_file_path='', filename=None,
                      default_filename='out', ellipsoid_prefs=None,
                      user_supplied_edges=False, ball_scale=1.0,
                      arrow_colors={
                          'line_color': 'white',
                          'head_color': 'red'},
                      add_vectors=None, plot_scaled=False, plot_unscaled=True):
    """Make 3d plots according to coloring options in prefs."""

    if filename is None:
        filename = default_filename
    kinpath = os.path.join(data_file_path, filename) + ".kin"
    data_folder = os.path.split(data_file_path)[-1]
    kinlink = os.path.join('./', data_folder, filename) + ".kin"
    htmlpath = dir_path

    coord_header, coords, eigvals, pct_var, coords_low, coords_high = \
        data['coord']
    mapping = data['map']

    edges = None
    if 'edges' in data:
        edges = data['edges']

    taxa = None
    if 'taxa' in data:
        taxa = data['taxa']

    res = make_3d_plots(coord_header, coords, pct_var, mapping, prefs,
                        background_color, label_color,
                        taxa, custom_axes=custom_axes, edges=edges,
                        coords_low=coords_low, coords_high=coords_high,
                        ellipsoid_prefs=ellipsoid_prefs,
                        user_supplied_edges=user_supplied_edges,
                        ball_scale=ball_scale, arrow_colors=arrow_colors,
                        add_vectors=add_vectors, plot_scaled=plot_scaled,
                        plot_unscaled=plot_unscaled)

    # Validating if we should add specific values for the vectors option
    if add_vectors and add_vectors['vectors_algorithm'] and add_vectors['vectors'] and add_vectors['vectors_path']:

        # Two files are written one containing the ANOVA information for each
        # group and the other one containing the raw vectors
        f_vectors = open(
            os.path.join(htmlpath,
                         add_vectors['vectors_path']),
            'w')
        f_vectors_raw_values = open(
            os.path.join(htmlpath, 'vectors_raw_values.txt'), 'w')

        f_vectors.write(
            'Vectors algorithm: %s\n' %
            add_vectors['vectors_algorithm'])
        f_vectors_raw_values.write(
            'Vectors algorithm: %s\n' %
            add_vectors['vectors_algorithm'])

        if add_vectors['weight_by_vector']:
            f_vectors.write('** This output is weighted **\n')
            f_vectors_raw_values.write('** This output is weighted **\n')

        # Each group has different categories, output the
        # information per group and per category
        for group in add_vectors['vectors_output']:

            # buffer all the extra info per category of group (rms, avg, diff, etc)
            # see make_subgroup_vectors to see how this values are calculated
            per_category_information = []

            if len(add_vectors['vectors_output'][group].keys()) == 1:
                f_vectors.write('\nGrouped by %s: Only one value in the group.\n'
                                % (group))
            else:
                to_test = {}
                for cat in add_vectors['vectors_output'][group]:
                    if len(add_vectors['vectors_output'][group][cat]['vectors_vector']) != 0:
                        to_test[
                            cat] = add_vectors[
                            'vectors_output'][
                            group][
                            cat][
                            'vectors_vector']
                        per_category_information.append(
                            add_vectors[
                                'vectors_output'][
                                group][
                                cat][
                                'vectors_result'])
                    else:
                        add_vectors[
                            'vectors_output'][
                            group][
                            cat][
                            'vectors_result'] = nan

                # Not all dicts can be tested
                if can_run_ANOVA_trajectories(to_test):
                    labels, group_means, prob = run_ANOVA_trajetories(to_test)
                    if (len(labels) == len(group_means)):
                        f_vectors.write(
                            '\nGrouped by "%s", probability: %f\n' %
                            (group, prob))
                        f_vectors_raw_values.write(
                            '\nGrouped by "%s"\n' %
                            group)

                        # Each category of the current group has its own values
                        for i in range(len(labels)):

                            # print the per group header
                            f_vectors.write('For group: "{0}", the group means is: {1}\n'
                                            .format(labels[i], group_means[i]))
                            f_vectors_raw_values.write(
                                'For group: "%s"\n' %
                                labels[i])

                            # the corresponding message if there is one
                            if add_vectors['vectors_output'][group][cat]['vectors_message'] is not None:
                                f_vectors.write(
                                    '%s' %
                                    add_vectors[
                                        'vectors_output'][
                                        group][
                                        cat][
                                        'vectors_message'])
                                f_vectors_raw_values.write(
                                    '%s' %
                                    add_vectors[
                                        'vectors_output'][
                                        group][
                                        cat][
                                        'vectors_message'])

                            # print the extra information (diff, avg, rms ...)
                            f_vectors.write(
                                'The info is: {0}\n'.format(per_category_information[i]))
                            f_vectors_raw_values.write(
                                'The vector is:\n %s\n' %
                                to_test[labels[i]])

                    else:
                        f_vectors.write(
                            '\nGrouped by %s: this value can not be used\n' %
                            (group))
                else:
                    f_vectors.write(
                        '\nGrouped by %s: this value can not be used\n' %
                        (group))

        f_vectors.close()
        f_vectors_raw_values.close()

    # Write kinemage file
    f = open(kinpath, 'w')
    f.write('\n'.join(res))
    f.close()
    basename, extension = os.path.splitext(filename)
    filename = '%s_3D_PCoA_plots' % (basename)

    # Write html page with the kinemage embedded
    f2 = open(os.path.join(htmlpath, filename) + '.html', 'w')
    f2.write("<html><head></head><body><applet code='king/Kinglet.class' \
            archive='./jar/king.jar' width=800 height=600> \
            <param name='kinSource' value='%s'></body></html>" % (kinlink))
    f2.write('\n'.join(res))
    f2.close()


def can_run_ANOVA_trajectories(input_dict):
    """
        can_run_ANOVA_trajectories: checks if input_dict can be tested using
        ANOVA, without raising any RuntimeWarning. ANOVA testing requires
        for at least one element to have a size different to one. When a
        dictionary of size N by 1 is passed, a division by zero will happen.

        input_dict: dict with arrays as values for each of the keys
    """
    # Obtain the total number of elements in the dict
    total_values = sum([len(value) for value in input_dict.values()])
    return not (total_values == len(input_dict.keys()))


def avg_vector_for_group(
        group_name, ids, coord_dict, custom_axes, add_vectors):
    """
        insert_average_vectors: adds a new group that represents the average
        of the elements of the elements with the same first value (the added
        vector) in the coordinates matrix. A new group will be created for
        each of the available groups, named with the prefix avg.

        Input:
        group_name: name of the category that will be used to obtain the mean
        of the coordinates.

        ids: list of ids that belong to the group, the coordinates belonging to
        this list, are the ones to be used for the average.

        coord_dict: dictionary containing the coordinates values for each of
        the unique ids. New ids will be added with the prefix avg, for each
        of the elements of each group's average representation.

        custom_axes: added column to the coordinates matrix (usually a time
        vector). This first column is used to group together and average the
        coordinates, if none is provided all the coordinates of the group are
        averaged.

        add_vectors: reference dictionary used in other functions.

        Output:
        avg_coord_dict: dictionary of coordinates corresponding to the average
        of the COORD_DICT input argument. If there is an added vector, the
        average will be taken over that vector, else all the points will be
        averaged togther.

        avg_add_vectors: dictionary, modified properly so that the other
        functions can it as an input.
    """
    avg_add_vectors = add_vectors

    grouped_coords = {}

    # For each of the key-value pairs in the add_vectors dictionary, search
    # through each of the values, each of these values is a tuple containing
    # in the first element either a value for the custom axes or another value
    # to order them (if no custom axes provided), the second element is the id
    # of a list of coordinates, if this id is included in the input ids, then
    # add it to a group corresponding to the first element of this tuple
    for keys, values in add_vectors['vectors'].iteritems():
        # Be sure to sort the values otherwise, lines won't make much sense
        for value in signed_natsort(values):
            if value[1] in ids:
                # If the list hasn't been created yet, then create an empty one
                if value[0] not in grouped_coords.keys():
                    grouped_coords[value[0]] = []
                grouped_coords[value[0]].append(coord_dict[value[1]])

    avg_coord_dict = {}

    # Create a dictionary keyed by the custom axes (or the other value provided)
    # and valued with the means of all the lists of coordinates with the same
    # id
    for per_point_id, per_point_coords in grouped_coords.iteritems():
        avg_coord_dict[per_point_id] = mean(per_point_coords, axis=0)

    # Add vectors is used by other functions, so create a mock dict, be aware
    # of the order in the ids, otherwise the visualization will be confusing
    # the sorting of the keys must be casting the values to floating point type
    avg_add_vectors['vectors'] = {group_name: [(str(i), avg_id)
                                               for i, avg_id in enumerate(signed_natsort(avg_coord_dict.keys()))]}

    return avg_coord_dict, avg_add_vectors


def weight_by_vector(vector, w_vector):
    """
    weight_by_vector: weights the values of 'vector' given a weighting vector
    'w_vector'. Each value in 'vector' will be weighted by the 'rate of change'
    to 'optimal rate of change' ratio, meaning that when calling this function
    over evenly spaced 'w_vector' values, no change will be reflected on the
    output.

    Input:
    vector: numpy array of values to weight
    w_vector: numpy array used to weight 'vector'.

    Output:
    op_vector: numpy array representing a weighted version of 'vector'.
    """

    try:
        if len(vector) != len(w_vector):
            raise ValueError("vector (%d) & w_vector (%d) must be equal lengths"
                             % (len(vector), len(w_vector)))
    except TypeError:
        raise TypeError("vector and w_vector must be iterables")

    # check no repeated values are passed in the weighting vector
    if len(list(set(w_vector))) != len(w_vector):
        raise ValueError("The weighting vector must be a gradient")

    # no need to weight in case of a one element vector
    if len(w_vector) == 1:
        return vector

    op_vector = []

    # Cast to float so divisions have a floating point resolution
    total_length = float(max(w_vector) - min(w_vector))

    # reflects the expected gradient between subsequent values in w_vector
    # the first value isn't weighted so subtract one from the number of
    # elements
    optimal_gradient = total_length / (len(w_vector) - 1)

    for n, vector_value in enumerate(vector):
        # for all elements apply the weighting function
        if n != 0:
            op_vector.append(
                vector_value * (optimal_gradient) / numpy_abs((w_vector[n] - w_vector[n - 1])))
        # if it's the first element, just return it as is, no weighting to do
        else:
            op_vector.append(vector_value)

    return op_vector


def windowed_diff(vector, window_size):
    """
    windowed_diff: perform the first difference algorithm between windows of
    values in a vector and each value.

    Input:
    vector: numpy array of values to calculate the windowed_diff
    window_size: size of the window

    Output:
    op_vector: a vector where the Nth value is the difference between the mean
    of vector[N+1:N+1+window_size] and vector[N]. By definition this vector will
    have 'window_size' less elements than 'vector'.
    """

    # check for consistency in window size and vector size
    if window_size < 1 or not isinstance(window_size, (long, int)):
        raise ValueError("The window_size must be a positive integer")

    if len(vector) <= window_size:
        raise ValueError("The window_size must be smaller than the vector")

    # replicate the last element as many times as required
    for index in range(0, window_size):
        vector = append(vector, vector[-1:], axis=0)

    op_vector = []

    for index in range(0, len(vector) - window_size):
        # mean has to be over axis 0 so it handles vectors of vectors
        element = mean(vector[(index + 1):(index + 1 + window_size)], axis=0)
        op_vector.append(element - vector[index])

    return op_vector
