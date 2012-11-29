#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Contains functions used in the make_distance_boxplots.py script."""

from collections import defaultdict
from qiime.colors import data_colors, data_color_order
from qiime.make_distance_histograms import matplotlib_rgb_color
from qiime.util import MetadataMap

def color_field_states(mapping_f, samp_ids, field, field_states,
                       color_by_field):
    """Colors one field by another.

    Returns a list of matplotlib-compatible colors, one for each of the input
    field_states. Also returns a dictionary mapping color_by_field states to
    colors (useful for building a legend, for example).

    If there are not enough colors available (they are drawn from
    qiime.colors.data_colors), an error will be raised as the color mapping
    (and legend) will be ambiguous.

    A one-to-one mapping must exist between each field_state and its
    corresponding color_by field state (otherwise it is unclear which
    corresponding color_by field state should be used to color it by). An error
    will be raised if this one-to-one mapping does not exist.

    Arguments:
        mapping_f - the mapping file (file-like object)
        samp_ids - a list of sample IDs to consider in the mapping file. Only
            these sample IDs will be used when coloring field states
        field - the field in the mapping file to color
        field_states - the field states in field to color
        color_by_field - the field in the mapping file to color field_states by
    """
    colors = []
    color_pool = [matplotlib_rgb_color(data_colors[color].toRGB())
                  for color in data_color_order]
    metadata_map = MetadataMap.parseMetadataMap(mapping_f)

    for field_to_check in field, color_by_field:
        if field_to_check not in metadata_map.CategoryNames:
            raise ValueError("The field '%s' is not in the metadata mapping "
                             "file's column headers." % field_to_check)

    field_mapping = defaultdict(list)
    all_field_states = metadata_map.getCategoryValues(samp_ids, field)
    all_color_by_states = metadata_map.getCategoryValues(samp_ids,
                                                         color_by_field)

    if len(set(field_states) - set(all_field_states)) != 0:
        raise ValueError("Encountered unrecognizable field state(s) in %r "
                         "for field '%s'." % (field_states, field))

    for field_state, color_by_state in zip(all_field_states,
                                           all_color_by_states):
        if field_state in field_states:
            field_mapping[field_state].append(color_by_state)

    color_mapping = {}
    for field_state in field_states:
        color_by_states = set(field_mapping[field_state])

        if len(color_by_states) != 1:
            raise ValueError("The field '%s' to color by does not have a "
                             "one-to-one mapping with field '%s'. Coloring "
                             "would be ambiguous." % (color_by_field, field))

        color_by_state = list(color_by_states)[0]
        if color_by_state not in color_mapping:
            if len(color_pool) > 0:
                color_mapping[color_by_state] = color_pool.pop(0)
            else:
                raise ValueError("There are not enough available QIIME colors "
                                 "to color each of the field states in field "
                                 "'%s'. Coloring would be ambiguous." %
                                 color_by_field)

        colors.append(color_mapping[color_by_state])

    return colors, color_mapping
