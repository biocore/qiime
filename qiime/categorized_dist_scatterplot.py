#!/usr/bin/env python
# File created on 19 Jan 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


from qiime.parse import parse_metadata_state_descriptions
from qiime.filter import get_sample_ids
import numpy


def get_avg_dists(state1_samids, state2_samids, distdict):
    """ foreach sample in stat1_sams, return average dist to all state2_sams

    doesn't include distance to self in average, but usually state2 doesn't
    overlap with state1, so it dist to self doesn't apply
    returns a list of length = len(state1_samids)
    """
    # go through dmtx
    state1_avg_dists = []
    for sam1 in state1_samids:
        dists = []
        for sam2 in state2_samids:
            if sam1 == sam2:
                continue
            dists.append(distdict[sam1][sam2])
        state1_avg_dists.append(numpy.mean(dists))
    return state1_avg_dists


def get_sam_ids(map_data, map_header, colorby, cat,
                primary_state, secondary_state):
    """ returns all sample ids matching the state strings and colorby:cat

    colorby: eg: 'Country', or pass None to not filter only colorby:cat samples
    cat: e.g.: 'USA'
    primary_state: e.g.: 'AgeCategory:Child'
    secondary state can be None, or like primary state

    returns uniquified lists in randomized order
    """
    if colorby is None:
        sample_ids = [sam[0] for sam in map_data]
    else:

        sample_ids = get_sample_ids(
            map_data, map_header, {colorby: [cat]})
    # primary key is the category label, e.g. AgeCategory
    # value is the val for that category, e.g. Adult

    # go through age1/age2
    primary_states = parse_metadata_state_descriptions(primary_state)
    if colorby is not None:
        primary_states[colorby] = [cat]
    state1_samids = get_sample_ids(
        map_data, map_header, primary_states)

    if secondary_state is None:
        state2_samids = set(sample_ids).difference(set(state1_samids))
    else:
        secondary_states =\
            parse_metadata_state_descriptions(secondary_state)
        if colorby is not None:
            secondary_states[colorby] = [cat]
        state2_samids = get_sample_ids(
            map_data, map_header, secondary_states)

    return list(set(state1_samids)), list(set(state2_samids))
