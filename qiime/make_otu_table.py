#!/usr/bin/env python
# make_otu_table: makes sample x OTU table
__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

"""Makes sample x OTU table from OTU map and taxonomy.

Assumes that in the OTU map, the ids are in the format lib_seq, e.g.
M3FclSwb_1023. Will not work if this assumption is not met. Splits on last
underscore only so should be relatively robust to underscore in sample id.
"""

from collections import defaultdict
from string import strip
from sys import stderr
from numpy import array, zeros
from cogent.util.misc import flatten
from qiime.format import format_otu_table
from qiime.parse import parse_otu_map
from qiime.format import format_biom_table
from biom.table import Table, table_factory


def libs_from_seqids(seq_ids, delim='_'):
    """Returns set of libraries."""
    all_libs = set([i.rsplit(delim, 1)[0] for i in seq_ids])
    return all_libs


def seqids_from_otu_to_seqid(otu_to_seqid):
    """Returns set of all seq ids from libs"""
    return set(flatten(otu_to_seqid.values()))


def make_otu_table(otu_map_f,
                   otu_to_taxonomy=None,
                   delim='_',
                   table_id=None,
                   sample_metadata=None,
                   constructor=Table):
                   input_is_dense=False):

    data, sample_ids, otu_ids = parse_otu_map(otu_map_f, delim)

    if otu_to_taxonomy is not None:
        otu_metadata = []
        for o in otu_ids:
            try:
                otu_metadata.append({'taxonomy': otu_to_taxonomy[o]})
            except KeyError:
                otu_metadata.append({'taxonomy': ["None"]})
    else:
        otu_metadata = None

    if sample_metadata is not None:
        raise NotImplementedError(
            "Passing of sample metadata to make_otu_table is not currently supported.")
    try:
        otu_table = table_factory(data, sample_ids, otu_ids,
                                  sample_metadata=sample_metadata,
                                  observation_metadata=otu_metadata,
                                  table_id=table_id,
                                  constructor=constructor,
                                  dtype=int, input_is_dense=input_is_dense)
    except ValueError as e:
        raise ValueError("Couldn't create OTU table. Is your OTU map empty?"
                         " Original error message: %s" % (str(e)))
    return format_biom_table(otu_table)
