#!/usr/bin/env python
# make_otu_table: makes sample x OTU table
__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski", "Adam Robbins-Pianka",
               "Sami Pietila"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

"""Makes sample x OTU table from OTU map and taxonomy.

Assumes that in the OTU map, the ids are in the format lib_seq, e.g.
M3FclSwb_1023. Will not work if this assumption is not met. Splits on last
underscore only so should be relatively robust to underscore in sample id.
"""

from datetime import datetime
from collections import defaultdict
from string import strip
from sys import stderr

from numpy import array, zeros
from cogent.util.misc import flatten
from biom.table import Table

from qiime.parse import parse_otu_map
from qiime.util import get_generated_by_for_biom_tables


def libs_from_seqids(seq_ids, delim='_'):
    """Returns set of libraries."""
    all_libs = set([i.rsplit(delim, 1)[0] for i in seq_ids])
    return all_libs


def seqids_from_otu_to_seqid(otu_to_seqid):
    """Returns set of all seq ids from libs"""
    return set(flatten(otu_to_seqid.values()))


def make_otu_table(otu_map_f, otu_to_taxonomy=None, delim='_', table_id=None,
                   otu_ids_to_exclude=None, sample_metadata=None):
    """Generate a BIOM table from an OTU map

    Parameters
    ----------
    otu_map_f : file-like object
        The OTU map. Jagged tab-separated file where the first column contains
        the OTU ID and subsequent columns contain sequence IDs belonging to
        that OTU
    otu_to_taxonomy : dict, optional
        Defaults to ``None``. If supplied, the dict maps OTU IDs to taxonomies
    delim : str, optional
        Defaults to "_". The delimiter that is used in the sequence IDs to join
        the sample ID to the sequence number
    table_id : object, optional
        Defaults to ``None``. The identifier that will be given to the
        generated BIOM table
    otu_ids_to_exclude : iterable, optional
        Defaults to ``None``. If present, these OTUs will not be added to the
        OTU table from the OTU map
    sample_metadata : dict of dicts, optional
        Defaults to ``None``. If supplied, keys in the outer dict should be
        sample IDs, and keys in the inner dicts should be column names.
    """
    data, sample_ids, otu_ids = parse_otu_map(
        otu_map_f, delim=delim, otu_ids_to_exclude=otu_ids_to_exclude)

    if otu_to_taxonomy is not None:
        otu_metadata = []
        for o in otu_ids:
            otu_metadata.append({'taxonomy': otu_to_taxonomy.get(o, ["None"])})
    else:
        otu_metadata = None

    # if sample_metadata is supplied, put in index-order with the OTU map's
    # sample_ids, and do not include samples that were in the mapping file
    # but NOT in the OTU map
    if sample_metadata is not None:
        try:
            sample_metadata = [sample_metadata[sample_id]
                               for sample_id in sample_ids]
        except KeyError:
            raise KeyError("Sample IDs found in OTU map without sample "
                           "metadata")

    try:
        return Table(data, otu_ids, sample_ids,
                     observation_metadata=otu_metadata, 
                     sample_metadata=sample_metadata, table_id=table_id,
                     generated_by=get_generated_by_for_biom_tables(),
                     create_date=datetime.now().isoformat())
    except ValueError as e:
        raise ValueError("Couldn't create OTU table. Is your OTU map empty?"
                         " Original error message: %s" % (str(e)))
