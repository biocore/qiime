#!/usr/bin/env python
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Justin Kuczynski",
               "Julia Goodrich", "Antonio Gonzalez Pena",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Development"

"""Contains code for summarizing OTU table based on metadata."""

from collections import defaultdict
from string import strip
from biom.table import SparseTaxonTable

ONE_TO_MANY_TYPES = ['first', 'add', 'divide']

def make_summary(otu_table, level, upper_percentage=0.0, lower_percentage=1.0,
                 md_as_string=False, md_identifier='taxonomy', delimiter=';',
                 one_to_many='first', absolute_abundance=False,
                 constructor=SparseTaxonTable):
    """Returns a taxa summary table in BIOM format.

    ``otu_table`` should be a BIOM table with absolute or relative abundances.

    WARNING: this function will work with relative abundance BIOM tables, but
    there is a special case where you may end up with a summarized BIOM table
    that is neither in relative or absolute abundances. This can happen if you
    supply a relative abundance table and one-to-many relationships exist in
    the metadata, and you specify ``one_to_many='add'``. With this type of
    one-to-many relationship handling, the counts in the table may increase, so
    the resulting summarized table will need to have its relative abundances
    recalculated (which will happen by default with
    ``absolute_abudance=False``). Thus, this problem only exists if you specify
    ``absolute_abundance=True``.

    WARNING: Specifying ``upper_percentage`` and ``lower_percentage`` may
    change what summarized taxa are filtered out depending on whether
    ``otu_table`` has absolute or relative abundances. In most cases, users
    will likely want to supply ``otu_table`` in absolute abundances.

    ``one_to_many='divide'`` and ``absolute_abundance=True`` are not supported
    because counts will not always evenly divide.
    """
    if otu_table.ObservationMetadata is None:
        raise ValueError("BIOM table does not contain any observation "
                         "metadata (e.g., taxonomy). You can add metadata to "
                         "it using add_metadata.py.")

    collapse_fn = _make_collapse_fn(level, md_identifier, md_as_string,
                                    delimiter=delimiter,
                                    one_to_many=one_to_many,
                                    absolute_abundance=absolute_abundance)

    if one_to_many in ['add', 'divide']:
        one_to_many_mode = one_to_many
    else:
        # If we're not adding or dividing, it doesn't matter which mode we
        # pick, so we'll choose to add since it should have slightly better
        # performance than divide mode.
        one_to_many_mode = 'add'

    ts_table = otu_table.collapseObservationsByMetadata(collapse_fn,
            norm=False, min_group_size=1, include_collapsed_metadata=False,
            constructor=constructor, one_to_many=True,
            one_to_many_mode=one_to_many_mode, strict=True)

    filter_fn = _make_abundance_filter_fn(ts_table, upper_percentage,
                                          lower_percentage)
    ts_table = ts_table.filterObservations(filter_fn)
    ts_table = ts_table.sortByObservationId(sorted)

    if not absolute_abundance:
        ts_table = ts_table.normObservationBySample()

    return ts_table

def add_summary_mapping(otu_table,
                        mapping, 
                        level,
                        md_as_string=False,
                        md_identifier='taxonomy',
                        delimiter=';',
                        one_to_many='first',
                        absolute_abundance=False):
    """Returns sample summary of sample counts by taxon
    
    Summary is keyed by sample_id, valued by otu counts for each taxon
    Taxon order is a list of taxons where idx n corresponds to otu count idx n
    """
    ts_table = make_summary(otu_table, level, md_as_string=md_as_string,
                            md_identifier=md_identifier, delimiter=delimiter,
                            one_to_many=one_to_many,
                            absolute_abundance=absolute_abundance)

    summary = defaultdict(list)
    for row in mapping:
        # grab sample idx if the sample exists, otherwise ignore it
        sample_id = row[0]
        if not ts_table.sampleExists(sample_id):
            continue
        sample_idx = ts_table.getSampleIndex(sample_id)

        for obs_v in ts_table.iterObservationData():
            summary[sample_id].append(obs_v[sample_idx])

    return summary, ts_table.ObservationIds

# Taken from PICRUSt's categorize_by_function.py script
# (http://picrust.github.io/picrust/) and modified.
def _make_collapse_fn(level, md_identifier='taxonomy', md_as_string=False,
                      missing_name='Other', delimiter=';',
                      one_to_many='first', absolute_abundance=False):
    """Returns a collapsing function for 1-1 and 1-M relationships."""
    if one_to_many not in ONE_TO_MANY_TYPES:
        raise ValueError('Encountered unrecognized method "%s" for handling '
                         'one-to-many relationships in metadata.' %
                         one_to_many)

    if md_as_string:
        # Strings will always be processed as one-to-one because they will be
        # split into a single-level list.
        one_to_many = 'first'
        def process_md(v):
            return v.split(delimiter)
    else:
        def process_md(v):
            return v

    def collapse(md):
        if md_identifier not in md:
            raise KeyError("An observation (e.g., OTU) does not have the "
                           "metadata identifier '%s'. Did you pass the "
                           "correct metadata identifier?" % md_identifier)

        md_val = process_md(md[md_identifier])

        is_single_level = False
        for md_item in md_val:
            if isinstance(md_item, basestring):
                # If we have a list of strings, we want the whole thing (only).
                md_item = md_val
                is_single_level = True

            num_ranks = len(md_item)
            if num_ranks > level:
                md_item = md_item[:level]
            elif num_ranks < level:
                md_item.extend([missing_name
                                for i in range(level - num_ranks)])
            else:
                # md_item is the correct number of levels.
                pass

            md_item = delimiter.join(md_item)

            yield md_item, md_item

            # If we only have one list of strings or we were told to only use
            # the first metadata item, we're done.
            if is_single_level or one_to_many == 'first':
                break
            elif one_to_many == 'divide' and absolute_abundance:
                raise ValueError("Encountered a one-to-many relationship "
                                 "between an observation (e.g., OTU) and its "
                                 "metadata. The one-to-many mode '%s' cannot "
                                 "be used when absolute abundance is "
                                 "specified as the output type because "
                                 "observation counts will not always divide "
                                 "evenly between the multiple metadata "
                                 "values. The output type must be in relative "
                                 "abundances." % one_to_many)

    return collapse

def _make_abundance_filter_fn(table, upper_percentage, lower_percentage):
    """Returns a filtering function for abundance-based filtering."""
    total = table.sum()
    lower = lower_percentage * total
    upper = upper_percentage * total

    def filter_fn(obs_val, obs_id, obs_md):
        if obs_val.sum() > lower or obs_val.sum() < upper:
            keep_obs = False
        else:
            keep_obs = True

        return keep_obs

    return filter_fn
