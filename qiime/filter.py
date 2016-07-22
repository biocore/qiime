#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Will Van Treuren", "Daniel McDonald",
               "Jai Ram Rideout", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from collections import defaultdict
from random import shuffle, sample
from numpy import array, inf

from skbio.parse.sequences import parse_fasta, parse_fastq
from skbio.format.sequences import format_fastq_record
from biom import load_table

from qiime.parse import (parse_distmat, parse_mapping_file,
                         parse_metadata_state_descriptions)
from qiime.format import format_distance_matrix, format_mapping_file
from qiime.util import MetadataMap


def get_otu_ids_from_taxonomy_f(positive_taxa=None,
                                negative_taxa=None,
                                metadata_field="taxonomy"):
    """ return function to pass to Table.filter_observations for taxon-based filtering
        positive_taxa : a list of strings that will be compared to each
         taxonomy level in an observation's (i.e., OTU's) metadata_field. If
         one of the levels matches exactly (except for case) to an item in
         positive_taxa, that OTU will be marked for retention. Default: All
         OTUs are retained.
        negative_taxa : a list of strings that will be compared to each
         taxonomy level in an observation's (i.e., OTU's) metadata_field. If
         one of the levels matches exactly (except for case)  to an item in
         negative_taxa, that OTU will be marked for removal. Default: All
         OTUs are retained.
        metadata_field : the metadata field to look up in the
         observation metadata

        Note: string matches are case insensitive.
    """
    # define a positive screening function - if the user doesn't pass
    # positive_taxa, all OTUs will pass this filter
    # (i.e., be marked for retention)
    if positive_taxa is None:
        positive_taxa = set()

        def positive_screen(e):
            return True
    else:
        positive_taxa = set([t.strip().lower() for t in positive_taxa])

        def positive_screen(e):
            return e in positive_taxa

    # define a negative screening function - if the user doesn't pass
    # negative_taxa, all OTUs will pass this filter
    # (i.e., be marked for retention)
    if negative_taxa is None:
        negative_taxa = set()

        def negative_screen(e):
            return False
    else:
        negative_taxa = set([t.strip().lower() for t in negative_taxa])

        def negative_screen(e):
            return e in negative_taxa

    # The positive_taxa and negative_taxa lists must be mutually exclusive.
    if len(positive_taxa & negative_taxa) != 0:
        raise ValueError("Your positive and negative taxa lists contain "
                         "overlapping values. These lists must be mutually "
                         "exclusive.\nOffending values are: %s" %
                         ' '.join(positive_taxa & negative_taxa))

    # Define the function that can be passed to Table.filter_observations
    def result(v, oid, md):
        positive_hit = False
        negative_hit = False
        for e in md[metadata_field]:
            if positive_screen(e.strip().lower()):
                # Note that we don't want to just do
                # positive_hit = positive_screen(e.strip())
                # we're checking whether any e hits the positive taxa
                # and doing that be the same as
                # positive_hit = md[metadata_field][-1]
                positive_hit = True
            if negative_screen(e.strip().lower()):
                # see note in previous if statement for why we don't use
                # negative_hit = negative_screen(e.strip())
                negative_hit = True
        return positive_hit and not negative_hit

    return result


def sample_ids_from_metadata_description(mapping_f, valid_states_str):
    """ Given a description of metadata, return the corresponding sample ids
    """
    map_data, map_header, map_comments = parse_mapping_file(mapping_f)
    valid_states = parse_metadata_state_descriptions(valid_states_str)
    sample_ids = get_sample_ids(map_data, map_header, valid_states)

    if len(sample_ids) < 1:
        raise ValueError("All samples have been filtered out for the criteria"
                         " described in the valid states")

    return sample_ids


def get_sample_ids(map_data, map_header, states):
    """Takes col states in {col:[vals]} format.

    If val starts with !, exclude rather than include.

    Combines cols with and, states with or.

    For example, Study:Dog,Hand will return rows where Study is Dog or Hand;
    Study:Dog,Hand;BodySite:Palm,Stool will return rows where Study is Dog
    or Hand _and_ BodySite is Palm or Stool; Study:*,!Dog;BodySite:*,!Stool
    will return all rows except the ones where the Study is Dog or the BodySite
    is Stool.
    """

    name_to_col = dict([(s, map_header.index(s)) for s in states])
    good_ids = []
    for row in map_data:  # remember to exclude header
        include = True
        for s, vals in states.items():
            curr_state = row[name_to_col[s]]
            include = include and (curr_state in vals or '*' in vals) \
                and not '!' + curr_state in vals
        if include:
            good_ids.append(row[0])
    return good_ids


def sample_ids_from_category_state_coverage(mapping_f,
                                            coverage_category,
                                            subject_category,
                                            min_num_states=None,
                                            required_states=None,
                                            considered_states=None,
                                            splitter_category=None):
    """Filter sample IDs based on subject's coverage of a category.

    Given a category that groups samples by subject (subject_category), samples
    are filtered by how well a subject covers (i.e. has at least one sample
    for) the category states in coverage_category.

    Two filtering criteria are provided (min_num_states and required_states).
    At least one must be provided. If both are provided, the subject must meet
    both criteria to pass the filter (i.e. providing both filters is an AND,
    not an OR, operation).

    A common use case is to provide a 'time' category for coverage_category and
    an 'individual' category for subject_category in order to filter out
    individuals from a study that do not have samples for some minimum number
    of timepoints (min_num_states) and that do not have samples for certain
    timepoints (required_states). For example, this could be the first and last
    timepoints in the study.

    Returns a set of sample IDs to keep, the number of subjects that were
    kept, and a set of the unique category states in coverage_category that
    were kept. The set of sample IDs is not guaranteed to be in any specific
    order relative to the order of sample IDs or subjects in the mapping file.

    Arguments:
        mapping_f - metadata mapping file (file-like object)
        coverage_category - category to test subjects' coverage (string)
        subject_category - category to group samples by subject (string)
        min_num_states - minimum number of category states in coverage_category
            that a subject must cover (i.e. have at least one sample for) to be
            included in results (integer)
        required_states - category states in coverage_category that must be
            covered by a subject's samples in order to be included in results
            (list of strings or items that can be converted to strings)
        considered_states - category states that are counted toward the
            min_num_states (list of strings or items that can be converted to
            strings)
        splitter_category - category to split input mapping file on prior to
            processing. If not supplied, the mapping file will not be split. If
            supplied, a dictionary mapping splitter_category state to results
            will be returned instead of the three-element tuple. The supplied
            filtering criteria will apply to each split piece of the mapping
            file independently (e.g. if an individual passes the filters for
            the tongue samples, his/her tongue samples will be included for
            the tongue results, even if he/she doesn't pass the filters for the
            palm samples)
    """
    metadata_map = MetadataMap.parseMetadataMap(mapping_f)

    # Make sure our input looks sane.
    categories_to_test = [coverage_category, subject_category]
    if splitter_category is not None:
        categories_to_test.append(splitter_category)

    if 'SampleID' in categories_to_test:
        raise ValueError("The 'SampleID' category is not suitable for use in "
                         "this function. Please choose a different category "
                         "from the metadata mapping file.")

    for category in categories_to_test:
        if category not in metadata_map.CategoryNames:
            raise ValueError("The category '%s' is not in the metadata "
                             "mapping file." % category)

    if len(set(categories_to_test)) < len(categories_to_test):
        raise ValueError("The coverage, subject, and (optional) splitter "
                         "categories must all be unique.")

    if required_states is not None:
        # required_states must be in coverage_category's states in the mapping
        # file.
        required_states = set(map(str, required_states))
        valid_coverage_states = set(metadata_map.getCategoryValues(
            metadata_map.sample_ids, coverage_category))
        invalid_coverage_states = required_states - valid_coverage_states

        if invalid_coverage_states:
            raise ValueError("The category state(s) '%s' are not in the '%s' "
                             "category in the metadata mapping file." %
                             (', '.join(invalid_coverage_states),
                              coverage_category))

    if considered_states is not None:
        # considered_states is not as restrictive as required_states - we don't
        # require that these are present, so it's OK if some of the states
        # listed here don't actually show up in the mapping file (allowing
        # the user to pass something like range(100) to consider only states
        # that fall in some range)
        considered_states = set(map(str, considered_states))
        # define a function to determine if a state should be considered
        consider_state = lambda s: s in considered_states
    else:
        # define a dummy function to consider all states (the default
        # if the user does not provide a list of considered_states)
        consider_state = lambda s: True

    if min_num_states is None and required_states is None:
        raise ValueError("You must specify either the minimum number of "
                         "category states the subject must have samples for "
                         "(min_num_states), or the minimal category states "
                         "the subject must have samples for "
                         "(required_states), or both. Supplying neither "
                         "filtering criteria is not supported.")

    if splitter_category is None:
        results = _filter_sample_ids_from_category_state_coverage(
            metadata_map, metadata_map.sample_ids, coverage_category,
            subject_category, consider_state, min_num_states,
            required_states)
    else:
        # "Split" the metadata mapping file by extracting only sample IDs that
        # match the current splitter category state and using those for the
        # actual filtering.
        splitter_category_states = defaultdict(list)
        for samp_id in metadata_map.sample_ids:
            splitter_category_state = \
                metadata_map.getCategoryValue(samp_id, splitter_category)
            splitter_category_states[splitter_category_state].append(samp_id)

        results = {}
        for splitter_category_state, sample_ids in \
                splitter_category_states.items():
            results[splitter_category_state] = \
                _filter_sample_ids_from_category_state_coverage(
                    metadata_map, sample_ids, coverage_category,
                    subject_category, consider_state, min_num_states,
                    required_states)

    return results


def _filter_sample_ids_from_category_state_coverage(metadata_map,
                                                    sample_ids,
                                                    coverage_category,
                                                    subject_category,
                                                    consider_state_fn,
                                                    min_num_states=None,
                                                    required_states=None):
    """Helper function to perform filtering based on category state coverage.

    Not explicitly unit-tested because it is implicitly tested by
    sample_ids_from_category_state_coverage's unit tests.
    """
    # Build mapping from subject to sample IDs.
    subjects = defaultdict(list)
    for samp_id in sample_ids:
        subject = metadata_map.getCategoryValue(samp_id, subject_category)
        subjects[subject].append(samp_id)

    # Perform filtering.
    samp_ids_to_keep = []
    num_subjects_kept = 0
    states_kept = []
    for subject, samp_ids in subjects.items():
        subject_covered_states = set(
            metadata_map.getCategoryValues(samp_ids, coverage_category))

        # Short-circuit evaluation of ANDing filters.
        keep_subject = True
        if min_num_states is not None:
            # note: when summing a list of boolean values, True == 1 and
            # False == 0
            if sum([consider_state_fn(s) for s in subject_covered_states]) < \
               min_num_states:
                keep_subject = False
        if keep_subject and required_states is not None:
            if len(subject_covered_states & required_states) != \
               len(required_states):
                keep_subject = False

        if keep_subject:
            samp_ids_to_keep.extend(samp_ids)
            states_kept.extend(subject_covered_states)
            num_subjects_kept += 1

    return set(samp_ids_to_keep), num_subjects_kept, set(states_kept)


def filter_fasta(input_seqs_f, output_seqs_f, seqs_to_keep, negate=False,
                 seqid_f=None):
    """ Write filtered input_seqs to output_seqs_f which contains only seqs_to_keep

        input_seqs can be the output of parse_fasta or parse_fastq
    """
    if seqid_f is None:
        seqs_to_keep_lookup = {}.fromkeys([seq_id.split()[0]
                                           for seq_id in seqs_to_keep])

        # Define a function based on the value of negate
        if not negate:
            def keep_seq(seq_id):
                return seq_id.split()[0] in seqs_to_keep_lookup
        else:
            def keep_seq(seq_id):
                return seq_id.split()[0] not in seqs_to_keep_lookup

    else:
        if not negate:
            keep_seq = seqid_f
        else:
            keep_seq = lambda x: not seqid_f(x)

    for seq_id, seq in parse_fasta(input_seqs_f):
        if keep_seq(seq_id):
            output_seqs_f.write('>%s\n%s\n' % (seq_id, seq))
    output_seqs_f.close()


def filter_fastq(input_seqs_f, output_seqs_f, seqs_to_keep, negate=False,
                 seqid_f=None):
    """ Write filtered input_seqs to output_seqs_f which contains only seqs_to_keep

        input_seqs can be the output of parse_fasta or parse_fastq
    """
    if seqid_f is None:
        seqs_to_keep_lookup = {}.fromkeys([seq_id.split()[0]
                                           for seq_id in seqs_to_keep])
        # Define a function based on the value of negate
        if not negate:
            def keep_seq(seq_id):
                return seq_id.split()[0] in seqs_to_keep_lookup
        else:
            def keep_seq(seq_id):
                return seq_id.split()[0] not in seqs_to_keep_lookup

    else:
        if not negate:
            keep_seq = seqid_f
        else:
            keep_seq = lambda x: not seqid_f(x)

    for seq_id, seq, qual in parse_fastq(input_seqs_f,
                                         enforce_qual_range=False):
        if keep_seq(seq_id):
            output_seqs_f.write(format_fastq_record(seq_id, seq, qual))
    output_seqs_f.close()


def filter_mapping_file(map_data, map_header, good_sample_ids,
                        include_repeat_cols=False, column_rename_ids=None):
    """Filters map according to several criteria.

    - keep only sample ids in good_sample_ids
    - drop cols that are different in every sample (except id)
    - drop cols that are the same in every sample
    """
    # keeping samples
    to_keep = []
    to_keep.extend([i for i in map_data if i[0] in good_sample_ids])

    # keeping columns
    headers = []
    to_keep = zip(*to_keep)
    headers.append(map_header[0])
    result = [to_keep[0]]

    if column_rename_ids:
        # reduce in 1 as we are not using the first colum (SampleID)
        column_rename_ids = column_rename_ids - 1
        for i, l in enumerate(to_keep[1:-1]):
            if i == column_rename_ids:
                if len(set(l)) != len(result[0]):
                    raise ValueError(
                        "The column to rename the samples is not unique.")
                result.append(result[0])
                result[0] = l
                headers.append('SampleID_was_' + map_header[i + 1])
            elif include_repeat_cols or len(set(l)) > 1:
                headers.append(map_header[i + 1])
                result.append(l)
    else:
        for i, l in enumerate(to_keep[1:-1]):
            if include_repeat_cols or len(set(l)) > 1:
                headers.append(map_header[i + 1])
                result.append(l)
    headers.append(map_header[-1])
    result.append(to_keep[-1])

    result = map(list, zip(*result))

    return headers, result


def filter_mapping_file_from_mapping_f(
        mapping_f, sample_ids_to_keep, negate=False):
    """ Filter rows from a metadata mapping file """
    mapping_data, header, comments = parse_mapping_file(mapping_f)
    filtered_mapping_data = []
    sample_ids_to_keep = {}.fromkeys(sample_ids_to_keep)

    for mapping_datum in mapping_data:
        hit = mapping_datum[0] in sample_ids_to_keep
        if hit and not negate:
            filtered_mapping_data.append(mapping_datum)
        elif not hit and negate:
            filtered_mapping_data.append(mapping_datum)
        else:
            pass
    return format_mapping_file(header, filtered_mapping_data)


def filter_mapping_file_by_metadata_states(mapping_f, valid_states_str):
    sample_ids_to_keep = sample_ids_from_metadata_description(
        mapping_f,
        valid_states_str)
    mapping_f.seek(0)
    return filter_mapping_file_from_mapping_f(mapping_f, sample_ids_to_keep)


def filter_samples_from_distance_matrix(dm, samples_to_discard, negate=False):
    """ Remove specified samples from distance matrix

        dm: (sample_ids, dm_data) tuple, as returned from
         qiime.parse.parse_distmat; or a file handle that can be passed
         to qiime.parse.parse_distmat

    """
    try:
        sample_ids, dm_data = dm
    except ValueError:
        # input was provide as a file handle
        sample_ids, dm_data = parse_distmat(dm)

    sample_lookup = {}.fromkeys([e.split()[0] for e in samples_to_discard])
    temp_dm_data = []
    new_dm_data = []
    new_sample_ids = []

    if negate:
        def keep_sample(s):
            return s in sample_lookup
    else:
        def keep_sample(s):
            return s not in sample_lookup

    for row, sample_id in zip(dm_data, sample_ids):
        if keep_sample(sample_id):
            temp_dm_data.append(row)
            new_sample_ids.append(sample_id)
    temp_dm_data = array(temp_dm_data).transpose()

    for col, sample_id in zip(temp_dm_data, sample_ids):
        if keep_sample(sample_id):
            new_dm_data.append(col)
    new_dm_data = array(new_dm_data).transpose()

    return format_distance_matrix(new_sample_ids, new_dm_data)


def negate_tips_to_keep(tips_to_keep, tree):
    """ Return the list of tips in the tree that are not in tips_to_keep"""
    tips_to_keep = set(tips_to_keep)
    # trees can return node names in ways that have multiple quotes, e.g.
    # '"node_1"' or ''node_1''. remove them or it can cause problems with
    # tips_to_keep not matching
    tmp_tips = set([tip.Name for tip in tree.tips()])
    tips = set([t.strip('\'').strip('\"') for t in tmp_tips])
    return tips - tips_to_keep


def get_seqs_to_keep_lookup_from_biom(biom_f):
    otu_table = load_table(biom_f)
    return set(otu_table.ids(axis='observation'))


def get_seqs_to_keep_lookup_from_seq_id_file(id_to_keep_f):
    """generate a lookup dict of chimeras in chimera file."""
    return (
        set([l.split()[0].strip()
            for l in id_to_keep_f if l.strip() and not l.startswith('#')])
    )
get_seq_ids_from_seq_id_file = get_seqs_to_keep_lookup_from_seq_id_file


def get_seqs_to_keep_lookup_from_fasta_file(fasta_f):
    """return the sequence ids within the fasta file"""
    return (
        set([seq_id.split()[0] for seq_id, seq in parse_fasta(fasta_f)])
    )
get_seq_ids_from_fasta_file = get_seqs_to_keep_lookup_from_fasta_file

# start functions used by filter_samples_from_otu_table.py and
# filter_otus_from_otu_table.py


def get_filter_function(ids_to_keep, min_count, max_count,
                        min_nonzero, max_nonzero, negate_ids_to_keep=False):
    if negate_ids_to_keep:
        def f(data_vector, id_, metadata):
            return (id_ not in ids_to_keep) and \
                   (min_count <= data_vector.sum() <= max_count) and \
                   (min_nonzero <= (data_vector > 0).sum() <= max_nonzero)
    else:
        def f(data_vector, id_, metadata):
            return (id_ in ids_to_keep) and \
                   (min_count <= data_vector.sum() <= max_count) and \
                   (min_nonzero <= (data_vector > 0).sum() <= max_nonzero)
    return f


def filter_samples_from_otu_table(otu_table, ids_to_keep, min_count, max_count,
                                  negate_ids_to_keep=False):
    filter_f = get_filter_function({}.fromkeys(ids_to_keep),
                                   min_count,
                                   max_count,
                                   0, inf, negate_ids_to_keep)
    return otu_table.filter(filter_f, axis='sample', inplace=False)


def filter_otus_from_otu_table(otu_table, ids_to_keep, min_count, max_count,
                               min_samples, max_samples,
                               negate_ids_to_keep=False):
    filter_f = get_filter_function({}.fromkeys(ids_to_keep),
                                   min_count,
                                   max_count,
                                   min_samples, max_samples,
                                   negate_ids_to_keep)
    return otu_table.filter(filter_f, axis='observation', inplace=False)

# end functions used by filter_samples_from_otu_table.py and
# filter_otus_from_otu_table.py


def filter_otu_table_to_n_samples(otu_table, n):
    """ Filter OTU table to n random samples.

        If n is greater than the number of samples or less than zero a
         ValueError will be raised.
    """
    if not (0 < n <= len(otu_table.ids())):
        raise ValueError("Number of samples to filter must be between 0 and "
                         "the number of samples.")
    return otu_table.subsample(n, axis='sample', by_id=True)


def filter_otus_from_otu_map(input_otu_map_fp,
                             output_otu_map_fp,
                             min_count,
                             min_sample_count=1):
    """ Filter otus with fewer than min_count sequences from input_otu_map_fp

        With very large data sets the number of singletons can be very large,
        and it becomes more efficent to filter them at the otu map stage than
        the otu table stage.

        There are two outputs from this function: the output file (which is the
        filtered otu map) and the list of retained otu ids as a set. Since I
        need to return the retained ids for pick_open_reference_otus, this
        takes filepaths instead of file handles (since it can't be a generator
        and return something).

    """
    results = set()
    output_otu_map_f = open(output_otu_map_fp, 'w')
    for line in open(input_otu_map_fp, 'U'):
        fields = line.strip().split('\t')
        sample_ids = set([e.split('_')[0] for e in fields[1:]])
        # only write this line if the otu has more than n sequences (so
        # greater than n tab-separated fields including the otu identifier)
        if (len(fields) > min_count) and (len(sample_ids) >= min_sample_count):
            output_otu_map_f.write(line)
            results.add(fields[0].split('\t')[0])
    output_otu_map_f.close()
    return results


def filter_tree(tree, tips_to_keep):
    result = tree.copy()
    # don't use this, it doesn't eliminate tips!
    # result = tree.getSubTree(tips_to_keep,ignore_missing=True)

    def f(node):
        if node.istip() and\
           node.Name is not None and\
           node.Name not in tips_to_keep and\
           node.Name.strip().strip('"').strip("'") not in tips_to_keep:
            return True
        return False
    result.removeDeleted(f)
    result.prune()
    return result
