#!/usr/bin/env python
# File created on 19 Mar 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight",
               "Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from os.path import join
from operator import add
from datetime import datetime
from random import choice

from numpy import zeros
from biom.table import Table

from qiime.format import format_mapping_file
from qiime.parse import parse_mapping_file
from qiime.util import (make_option, create_dir, parse_command_line_parameters,
                        write_biom_table, get_generated_by_for_biom_tables)
from qiime.sort import natsort


def sim_otu_table(sample_ids, otu_ids, samples, otu_metadata, tree,
                  num_replicates, dissimilarity):
    """ make n samples related to each sample in an input otu table

    input:
    the constituents of an otu table
     * sample_ids: list
     * otu_ids: list
     * samples: iterable, each element must have sample_vector = elem[0]
       where sample_vector looks like [1,0,7,9...]
     * otu_metadata: list, either empty or of length len(otu_ids)
    tree: a PhyloNode tree or compatible,
    num_replicates: how many replicates for each input sample
    dissimilarity: how phylogenetically dissimilar each replicate should be
    relative to the original

    output is a tuple with the constituents of an otu table, possibly missing
    some otu metadata:
    (res_sam_names, res_otus, res_otu_mtx, res_otu_metadata)
    """
    # hold sample abundance vector in a dict (sample_dict) temporarily
    # format of sample_dict: otu_id: num_seqs
    sample_dicts = []
    res_sam_names = []
    for i, sample_info in enumerate(samples):
        sample_vector = sample_info[0]
        # sample_vector = otu_observations.next()
        for j in range(num_replicates):
            sample_dict = {}
            for k in range(len(otu_ids)):
                otu_abundance = sample_vector[k]
                if otu_abundance == 0:
                    continue
                new_otu_id = get_new_otu_id(otu_ids[k], tree, dissimilarity)
                # beware, get_new_otu_id may return something we already
                # have in the table
                if new_otu_id in sample_dict:
                    sample_dict[new_otu_id] += otu_abundance
                else:
                    sample_dict[new_otu_id] = otu_abundance
            sample_dicts.append(sample_dict)
            res_sam_names.append(sample_ids[i] + '.' + str(j))

    res_otu_mtx, res_otus = combine_sample_dicts(sample_dicts)

    res_otu_metadata = []
    if otu_metadata is None or otu_metadata == []:
        res_otu_metadata = None
    else:
        for otu_id in res_otus:
            # if otu was in original table, just copy its metadata
            try:
                res_otu_metadata.append(otu_metadata[otu_ids.index(otu_id)])
            except ValueError:
                # else just append None since we don't have its metadata
                res_otu_metadata.append(None)

    return res_sam_names, res_otus, res_otu_mtx, res_otu_metadata


def create_tip_index(tree):
    """Create a tip lookup index on the tree"""
    if hasattr(tree, '_tip_index'):
        return
    else:
        tree._tip_index = {n.Name: n for n in tree.tips()}

def cache_tip_names(tree):
    """Cache tip names"""
    if hasattr(tree, '_tip_names'):
        return
    else:
        for n in tree.postorder():
            if n.isTip():
                n._tip_names = [n.Name]
            else:
                n._tip_names = reduce(add, [c._tip_names for c in n.Children])


def get_new_otu_id(old_otu_id, tree, dissim):
    """ simulates an otu switching to related one

    input a tipname, a tree, and a distance to walk up the tree
    ouputs the name of the new, randomly chosen, tree tip
    output tip name may be the same as
    """
    create_tip_index(tree)
    cache_tip_names(tree)

    node = tree._tip_index[old_otu_id]
    distance_up_tree = 0
    while (not node.isRoot()) and (distance_up_tree + node.Length < dissim):
        distance_up_tree += node.Length
        node = node.Parent

    # another option is to do 50-50 descending each node,
    # so we don't bias for overrepresented clades
    if node.isTip():
        return node.Name
    else:
        return choice(node._tip_names)


def combine_sample_dicts(sample_dicts):
    """ combines a list of sample_dicts into one otu table

    sample dicts is a list of dicts, each one {otu_id:num_seqs}

    output is a tuple:
    (otu_mtx (rows are otus), otu_ids (list))
    * otu_mtx has samples in order of dicts, otus sorted with natsort
    / human sort
    * otu_mtx will have all otus mentioned as keys in sample_dicts, even if
    they are abundance 0  ({otu_id:0,...})
    such otus will simply be rows of zeros
    """
    all_otu_ids = []
    for s in sample_dicts:
        all_otu_ids.extend(s.keys())

    all_otu_ids = list(set(all_otu_ids))
    all_otu_ids = natsort(all_otu_ids)

    # get index once now, for all samples, instead of all_otu_ids.index()
    indices = {}
    for i in range(len(all_otu_ids)):
        indices[all_otu_ids[i]] = i

    otu_mtx = zeros((len(all_otu_ids), len(sample_dicts)), int)
    # otus (rows) by samples (cols)
    for i, sample_dict in enumerate(sample_dicts):
        for otu, abund in sample_dict.items():
            otu_mtx[indices[otu], i] = abund

    return otu_mtx, all_otu_ids


def create_replicated_mapping_file(map_f, num_replicates, sample_ids):
    """Returns a formatted mapping file with replicated sample IDs.

    Each sample ID will have an ascending integer appended to it from the range
    [0, num_replicates - 1]. For example, if there are two input sample IDs, S1
    and S2, with 3 replicates each, the output will be:
        S1.0
        S1.1
        S1.2
        S2.0
        S2.1
        S2.2

    All other metadata columns will simply be copied to the output mapping
    file. The order of input sample IDs is preserved.

    Arguments:
        map_f - input mapping file to replicate (file-like object)
        num_replicates - number of replicates at each sample
        sample_ids - only sample IDs in the mapping file that are in this list
            will be replicated. Sample IDs in the mapping file that are not
            found in this list will not be added to the resulting mapping file
    """
    if num_replicates < 1:
        raise ValueError("Must specify at least one sample replicate (was "
                         "provided %d)." % num_replicates)
    map_data, header, comments = parse_mapping_file(map_f)

    rep_map_data = []
    for row in map_data:
        if row[0] in sample_ids:
            for rep_num in range(num_replicates):
                rep_map_data.append(['%s.%i' % (row[0], rep_num)] + row[1:])

    return format_mapping_file(header, rep_map_data, comments)


def simsam_range(table,
                 tree,
                 simulated_sample_sizes,
                 dissimilarities,
                 mapping_f=None):
    """Applies sim_otu_table over a range of parameters

     table: the input table to simulate samples from
     tree: tree related OTUs in input table
     simulated_sample_sizes: a list of ints defining how many
      output samples should be create per input sample
     dissimilarities: a list of floats containing the
      dissimilarities to use in simulating tables
     mapping_f: file handle for metadata mapping file, if
      a mapping file should be created with the samples from
      each simulated table

     This function will yield tuples with the following form:
      (output table, output mapping lines, simulated_sample_size, dissimilarity)

     If the user does not provide mapping_f, the tuples will look like:
      (output table, None, simulated_sample_size, dissimilarity)

    """
    if mapping_f is not None:
        # if the user provided a mapping file, load it into
        # a list for repeated use, and define the function for
        # processing the mapping file
        mapping_lines = list(mapping_f)
        process_map = create_replicated_mapping_file
    else:
        # otherwise create a dummy function for processing the
        # mapping file so we don't have to check whether it
        # exists on every iteration
        mapping_lines = None

        def process_map(mapping_lines, simulated_sample_size, sample_ids):
            return None

    for simulated_sample_size in simulated_sample_sizes:
        # create the output mapping file data
        output_mapping_lines = \
            process_map(mapping_lines, simulated_sample_size, table.ids())
        for dissimilarity in dissimilarities:
            # create the simulated otu table
            output_sample_ids, output_otu_ids, output_data, output_metadata = \
                sim_otu_table(table.ids(),
                              table.ids(axis='observation').tolist(),
                              table.iter(),
                              table.metadata(axis='observation'),
                              tree,
                              simulated_sample_size,
                              dissimilarity)
            output_table = Table(
                output_data, output_otu_ids, output_sample_ids,
                observation_metadata=output_metadata,
                generated_by=get_generated_by_for_biom_tables(),
                create_date=datetime.now().isoformat())
            yield (output_table,
                   output_mapping_lines,
                   simulated_sample_size,
                   dissimilarity)


def simsam_range_to_files(table,
                          tree,
                          simulated_sample_sizes,
                          dissimilarities,
                          output_dir,
                          mapping_f=None,
                          output_table_basename="table",
                          output_map_basename="map"):
    """Applies sim_otu_table over a range of parameters, writing output to file

     table: the input table to simulate samples from
     tree: tree related OTUs in input table
     simulated_sample_sizes: a list of ints defining how many
      output samples should be create per input sample
     dissimilarities: a list of floats containing the
      dissimilarities to use in simulating tables
     output_dir: the directory where all output tables and
      mapping files should be written
     mapping_f: file handle for metadata mapping file, if
      a mapping file should be created with the samples from
      each simulated table
     output_table_basename: basename for output table files
      (default: table)
     output_map_basename: basename for output mapping files
      (default: map)
    """
    create_dir(output_dir)
    for e in simsam_range(table, tree, simulated_sample_sizes, dissimilarities, mapping_f):
        output_table = e[0]
        output_mapping_lines = e[1]
        simulated_sample_size = e[2]
        dissimilarity = e[3]

        output_table_fp = join(output_dir, '%s_n%d_d%r.biom' %
                               (output_table_basename, simulated_sample_size, dissimilarity))
        write_biom_table(output_table, output_table_fp)

        if output_mapping_lines is not None:
            output_map_fp = join(output_dir, '%s_n%d_d%r.txt' %
                                 (output_map_basename, simulated_sample_size, dissimilarity))
            output_map_f = open(output_map_fp, 'w')
            output_map_f.write(''.join(output_mapping_lines))
            output_map_f.close()
