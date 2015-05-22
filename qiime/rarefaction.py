#!/usr/bin/env python
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Jose Carlos Clemente Litran", "Rob Knight",
               "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""Contains code for generating rarefied OTU tables at varying depth

this takes an otu table and generates a series of subsampled (without
replacement) otu tables.
"""
import os.path

import numpy
from numpy import inf
from skbio.stats import subsample
from biom.err import errstate

from qiime.util import FunctionWithParams, write_biom_table
from qiime.filter import (filter_samples_from_otu_table,
                          filter_otus_from_otu_table)


class SingleRarefactionMaker(FunctionWithParams):

    def __init__(self, otu_path, depth):
        """ init a singlerarefactionmaker

        otu_path has to be a parseable BIOM table,
        we just ignore any rarefaction levels beyond any sample in the data
        """
        self.depth = depth
        self.otu_table = self.getBiomData(otu_path)

        self.max_num_taxa = -1
        for val in self.otu_table.iter_data(axis='observation'):
            self.max_num_taxa = max(self.max_num_taxa, val.sum())

    def rarefy_to_file(self, output_fname, small_included=False,
                       include_lineages=False, empty_otus_removed=False, subsample_f=subsample):
        """ computes rarefied otu tables and writes them, one at a time

        this prevents large memory usage

        for depth in self.rare_depths:
            for rep in range(self.num_reps):"""

        if not include_lineages:
            for (val, id, meta) in self.otu_table.iter(axis='observation'):
                del meta['taxonomy']

        sub_otu_table = get_rare_data(self.otu_table,
                                      self.depth,
                                      small_included,
                                      subsample_f=subsample_f)

        if empty_otus_removed:
            sub_otu_table = filter_otus_from_otu_table(
                sub_otu_table, sub_otu_table.ids(axis='observation'),
                1, inf, 0, inf)

        self._write_rarefaction(output_fname, sub_otu_table)

    def _write_rarefaction(self, fname, sub_otu_table):
        """ depth and rep can be numbers or strings
        """
        if sub_otu_table.is_empty():
            return
        write_biom_table(sub_otu_table, fname)


class RarefactionMaker(FunctionWithParams):

    def __init__(self, otu_path, min, max, step, num_reps):
        """ init a rarefactionmaker

        otu_path is path to .biom otu table
        we just ignore any rarefaction levels beyond any sample in the data
        """
        self.rare_depths = range(min, max + 1, step)
        self.num_reps = num_reps
        self.otu_table = self.getBiomData(otu_path)
        self.max_num_taxa = -1
        tmp = -1
        for val in self.otu_table.iter_data(axis='observation'):
            if val.sum() > tmp:
                tmp = val.sum()
        self.max_num_taxa = tmp

    def rarefy_to_files(self, output_dir, small_included=False,
                        include_full=False, include_lineages=False,
                        empty_otus_removed=False, subsample_f=subsample):
        """ computes rarefied otu tables and writes them, one at a time

        this prevents large memory usage"""
        if not include_lineages:
            for (val, id, meta) in self.otu_table.iter(axis='observation'):
                try:
                    del meta['taxonomy']
                except (TypeError, KeyError) as e:
                    # no meta or just no taxonomy present
                    pass

        self.output_dir = output_dir
        for depth in self.rare_depths:
            for rep in range(self.num_reps):
                sub_otu_table = get_rare_data(self.otu_table,
                                              depth,
                                              small_included,
                                              subsample_f=subsample_f)
                if empty_otus_removed:
                    sub_otu_table = filter_otus_from_otu_table(
                        sub_otu_table, sub_otu_table.ids(axis='observation'),
                        1, inf, 0, inf)

                self._write_rarefaction(depth, rep, sub_otu_table)

        if include_full:
            self._write_rarefaction('full', 0, self.otu_table)

    def rarefy_to_list(self, small_included=False, include_full=False,
                       include_lineages=False):
        """ computes rarefied otu tables and returns a list

        each element
        is (depth, rep, sample_ids, taxon_ids, otu_table)
        depth is string "full" for one instance
        """
        if include_lineages:
            otu_lineages = self.lineages
        else:
            otu_lineages = None
        res = []
        for depth in self.rare_depths:
            for rep in range(self.num_reps):
                sub_otu_table = get_rare_data(self.otu_table,
                                              depth, small_included)
                res.append([depth, rep, sub_otu_table])

        if include_full:
            res.append(['full', 0, self.otu_table])
        return res

    def _write_rarefaction(self, depth, rep, sub_otu_table):
        """ depth and rep can be numbers or strings
        """
        if sub_otu_table.is_empty():
            return

        fname = 'rarefaction_' + str(depth) + '_' + str(rep) + '.biom'
        fname = os.path.join(self.output_dir, fname)
        write_biom_table(sub_otu_table, fname)


def get_rare_data(otu_table,
                  seqs_per_sample,
                  include_small_samples=False,
                  subsample_f=subsample):
    """Filter OTU table to keep only desired sample sizes.

    - include_small_sampes=False => do not write samples with < seqs_per_sample
    total sequecnes
    - otu_table (input and out) is otus(rows) by samples (cols)
    - no otus are removed, even if they are absent in the rarefied table"""

    with errstate(empty='raise'):
        if not include_small_samples:
            otu_table = filter_samples_from_otu_table(
                otu_table,
                otu_table.ids(),
                seqs_per_sample,
                inf)

        # subsample samples that have too many sequences
        def func(x, s_id, s_md):
            if x.sum() < seqs_per_sample:
                return x
            else:
                return subsample_f(x.astype(int), seqs_per_sample)

        subsampled_otu_table = otu_table.transform(func, axis='sample')

        return subsampled_otu_table
