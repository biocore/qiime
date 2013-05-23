#!/usr/bin/env python
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Justin Kuczynski", 
               "Jose Carlos Clemente Litran", 
               "Rob Knight", 
               "Greg Caporaso"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

"""Contains code for generating rarefied OTU tables at varying depth

this takes an otu table and generates a series of subsampled (without 
replacement) otu tables.
"""
import os.path
import numpy
from numpy import inf
from cogent.maths.stats.rarefaction import (subsample,
                                                  subsample_freq_dist_nonzero,
                                                  subsample_random,
                                                  subsample_multinomial)
from qiime.util import FunctionWithParams
from qiime.filter import filter_samples_from_otu_table, filter_otus_from_otu_table
from qiime.format import format_biom_table

class SingleRarefactionMaker(FunctionWithParams):
    def __init__(self, otu_path, depth):
        """ init a singlerarefactionmaker
        
        otu_path has to be parseable when opened by parse_biom_table,
        we just ignore any rarefaction levels beyond any sample in the data
        """
        self.depth = depth
        #self.otu_table = parse_biom_table(open(otu_path,'U'))
        self.otu_table = self.getBiomData(otu_path)

        self.max_num_taxa = -1
        for val in self.otu_table.iterObservationData():
            self.max_num_taxa = max(self.max_num_taxa, val.sum())

    def rarefy_to_file(self, output_fname, small_included=False,
        include_lineages=False,empty_otus_removed=False,subsample_f=subsample):
        """ computes rarefied otu tables and writes them, one at a time
        
        this prevents large memory usage
        
        for depth in self.rare_depths:
            for rep in range(self.num_reps):"""
        
        if not include_lineages:
            for (val, id, meta) in self.otu_table.iterObservations():
                del meta['taxonomy']

        sub_otu_table = get_rare_data(self.otu_table,
                                      self.depth, 
                                      small_included, 
                                      subsample_f=subsample_f)

        if empty_otus_removed:
            sub_otu_table = filter_otus_from_otu_table(sub_otu_table,
                                                       sub_otu_table.ObservationIds,
                                                       1, inf, 0, inf)

        self._write_rarefaction(output_fname, sub_otu_table)
    
    def _write_rarefaction(self, fname, sub_otu_table):
        """ depth and rep can be numbers or strings
        """
        if sub_otu_table.isEmpty():
            return
        f = open(fname, 'w')
        f.write(format_biom_table(sub_otu_table))
        f.close()


class RarefactionMaker(FunctionWithParams):
    def __init__(self, otu_path, min, max, step, num_reps):
        """ init a rarefactionmaker
        
        otu_path is path to .biom otu table
        we just ignore any rarefaction levels beyond any sample in the data
        """
        self.rare_depths = range(min,max+1, step)
        self.num_reps = num_reps
        #self.otu_table = parse_biom_table(open(otu_path,'U'))
        self.otu_table = self.getBiomData(otu_path)
        self.max_num_taxa = -1
        tmp = -1
        for val in self.otu_table.iterObservationData():
            if val.sum() > tmp:
                tmp = val.sum()
        self.max_num_taxa = tmp

    def rarefy_to_files(self, output_dir, small_included=False, 
        include_full=False, include_lineages=False,
        empty_otus_removed=False,subsample_f=subsample):
        """ computes rarefied otu tables and writes them, one at a time
        
        this prevents large memory usage"""
        if not include_lineages:
            for (val, id, meta) in self.otu_table.iterObservations():
                try:
                    del meta['taxonomy']
                except (TypeError,KeyError) as e:
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
                    sub_otu_table = filter_otus_from_otu_table(\
                        sub_otu_table,
                        sub_otu_table.ObservationIds, 1, inf, 0, inf)

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
        if sub_otu_table.isEmpty():
            return

        fname = 'rarefaction_'+str(depth)+'_'+str(rep)+'.biom'
        f = open(os.path.join(self.output_dir,fname), 'w')
        f.write(format_biom_table(sub_otu_table))
        f.close()

def get_rare_data(otu_table, 
                  seqs_per_sample, 
                  include_small_samples=False, 
                  subsample_f=subsample):
    """Filter OTU table to keep only desired sample sizes.
    
    - include_small_sampes=False => do not write samples with < seqs_per_sample
    total sequecnes
    - otu_table (input and out) is otus(rows) by samples (cols)
    - no otus are removed, even if they are absent in the rarefied table"""

    if not include_small_samples:
        otu_table = filter_samples_from_otu_table(otu_table, otu_table.SampleIds, seqs_per_sample, inf)

    # subsample samples that have too many sequences
    def func(x, s_id, s_md):
        if x.sum() < seqs_per_sample:
            return x
        else:
            return subsample_f(x, seqs_per_sample)

    subsampled_otu_table = otu_table.transformSamples(func)
    
    # remove small samples if required

    return subsampled_otu_table

# Not necessary anymore; tagged for removal
def remove_empty_otus(otu_mtx, otu_ids, otu_lineages=None):
    """ return matrix and otu_ids with otus of all 0's removed
    
    otu_mtx (in and out) is otus (rows) by samples (cols)"""
    nonempty_otu_idxs = []
    res_otu_ids = []
    for i in range(len(otu_ids)):
        if otu_mtx[i].sum() != 0:
            nonempty_otu_idxs.append(i)
            res_otu_ids.append(otu_ids[i])
    res_otu_mtx = otu_mtx[nonempty_otu_idxs,:]

    if otu_lineages == None or otu_lineages == []:
        res_otu_lineages = []
    else:
        res_otu_lineages = [otu_lineages[i] for i in nonempty_otu_idxs]

    return res_otu_mtx, res_otu_ids, res_otu_lineages

