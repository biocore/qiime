#!/usr/bin/env python
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

"""Contains code for generating rarefied OTU tables at varying depth

this takes an otu table and generates a series of subsampled (without 
replacement) otu tables.
"""
from qiime.format import format_otu_table
from qiime.util import FunctionWithParams
import os.path
import numpy
from cogent.maths.stats.rarefaction import subsample

class SingleRarefactionMaker(FunctionWithParams):
    def __init__(self, otu_path, depth):
        """ init a singlerarefactionmaker
        
        otu_path can be path or otu tuple, defined by util.py's getOtuTable
        we just ignore any rarefaction levels beyond any sample in the data
        """
        
        self.depth = depth

        self.sample_names, self.taxon_names, self.otu_table, self.lineages = \
            self.getOtuTable(otu_path)
        self.max_num_taxa = (self.otu_table.sum(1)).max()


    def rarefy_to_file(self, output_fname, small_included=False,
        include_lineages=False):
        """ computes rarefied otu tables and writes them, one at a time
        
        this prevents large memory usage
        
        for depth in self.rare_depths:
            for rep in range(self.num_reps):"""
            
        if include_lineages:
            otu_lineages = self.lineages
        else:
            otu_lineages = None
        sub_sample_ids, sub_otu_table = get_rare_data(self.sample_names,
            self.otu_table, self.depth, small_included)
        self._write_rarefaction(output_fname, sub_sample_ids, self.taxon_names,
            sub_otu_table, otu_lineages)
    
    def _write_rarefaction(self, fname, sub_sample_ids, sub_otu_ids,\
        sub_otu_table, otu_lineages):
        """ depth and rep can be numbers or strings
        """
        if min(numpy.shape(sub_otu_table)) == 0: # no data to write
            return
        f = open(fname, 'w')
        f.write(format_otu_table(sub_sample_ids, sub_otu_ids,\
            sub_otu_table, otu_lineages, comment=fname))
        f.close()


class RarefactionMaker(FunctionWithParams):
    def __init__(self, otu_path, min, max, step, num_reps):
        """ init a rarefactionmaker
        
        otu_path can be path or otu tuple, defined by util.py's getOtuTable
        we just ignore any rarefaction levels beyond any sample in the data
        """
        self.rare_depths = range(min,max+1, step)
        self.num_reps = num_reps
        self.sample_names, self.taxon_names, self.otu_table, self.lineages = \
            self.getOtuTable(otu_path)
        self.max_num_taxa = (self.otu_table.sum(1)).max()


    def rarefy_to_files(self, output_dir, small_included=False, 
        include_full=False, include_lineages=False):
        """ computes rarefied otu tables and writes them, one at a time
        
        this prevents large memory usage"""
        if include_lineages:
            otu_lineages = self.lineages
        else:
            otu_lineages = None
        self.output_dir = output_dir
        for depth in self.rare_depths:
            for rep in range(self.num_reps):
                sub_sample_ids, sub_otu_table = \
                get_rare_data(self.sample_names, self.otu_table, depth, 
                  small_included)
                self._write_rarefaction(depth, rep, sub_sample_ids, 
                    self.taxon_names, sub_otu_table, otu_lineages)

        if include_full:
            self._write_rarefaction('full', 0, self.sample_names, \
                self.taxon_names,self.otu_table, otu_lineages)
    
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
                sub_sample_ids, sub_otu_table = \
                get_rare_data(self.sample_names,
                    self.otu_table, depth, small_included)
                res.append([depth, rep, sub_sample_ids, self.taxon_names, \
                    sub_otu_table, otu_lineages])

        if include_full:
            res.append(['full', 0, self.sample_names, self.taxon_names,\
                self.otu_table, otu_lineages])
        return res
    
    def _write_rarefaction(self, depth, rep, sub_sample_ids, sub_otu_ids,\
        sub_otu_table, otu_lineages):
        """ depth and rep can be numbers or strings
        """
        if min(numpy.shape(sub_otu_table)) == 0: # no data to write
            return
        fname = 'rarefaction_'+str(depth)+'_'+str(rep)+'.txt'
        f = open(os.path.join(self.output_dir,fname), 'w')
        f.write(format_otu_table(sub_sample_ids, sub_otu_ids,\
            sub_otu_table, otu_lineages, comment=fname))
        f.close()

def get_rare_data(sample_ids, otu_table, \
    seqs_per_sample, include_small_samples=False):
    """Filter OTU table to keep only desired sample sizes.
    
    - include_small_sampes=False => do not write samples with < seqs_per_sample
    total sequecnes
    - otu_table (input and out) is otus(rows) by samples (cols)
    - no otus are removed, even if they are absent in the rarefied table"""
    res_otu_table = otu_table.copy()
    res_sample_ids = sample_ids
    #figure out which samples will be dropped because too small
    too_big_samples = (otu_table.sum(0)>seqs_per_sample).nonzero()[0]
    if too_big_samples.shape[0]:    #means that there were some
        for i in too_big_samples:
            res_otu_table[:,i] = subsample(otu_table[:,i].ravel(), 
                seqs_per_sample)
    if not include_small_samples:
        big_enough_samples = (res_otu_table.sum(0)>=seqs_per_sample).nonzero()
        res_otu_table = res_otu_table[:,big_enough_samples[0]]
        res_sample_ids = map(sample_ids.__getitem__, big_enough_samples[0])
    #figure out which samples will be reduced because too big
    return res_sample_ids, res_otu_table

