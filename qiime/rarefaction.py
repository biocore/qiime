#!/usr/bin/env python
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Justin Kuczynski", "Jose Carlos Clemente Litran"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.4.0-dev"
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
from numpy import inf
from qiime.filter import filter_samples_from_otu_table, filter_otus_from_otu_table
from qiime.pycogent_backports.parse_biom import parse_biom_table

class SingleRarefactionMaker(FunctionWithParams):
    def __init__(self, otu_path, depth):
        """ init a singlerarefactionmaker
        
        otu_path has to be parseable when opened by parse_biom_table,
        we just ignore any rarefaction levels beyond any sample in the data
        """
#        """ init a singlerarefactionmaker
#        
#        otu_path has to be parseable when opened by parse_otu_table,
#        or it can be a 4-tuple.
#        we just ignore any rarefaction levels beyond any sample in the data
#        """

        # No tuples anymore: file paths only

        self.depth = depth

        #if type(otu_path) == type(('a',1)):
        #    self.sample_names, self.taxon_names, \
        #    self.otu_table, self.lineages = otu_path

        #else:
        #    self.sample_names, self.taxon_names, \
        #    self.otu_table, self.lineages = parse_otu_table(open(otu_path,'U'))
        self.otu_table = parse_biom_table(open(otu_path,'U'))

        #self.max_num_taxa = (self.otu_table.sum(1)).max()
        self.max_num_taxa = -1
        for val in self.otu_table.iterObservationData():
            self.max_num_taxa = max(self.max_num_taxa, val.sum())

    def rarefy_to_file(self, output_fname, small_included=False,
        include_lineages=False,empty_otus_removed=False):
        """ computes rarefied otu tables and writes them, one at a time
        
        this prevents large memory usage
        
        for depth in self.rare_depths:
            for rep in range(self.num_reps):"""
            
        #if include_lineages:
        #    sub_otu_lineages = self.otu_table.ObservationMetadata['taxonomy']
        #else:
        #    sub_otu_lineages = None
        if not include_lineages:
            for (val, id, meta) in self.otu_table.iterObservations():
                del meta['taxonomy']

        #sub_sample_ids, sub_otu_table = get_rare_data(self.sample_names,
        #    self.otu_table, self.depth, small_included)
        sub_otu_table = get_rare_data(self.otu_table,
                                      self.depth, small_included)
        #sub_otu_ids = self.taxon_names
        #if empty_otus_removed:
        #    sub_otu_table, sub_otu_ids, sub_otu_lineages = \
        #        remove_empty_otus(sub_otu_table, 
        #        sub_otu_ids, sub_otu_lineages) 
        if empty_otus_removed:
            sub_otu_table = filter_otus_from_otu_table(sub_otu_table,
                                                       sub_otu_table.ObservationIds,
                                                       1, inf, 0, inf)
                # sub_otu_lineages can be None or []
        #self._write_rarefaction(output_fname, sub_otu_table.SampleIds, sub_otu_ids,
        #    sub_otu_table, sub_otu_lineages)
        self._write_rarefaction(output_fname, sub_otu_table)
    
    #def _write_rarefaction(self, fname, sub_sample_ids, sub_otu_ids,\
    #    sub_otu_table, otu_lineages):
    def _write_rarefaction(self, fname, sub_otu_table):
        """ depth and rep can be numbers or strings
        """
        #if min(numpy.shape(sub_otu_table)) == 0: # no data to write
        #    return
        #f = open(fname, 'w')
        #f.write(format_otu_table(sub_sample_ids, sub_otu_ids,\
        #    sub_otu_table, otu_lineages, comment=fname))

        if sub_otu_table.isEmpty():
            return
        f = open(fname, 'w')
        f.write(sub_otu_table.getBiomFormatJsonString())
        f.close()


class RarefactionMaker(FunctionWithParams):
    def __init__(self, otu_path, min, max, step, num_reps):
        """ init a rarefactionmaker
        
        otu_path is path to .biom otu table
        we just ignore any rarefaction levels beyond any sample in the data
        """
        self.rare_depths = range(min,max+1, step)
        self.num_reps = num_reps
        self.otu_table = parse_biom_table(open(otu_path,'U'))
        self.max_num_taxa = -1
        tmp = -1
        for val in self.otu_table.iterObservationData():
            if val.sum() > tmp:
                tmp = val.sum()
            #self.max_num_taxa = max(self.max_num_taxa, val.sum())
        self.max_num_taxa = tmp

    def rarefy_to_files(self, output_dir, small_included=False, 
        include_full=False, include_lineages=False,
        empty_otus_removed=False):
        """ computes rarefied otu tables and writes them, one at a time
        
        this prevents large memory usage"""
        #if include_lineages:
        #    all_otu_lineages = self.lineages
        #else:
        #    all_otu_lineages = None
        if not include_lineages:
            for (val, id, meta) in self.otu_table.iterObservations():
                del meta['taxonomy']

        self.output_dir = output_dir
        for depth in self.rare_depths:
            for rep in range(self.num_reps):
                sub_otu_table = get_rare_data(self.otu_table, depth, 
                                              small_included)
                #sub_otu_ids = self.taxon_names
                #if empty_otus_removed:
                #    sub_otu_table, sub_otu_ids, sub_otu_lineages =\
                #        remove_empty_otus(sub_otu_table,sub_otu_ids,
                #        all_otu_lineages)
                if empty_otus_removed:
                    sub_otu_table = filter_otus_from_otu_table(\
                        sub_otu_table,
                        sub_otu_table.ObservationIds, 1, inf, 0, inf)

                self._write_rarefaction(depth, rep, sub_otu_table)
                #self._write_rarefaction(depth, rep, sub_otu_table.SampleIds, 
                #    sub_otu_ids, sub_otu_table, sub_otu_lineages)

        if include_full:
            #self._write_rarefaction('full', 0, self.otu_table.SampleIds, \
            #    self.ObservationIds,self.otu_table, sub_otu_lineages)
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
                #res.append([depth, rep, sub_otu_table.SampleIds, \
                #            self.ObservationIds, sub_otu_table, otu_lineages])
                res.append([depth, rep, sub_otu_table])

        if include_full:
            #res.append(['full', 0, self.SampleIds, self.ObservationIds,\
            #    self.otu_table, otu_lineages])
            res.append(['full', 0, self.otu_table])
        return res
    
    #def _write_rarefaction(self, depth, rep, sub_sample_ids, sub_otu_ids,\
    #    sub_otu_table, otu_lineages):
    def _write_rarefaction(self, depth, rep, sub_otu_table):
        """ depth and rep can be numbers or strings
        """
        #if min(numpy.shape(sub_otu_table)) == 0: # no data to write
        #    return
        if sub_otu_table.isEmpty():
            return

        fname = 'rarefaction_'+str(depth)+'_'+str(rep)+'.txt'
        f = open(os.path.join(self.output_dir,fname), 'w')
        f.write(sub_otu_table.getBiomFormatJsonString())
        #f.write(format_otu_table(sub_sample_ids, sub_otu_ids,\
        #    sub_otu_table, otu_lineages, comment=fname))
        f.close()

#def get_rare_data(sample_ids, otu_table,
#    seqs_per_sample, include_small_samples=False):
def get_rare_data(otu_table, seqs_per_sample, include_small_samples=False):
    """Filter OTU table to keep only desired sample sizes.
    
    - include_small_sampes=False => do not write samples with < seqs_per_sample
    total sequecnes
    - otu_table (input and out) is otus(rows) by samples (cols)
    - no otus are removed, even if they are absent in the rarefied table"""
    #res_otu_table = otu_table.copy()
    #res_sample_ids = sample_ids
    #res_sample_ids = otu_table.ObservationIds

    # subsample samples that have too many sequences
    def func(x):
        if x.sum() < seqs_per_sample:
            return x
        else:
            return subsample(x, seqs_per_sample)

    subsampled_otu_table = otu_table.transformSamples(func)
    
    # remove small samples if required
    if not include_small_samples:
        subsampled_otu_table = filter_samples_from_otu_table(subsampled_otu_table, subsampled_otu_table.SampleIds, seqs_per_sample, inf)

    return subsampled_otu_table

    #figure out which samples will be dropped because too small
    #too_big_samples = (otu_table.sum(0)>seqs_per_sample).nonzero()[0]
    #if too_big_samples.shape[0]:    #means that there were some
    #    for i in too_big_samples:
    #        res_otu_table[:,i] = subsample(otu_table[:,i].ravel(), 
    #            seqs_per_sample)

    #if not include_small_samples:
    #    big_enough_samples = (res_otu_table.sum(0)>=seqs_per_sample).nonzero()
    #    res_otu_table = res_otu_table[:,big_enough_samples[0]]
    #    res_sample_ids = map(sample_ids.__getitem__, big_enough_samples[0])
    #figure out which samples will be reduced because too big

    #return res_sample_ids, res_otu_table

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

