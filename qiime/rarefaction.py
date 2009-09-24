#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

"""Contains code for generating rarefied OTU tables at varying depth

example:
python rarefaction.py TEST/otu_table.txt -m 700 -x 700 -s 1 -n 100 -d TEST/beta_rare/

this takes an otu table and generates a series of subsampled (without 
replacement) otu tables.
"""
from optparse import OptionParser
from qiime.parse import parse_otus
from qiime.format import format_otu_table
from qiime.util import FunctionWithParams
import os.path
import os
import numpy
from cogent.maths.stats.rarefaction import subsample


class RarefactionMaker(FunctionWithParams):
    def __init__(self, otu_path, min, max, step, num_reps):
        """ otu_path can be path or otu tuple, defined by util.py's getOtuTable
        we just ignore any rarefaction levels beyond any sample in the data
        """
        self.rare_depths = range(min,max+1, step)
        self.num_reps = num_reps
        self.sample_names, self.taxon_names, self.otu_table, self.lineages = \
            self.getOtuTable(otu_path)
        self.max_num_taxa = (self.otu_table.sum(1)).max()


    def rarefy_to_files(self, output_dir, small_included=False, include_full=False):
        """ computes rarefied otu tables and writes them, one at a time
        this prevents large memory usage"""
        self.output_dir = output_dir
        for depth in self.rare_depths:
            for rep in range(self.num_reps):
                sub_sample_ids, sub_otu_ids, sub_otu_table = \
                get_rare_data(self.sample_names, self.taxon_names,
                    self.otu_table, depth, small_included)
                self._write_rarefaction(depth, rep, sub_sample_ids, sub_otu_ids,\
                    sub_otu_table)

        if include_full:
            self._write_rarefaction('full', 0, self.sample_names, \
                self.taxon_names,\
                self.otu_table)
    
    def rarefy_to_list(self, small_included=False, include_full=False):
        """ computes rarefied otu tables and returns a list, each element
        is (depth, rep, sample_ids, taxon_ids, otu_table)
        depth is string "full" for one instance
        """
        res = []
        for depth in self.rare_depths:
            for rep in range(self.num_reps):
                sub_sample_ids, sub_otu_ids, sub_otu_table = \
                get_rare_data(self.sample_names, self.taxon_names,
                    self.otu_table, depth, small_included)
                res.append([depth, rep, sub_sample_ids, sub_otu_ids, \
                    sub_otu_table])

        if include_full:
            res.append(['full', 0, self.sample_names, self.taxon_names,\
                self.otu_table])
        return res
    
    def _write_rarefaction(self, depth, rep, sub_sample_ids, sub_otu_ids,\
        sub_otu_table):
        """ depth and rep can be numbers or strings
        """
        if min(numpy.shape(sub_otu_table)) == 0: # no data to write
            return
        fname = 'rarefaction_'+str(depth)+'_'+str(rep)+'.txt'
        f = open(os.path.join(self.output_dir,fname), 'w')
        f.write(format_otu_table(sub_sample_ids, sub_otu_ids,\
            sub_otu_table, None, comment=fname))
        f.close()

def get_rare_data(sample_ids, otu_ids, otu_table, \
	seqs_per_sample, include_small_samples=False):
    """Filter OTU table to keep only desired sample sizes.
    
    include_small_sampes=False => 
    sample has 100 seqs, rarefy """
    res_otu_table = otu_table.copy()
    res_sample_ids = sample_ids
    #figure out which samples will be dropped because too small
    too_big_samples = (otu_table.sum(0)>seqs_per_sample).nonzero()[0]
    if too_big_samples.shape[0]:    #means that there were some
        for i in too_big_samples:
            res_otu_table[:,i] = subsample(otu_table[:,i].ravel(), seqs_per_sample)
    if not include_small_samples:
        big_enough_samples = (res_otu_table.sum(0)>=seqs_per_sample).nonzero()
        res_otu_table = res_otu_table[:,big_enough_samples[0]]
        res_sample_ids = map(sample_ids.__getitem__, big_enough_samples[0])
    #figure out which samples will be reduced because too big
    return res_sample_ids, otu_ids, res_otu_table

def make_cmd_parser():
    """returns command-line options"""
    usage = '%prog input_filepath [options]\n' +\
    'use %prog -h for help.'
    description=\
        'input_filepath is an otu table, output will be a '+\
        'series of otu table files, named e.g.: rarefaction_4_2.txt '+\
        '(4 seqs per sample, iteration 2 (3rd such file written)\n'+\
        'example: %prog otu_table.txt -m 10 -x 100 -s 10 -n 2'
    parser = OptionParser(usage=usage, description=description)
    parser.add_option('-d', '--dest', dest='output_dir', 
        default=os.getcwd(),
        help='write output rarefied otu tables here (directory) '+\
            "default current dir.  makes dir if it doesn't exist")
    parser.add_option('-m', '--min', dest='min', default=1, type=int,
        help='min seqs/sample, default 1')
    parser.add_option('-x', '--max', dest='max', default=10, type=int,
        help='max seqs/sample (inclusive), default 10')
    parser.add_option('-s', '--step', dest='step', default=1, type=int,
        help='levels: min, min+step... for level <= max, default 1')
    parser.add_option('-n', '--num-reps', dest='num_reps', default=1, type=int,
        help='num reps at each seqs/sample level, default 1')
    parser.add_option('--small_included', dest='small_included', default=False, 
    	action="store_true",
        help="""samples containing fewer seqs than the rarefaction
level are included in the output but not rarefied""")
    options, args = parser.parse_args()
    if (len(args) != 1) :
        parser.error("incorrect number of arguments")
    return options, args


if __name__ == '__main__':
    #if called from command-line, should run the analysis.
    opts, args = make_cmd_parser()
    otu_path = args[0]
    if not os.path.exists(opts.output_dir):
    	os.mkdir(opts.output_dir)
    maker = RarefactionMaker(otu_path, opts.min, opts.max,
        opts.step, opts.num_reps)
    maker.rarefy_to_files(opts.output_dir, opts.small_included)
