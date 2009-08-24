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
from qiime.parse import parse_otus, filter_otus_by_lineage
from qiime.format import format_otu_table
from qiime.util import FunctionWithParams
import os.path
import os
import numpy


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


    def rarefy_to_files(self, output_dir, include_full=False):
        """ computes rarefied otu tables and writes them, one at a time
        this prevents large memory usage"""
        self.output_dir = output_dir

        for depth in self.rare_depths:
            for rep in range(self.num_reps):
                sub_sample_ids, sub_otu_ids, sub_otu_table, sub_lineages = \
                filter_otus_by_lineage(self.sample_names, self.taxon_names,
                    self.otu_table, lineages=None, wanted_lineage=None,
                    max_seqs_per_sample=depth, min_seqs_per_sample=depth)
                self._write_rarefaction(depth, rep, sub_sample_ids, sub_otu_ids,\
                    sub_otu_table)

        if include_full:
            self._write_rarefaction('full', 0, self.sample_names, \
                self.taxon_names,\
                self.otu_table)
    
    def rarefy_to_list(self, include_full=False):
        """ computes rarefied otu tables and returns a list, each element
        is (depth, rep, sample_ids, taxon_ids, otu_table)
        depth is string "full" for one instance
        """
        res = []
        for depth in self.rare_depths:
            for rep in range(self.num_reps):
                sub_sample_ids, sub_otu_ids, sub_otu_table, sub_lineages = \
                filter_otus_by_lineage(self.sample_names, self.taxon_names,
                    self.otu_table, lineages=None, wanted_lineage=None,
                    max_seqs_per_sample=depth, min_seqs_per_sample=depth)
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


def make_cmd_parser():
    """returns command-line options"""
    usage = '%prog input_filepath [options]\n' +\
    'use %prog -h for help.\n\n'+\
        'input_filepath is an otu table, output will be a '+\
        'series of otu table files, named e.g.: rarefaction_4_2.txt '+\
        '(4 seqs per sample, iteration 2 (3rd such file written)\n'+\
        'example: %prog otu_table.txt -m 10 -x 100 -s 10 -n 2'
    parser = OptionParser(usage=usage)
    parser.add_option('-d', '--dest', dest='output_dir', 
        default=os.getcwd(),
        help='write output rarefied otu tables here (directory) '+\
            'default current dir')
    parser.add_option('-m', '--min', dest='min', default=1, type=int,
        help='min seqs/sample, default 1')
    parser.add_option('-x', '--max', dest='max', default=10, type=int,
        help='max seqs/sample (inclusive), default 10')
    parser.add_option('-s', '--step', dest='step', default=1, type=int,
        help='levels: min, min+step... for level <= max, default 1')
    parser.add_option('-n', '--num-reps', dest='num_reps', default=1, type=int,
        help='num reps at each seqs/sample level, default 1')
    options, args = parser.parse_args()
    if (len(args) != 1) :
        parser.error("incorrect number of arguments")
    return options, args


if __name__ == '__main__':
    #if called from command-line, should run the analysis.
    options, args = make_cmd_parser()
    otu_path = args[0]
    maker = RarefactionMaker(otu_path, options.min, options.max,
        options.step, options.num_reps)
    maker.rarefy_to_files(options.output_dir)
