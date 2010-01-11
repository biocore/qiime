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
import sys
from cogent.maths.stats.rarefaction import subsample

class SingleRarefactionMaker(FunctionWithParams):
    def __init__(self, otu_path, depth):
        """ otu_path can be path or otu tuple, defined by util.py's getOtuTable
        we just ignore any rarefaction levels beyond any sample in the data
        """
        #self.rare_depths = range(min,max+1, step)
        
        self.depth = depth

        self.sample_names, self.taxon_names, self.otu_table, self.lineages = \
            self.getOtuTable(otu_path)
        self.max_num_taxa = (self.otu_table.sum(1)).max()


    def rarefy_to_file(self, output_fname, small_included=False):
        """ computes rarefied otu tables and writes them, one at a time
        this prevents large memory usage
        
        for depth in self.rare_depths:
            for rep in range(self.num_reps):"""
        sub_sample_ids, sub_otu_ids, sub_otu_table = \
        get_rare_data(self.sample_names, self.taxon_names,
            self.otu_table, self.depth, small_included)
        self._write_rarefaction(output_fname, sub_sample_ids, sub_otu_ids,\
            sub_otu_table)
    
    def _write_rarefaction(self, fname, sub_sample_ids, sub_otu_ids,\
        sub_otu_table):
        """ depth and rep can be numbers or strings
        """
        if min(numpy.shape(sub_otu_table)) == 0: # no data to write
            return
        f = open(fname, 'w')
        f.write(format_otu_table(sub_sample_ids, sub_otu_ids,\
            sub_otu_table, None, comment=fname))
        f.close()


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


    def rarefy_to_files(self, output_dir, small_included=False, 
        include_full=False):
        """ computes rarefied otu tables and writes them, one at a time
        this prevents large memory usage"""
        self.output_dir = output_dir
        for depth in self.rare_depths:
            for rep in range(self.num_reps):
                sub_sample_ids, sub_otu_ids, sub_otu_table = \
                get_rare_data(self.sample_names, self.taxon_names,
                    self.otu_table, depth, small_included)
                self._write_rarefaction(depth, rep, sub_sample_ids, 
                    sub_otu_ids,\
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
    
    include_small_sampes=False => do not write samples with < seqs_per_sample
    total sequecnes
    otu_table (input and out) is otus(rows) by samples (cols)"""
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
    return res_sample_ids, otu_ids, res_otu_table
    
    
usage_str = \
"""usage: %prog [options] {-i INPUT_PATH -o OUTPUT_PATH -m MIN -x MAX -s STEP}
or {-i INPUT_PATH -o OUTPUT_PATH -d DEPTH} (single output file)

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
single output file:
python %prog otu_table.txt -o rarefaction_20_17.txt -d 20
(subsample otu_table.txt w/o replacement at 20 seqs per sample, write resulting
otu table to 'rarefaction_20_17.txt'
or multiple output files:
python %prog otu_table.txt -m 100 -x 1200 -s 100 -n 2 -o mo/rare
(subsample otu_table.txt w/o replacement at 100 seqs per sample (twice),
200 seqs per sample (twice) ... 1200 seqs per sample (twice).
write 24 files total to mo/rare directory

Description:
input_filepath is an otu table, output will be a rarefied otu table
(with single output file syntax), or a series of otu table files, 
named e.g.: rarefaction_4_2.txt
(4 seqs per sample, iteration 2 (3rd such file written)
"""
def parse_command_line_parameters():
    """returns command-line options"""

    if len(sys.argv) == 1:
        sys.argv.append('-h')
    usage = usage_str
    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)
    
    parser.add_option('-i', '--input_path',
        help='input filepath, and otu table [REQUIRED]')
        
    parser.add_option('-o', '--output_path',
        help='write output rarefied otu tables here ' +\
            '(directory if generating multiple outputs, else filename) '+\
            "makes dir if it doesn't exist [REQUIRED]")

    parser.add_option('-m', '--min', type=int,
        help='min seqs/sample')
        
    parser.add_option('-x', '--max', type=int,
        help='max seqs/sample (inclusive)')
        
    parser.add_option('-s', '--step', type=int,
        help='levels: min, min+step... for level <= max')
   
    parser.add_option('-n', '--num-reps', dest='num_reps', default=1, type=int,
        help='num iterations at each seqs/sample level [default: %default]')

    parser.add_option('-d', '--depth', type=int,
        help=\
         'sequences per sample, required when generating a single output file')
    
    parser.add_option('--small_included', dest='small_included', default=False,
        action="store_true",
        help="""samples containing fewer seqs than the rarefaction
level are included in the output but not rarefied [default: %default]""")

    opts, args = parser.parse_args()
        
    if len(args) != 0:
        parser.error("positional argument detected.  make sure all"+\
         ' parameters are identified.' +\
         '\ne.g.: include the \"-m\" in \"-m MINIMUM_LENGTH\"')
         
    required_options1 = ['input_path', 'output_path', 'min', 'max', 'step']
    required_options2 = ['input_path', 'output_path', 'depth']
    if opts.depth != None:
        for option in required_options2:
            if eval('opts.%s' % option) == None:
                parser.error('Required option --%s omitted.' % option)
        if opts.min or opts.max or opts.step:
            parser.error('specifying option "depth" precludes use of min, '+\
                'max, or step, as only one output file is generated')
    else:
        for option in required_options1:
            if eval('opts.%s' % option) == None:
                parser.error('Required option --%s omitted.' % option)
    return opts, args


if __name__ == '__main__':
    #if called from command-line, should run the analysis.
    opts, args = parse_command_line_parameters()
    otu_path = opts.input_path
    if opts.depth != None:
        
        maker = SingleRarefactionMaker(otu_path, opts.depth)
        maker.rarefy_to_file(opts.output_path, opts.small_included)
    else:
        if not os.path.exists(opts.output_path):
            os.makedirs(opts.output_path)
        maker = RarefactionMaker(otu_path, opts.min, opts.max,
            opts.step, opts.num_reps)
        maker.rarefy_to_files(opts.output_path, opts.small_included)
