#!/usr/bin/env python

"""Converts BLAST results into an OTU table"""

from numpy import zeros
from os import listdir, path
from optparse import OptionParser, make_option
from operator import itemgetter
from cogent.parse.blast import BlastResult
from cogent.util.misc import unzip
from collections import defaultdict

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Pre-release"

options = [make_option('--blast-results-path',dest='blast_res_path',\
                       default=None),
           make_option('--sequence-to-otu-map',dest='seq_to_otu', \
                       default=None),
           make_option('--sequence-to-db-map',dest='seq_to_db',\
                       default=None)]

def parse_seq_to_otu(lines):
    """Mapping from representative_sequences"""
    result = {}
    for line in lines:
        fields = line.strip().split()
        result[fields[0]] = fields[1]
    return result

def parse_seq_to_db(lines):
    """Mapping from input sequences to microbiome database"""
    result = {}
    for line in lines:
        fields = line.strip().split()
        result[fields[0]] = '_'.join([fields[1], fields[2]])
    return result

def best_by_evalue(blast_result):
    """Get best hit by evalue"""
    vals = []
    for d in blast_result:
        vals.append((d['E-VALUE'], d))
    return sorted(vals, key=itemgetter(0))[0]

clean_queryid = lambda x: x[1:] # input sequences had an extra >

def get_best_hit(record):
    """Returns best hit sequence id and subject id"""
    queryid, hits = record
    evalue, best_hit = best_by_evalue(hits)
    return (clean_queryid(queryid), best_hit['SUBJECT ID'])

def convert_to_matrix(otu_table):
    """Build otu matrix"""
    all_sampleids, all_otus = unzip(otu_table.keys())
    all_sampleids = sorted(set(all_sampleids))
    all_otus = sorted(set(all_otus))
    matrix = zeros((len(all_otus), len(all_sampleids)), int)
    for row, otu in enumerate(all_otus):
        for col, sampleid in enumerate(all_sampleids):
            matrix[row, col] += otu_table.get((sampleid, otu), 0)
    return matrix, all_otus, all_sampleids

def main():
    parser = OptionParser(option_list=options)
    opts, args = parser.parse_args()

    seq_to_db = parse_seq_to_db(open(opts.seq_to_db))
    seq_to_otu = parse_seq_to_otu(open(opts.seq_to_otu))

    otu_table = defaultdict(int)

    for res_file in listdir(opts.blast_res_path):
        abs_res_file = path.join(opts.blast_res_path, res_file)
        results = BlastResult(open(abs_res_file))
        print abs_res_file 
        for hits in results.iterHitsByQuery():
            queryid, subjectid = get_best_hit(hits)
            sampleid = seq_to_db[queryid]
            otu = seq_to_otu[subjectid]
            otu_table[(sampleid, otu)] += 1

    matrix, otu_rows, sampleid_cols = convert_to_matrix(otu_table)
    output = open(opts.blast_res_path + '.otutable','w')
    header = '# OTU\t%s\n' % '\t'.join(sampleid_cols)
    otu_counts = '%s\t%s\n'
    output.write(header)
    for row_num, row in enumerate(matrix):
        output.write(otu_counts % (otu_rows[row_num], '\t'.join(map(str,row))))
    output.close()

if __name__ == '__main__':
    main()
