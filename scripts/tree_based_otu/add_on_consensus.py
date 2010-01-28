#!/usr/bin/env python

"""Adds on lineage information to an OTU table"""

from sys import argv

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Pre-release"

def parse_taxonomy(lines):
    """output from parse_greengenes_records.py"""
    result = {}
    for line in lines:
        fields = line.strip().split('\t')
        result[fields[1]] = fields[4]
    return result

def parse_otu_file(lines):
    """output from collapse_taxontree_nodes.py"""
    result = {}
    for line in lines:
        otu, taxonid, threshold = line.strip().split('\t')
        if float(threshold) > 0.01:
            break
        result[otu] = taxonid
    return result

def main():
    taxon_lookup = parse_taxonomy(open(argv[1]))
    otu_lookup = parse_otu_file(open(argv[2]))
    
    otu_table = open(argv[3]).readlines()

    new_table = []

    new_table.append(otu_table[0].strip() + '\tNCBI_LINEAGE')

    for line in otu_table[1:]:
        otu_in_table = line.split('\t',1)[0]
        taxonid = otu_lookup[otu_in_table]
        lineage = taxon_lookup.get(taxonid, 'UNCLASSIFIED')
        new_line = '\t'.join([line.strip(), lineage])
        new_table.append(new_line)

    output = open(argv[3] + '-with_lineage.txt','w')
    output.write('\n'.join(new_table))
    output.close()

if __name__ == '__main__':
    main()
