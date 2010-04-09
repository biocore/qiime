#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

"""Contains code for adding taxa to OTU table that lacks them.
"""
from sys import stdout
from string import strip
from qiime.parse import fields_to_dict

def fix_taxonomy_delimiters(taxonomy):
    """fixes delimiters in taxonomy (expect semicolons, but get commas)"""
    result = {}
    for k, vals in taxonomy.iteritems():
        v = vals[0]
        if ';' in v:
            result[k] = v.replace('"','')
        else:
            result[k] = v.replace(',',';').replace('"','')
    return result


def rewrite_otu_table_with_taxonomy(taxon_lines, otu_lines, id_map_lines=None,
    outfile=stdout):
    """Rewrites OTU table including taxonomy."""
    taxonomy = fields_to_dict(taxon_lines)
    #sometimes have extra fields after OTU id
    new_taxonomy = {}
    for k, v in taxonomy.items():
        new_taxonomy[k.split()[0]] = v
    taxonomy = new_taxonomy
    taxonomy = fix_taxonomy_delimiters(taxonomy)

    if id_map_lines:
        id_map = dict([map(strip, line.split('\t')) for line in
            id_map_lines])
        new_taxonomy = dict([(id_map[k], v) for k, v in taxonomy.items()
            if k in id_map])
        assert new_taxonomy != taxonomy
        taxonomy = new_taxonomy

    for line in otu_lines:
        if not line.endswith('\n'):
            line += '\n'
        if line.startswith('#OTU ID'):
            outfile.write(line[:-1]+'\tConsensus Lineage\n')
        elif line.startswith('#'):
            outfile.write(line)
        else:
            id_, rest = line.split('\t', 1)
            t = taxonomy.get(id_, 'None')
            outfile.write(line[:-1]+'\t'+t+'\n')

