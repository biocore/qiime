#!/usr/bin/env python
#file make_3d_plots.py

__author__ = "Dan Knights"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Dan Knights", "Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Development"

from qiime.parse import parse_otu_table
from numpy import array,apply_along_axis,dot,delete,argsort
import numpy as np

def get_taxa(taxa_fname,sample_ids_kept=None):
    """Opens and returns taxon summaries
       Parameters
        sample_ids, optional list of ids; all other ids are removed

       Returns lineages, counts
    """
    # future: pass in open file object instead
    taxa_f = open(taxa_fname, 'U')

    sample_ids, otu_ids, otu_table, lineages =\
        parse_otu_table(taxa_f,count_map_f=float)
    if sample_ids_kept:
        sam_idxs = [sample_ids.index(sam) for sam in sample_ids_kept]
        otu_table = otu_table[:,sam_idxs]
    return otu_ids, otu_table


def get_taxa_coords(tax_counts,sample_coords):
    """Returns the PCoA coords of each taxon based on the coords of the samples."""    
    # normalize taxa counts along each row/sample (i.e. to get relative abundance)
    tax_counts = apply_along_axis(lambda x: x/float(sum(x)), 0, tax_counts)
    # normalize taxa counts along each column/taxa (i.e. to make PCoA score contributions sum to 1)
    tax_ratios = apply_along_axis(lambda x: x/float(sum(x)), 1, tax_counts)
    return(dot(tax_ratios,sample_coords))

def get_taxa_prevalence(tax_counts):
    """Returns the each lineage's portion of the total count 
    
    takes an otu_table (rows = otus), normalizes samples to equal counts,
    and returns each otu's relative representation in this normalized otu table,
    scaled such that the rarest otu is 0, most prominent is 1
    """
    tax_ratios = apply_along_axis(lambda x: x/float(sum(x)), 0, tax_counts)
    lineage_sums = apply_along_axis(lambda x: sum(x), 1, tax_ratios)
    total_count = sum(lineage_sums)
    prevalence = lineage_sums / float(total_count)
    # scale prevalence from 0 to 1
    prevalence = (prevalence - min(prevalence)) / (max(prevalence) - min(prevalence))
    return prevalence

def remove_rare_taxa(taxdata,nkeep=-1):
    """Keeps only requested number of taxa. Removes empty taxa."""
    if nkeep > 0 and nkeep < len(taxdata['prevalence']):
        ixs = argsort(taxdata['prevalence'])
        ixs = ixs[::-1][:nkeep]
        ixs.sort()
        taxdata['counts'] = taxdata['counts'][ixs,:]
        taxdata['lineages'] = taxdata['lineages'][ixs]
        taxdata['prevalence'] = taxdata['prevalence'][ixs]

    # remove empty taxa
    tax_sums = taxdata['counts'].sum(1)
    for i in xrange(len(tax_sums)-1,-1,-1):
        if tax_sums[i] == 0:
            taxdata['counts'] = delete(taxdata['counts'], i, 0)
            taxdata['lineages'] = delete(taxdata['lineages'], i, 0)
            taxdata['prevalence'] = delete(taxdata['prevalence'], i, 0)


    #Write taxa points and labels if requested
def make_mage_taxa(taxa, num_coords, pct_var, scaled=False, scalars=None,
                   radius=1,
                   min_taxon_radius=0.5, max_taxon_radius=5,
                   taxon_alpha=.5):
    result = []

    # prepare taxonomy radius, coordinates if necessary
    ids = taxa['lineages']
    taxa_coords = taxa['coord']
    if scaled:
        taxa_coords = scale_taxa_data_matrix(taxa_coords,pct_var)
    taxa_radii = radius * (min_taxon_radius+(max_taxon_radius-min_taxon_radius)*taxa['prevalence'])
    radius_dict = dict(zip(ids,taxa_radii))  
#    if scaled:
#        taxa_coords = get_taxa_coords(taxa['tax_counts'],coords)
    coord_dict = dict(zip(ids, taxa_coords))
    result.append('@group {Taxa (n=%s)} collapsible' % (len(ids)))
    color = 'white'
    coord_lines = []
    for id_ in sorted(ids):
        if id_ in coord_dict:
            # note: we display id_[5:] to suppress 'Root;'
            coord_line = '{%s} %s' %(id_[5:], ' '.join(map(str,coord_dict[id_][:num_coords])))
            # each taxon has a different radius, so we have to create a new list
            result.append('@balllist color=%s radius=%s alpha=%s dimension=%s \
master={taxa_points} nobutton' % (color, radius_dict[id_], taxon_alpha, num_coords))
            result.append(coord_line)
            result.append('@labellist color=%s radius=%s alpha=%s dimension=%s \
master={taxa_labels} nobutton' % (color, radius_dict[id_], taxon_alpha, num_coords))
            result.append(coord_line)

    return result


def scale_taxa_data_matrix(coords, pct_var):
    """Scales pc data matrix by percent variation"""
    return coords[:,:len(pct_var)] * (pct_var / pct_var.max())
