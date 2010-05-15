#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight", "Justin Kuczynski","Jeremy Widmann"] 
#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

import numpy
from numpy import isnan
from StringIO import StringIO
from cogent import Sequence

"""Contains formatters for the files we expect to encounter in 454 workflow.

A lot of this might migrate into cogent at some point.
"""

def format_distance_matrix(labels, data):
    """Writes distance matrix as tab-delimited text, uses format_matrix"""
    return format_matrix(data, labels, labels)

def format_matrix(data, row_names, col_names):
    """Writes matrix as tab-delimited text.

    format is rows: samples, 
    cols: metrics/ confidence levels, etc

    data: array or 2d list.
    """
    len_col = len(col_names)
    try:
        if data.shape != (len(row_names), len_col):
            raise ValueError, "Data shape of %s doesn't match header sizes %s %s" %\
                (data.shape, len(row_names), len(col_names))
    except AttributeError:
        # must be list of list
        try:
            if not numpy.all([len_col==len(row) for row in data]) or\
                    len(row_names) != len(data):
                raise ValueError, "Data shape doesn't match header sizes %s %s" %\
                    (len(row_names), len(col_names))
        except:
            raise ValueError, "Unsupported data type for format_matrix"

    lines = []
    row_names = map(str, row_names)   
    col_names = map(str, col_names)   
    #just in case they weren't strings initially
    lines.append('\t'.join([''] + col_names))
    for sam, vals in zip(row_names, data):
        lines.append('\t'.join([sam] + map(str, vals)))
    return '\n'.join(lines)

def format_otu_table(sample_names, otu_names, data, taxonomy=None,
    comment='Full OTU Counts', skip_empty=False):
    """Writes OTU table as tab-delimited text.
    
    inputs: sample_names, otu_names are lists of strings
    data is numpy 2d array, num_otus x num_samples
    taxonomy is list of length = num_otus
    
    """
    if data.shape != (len(otu_names), len(sample_names)):
        raise ValueError, "Data shape of %s doesn't match %s OTUs, %s samples" \
            % (data.shape, len(otu_names), len(sample_names))
    lines = []
    #data = numpy.array(data, dtype='str') ##BAD! truncates some ints!
    sample_names = map(str, sample_names)
    otu_names = map(str, otu_names)
    lines.append('#'+comment)
    if taxonomy:
        lines.append('\t'.join(['#OTU ID'] + sample_names + 
            ['Consensus Lineage']))
        for otu_name, vals, taxon in zip(otu_names, data, taxonomy):
            if (skip_empty and filter(lambda a: a!=0, vals)==[]):
                #skip otu with zero counts
                continue
            if not isinstance(taxon, str):
                taxon = ';'.join(taxon)
            lines.append('\t'.join([otu_name] + map(str, vals.tolist()) + [taxon]))
    else:
        lines.append('\t'.join(['#OTU ID'] + sample_names))
        for otu_name, vals in zip(otu_names, data):
            lines.append('\t'.join([otu_name] + map(str,vals.tolist())))
    return '\n'.join(lines)

def format_coords(coord_header, coords, eigvals, pct_var):
    """formats coords given specified coords matrix etc."""
    result = []
    result.append('pc vector number\t' +
        '\t'.join(map(str, range(1,len(coord_header)+1))))
    for name, row in zip(coord_header, coords):
        result.append('\t'.join([name] + map(str, row)))
    result.append('')
    result.append('')
    result.append('eigvals\t' + '\t'.join(map(str,eigvals)))
    result.append('% variation explained\t' +
        '\t'.join(map(str, pct_var)))
    return '\n'.join(result)

def build_prefs_string(mapping_headers_to_use, background_color, monte_carlo_dist, headers, otu_ids):
    """Create a preferences file, which can be used for some of the \
    visualization scripts."""
    
    #Open up the prefs dictionary
    pref_lines=["{\n"]
    
    #Define and add the background_color dictionary to prefs dictionary
    bk_color="'background_color':'%s',\n" % (background_color)
    pref_lines.append(bk_color)
    
    #create a unique field dictionary for use with the FIELDS dictionary
    unique_dist_fields={}  
    
    #Iterate through the user-supplied fields or all the fields in the mapping
    #file, then validate that the fields exist
    if mapping_headers_to_use=='ALL':
        fields=headers
        for key in fields:
            unique_dist_fields[key]=monte_carlo_dist
    else:
        #If '&&' is used, split into multiple fields and valid them against
        #the mapping file
        fields=mapping_headers_to_use.split(',')
        for field_id in fields:
            f_str=field_id.split('&&')
            for id_ in f_str:
                validity=False
                for head_id in headers:
                    if id_==head_id:
                        unique_dist_fields[id_]=monte_carlo_dist
                        validity=True
                if not validity:
                    raise ValueError, \
                            "%s is not a header in your mapping file" % field_id
            
    #Syntax for sample_coloring dictionary
    sample_coloring = ["\n'sample_coloring':\n\t{"]
    sample_colors = \
    "\t\t'%s':" + \
    "\n\t\t{"+ \
    "\n\t\t\t'column':'%s',"+ \
    "\n\t\t\t'colors':(('red',(0,100,100)),('blue',(240,100,100)))"+ \
    "\n\t\t}"
    
    #Syntax for monte_carlo dictionary
    monte_carlo_main = ["\n'MONTE_CARLO_GROUP_DISTANCES':\n\t{\n"]
    monte_carlo=[]
    monte_carlo_distances = "\t\t'%s': %s"
    
    #Syntax for fields dictionary
    field_dict_main = ["\n'FIELDS':\n\t[\n"]
    field_dict=[]
    dist_fields = "\t\t'%s'"
    
    #This iterates through the fields and creates a sample_color dictionary     
    #values for each field
    first = True
    for field in fields:
        if first:
            first=False
            sample_coloring.append('\n')
        else:
            sample_coloring.append(',\n')
        sample_coloring.append(sample_colors % (field, field))
        monte_carlo.append(monte_carlo_distances % (field,monte_carlo_dist))
    
    #Syntax for taxonomy_coloring dictionary
    taxon_start = \
    "\n'taxonomy_coloring':\n\t{\n" + \
    "\t\t'Level_%s':" + \
    "\n\t\t{"+ \
    "\n\t\t\t'column':'%s',"+ \
    "\n\t\t\t'colors':\n\t\t\t{"
    taxon_colors = "\n\t\t\t\t'%s':('red%s',(%d,100,100))"
    taxon_coloring=[]
    
    if otu_ids:
        level = max([len(t.split(';')) - 1 for t in otu_ids])
        taxon_coloring.append(taxon_start % (str(level),str(level)))
        taxons=[]
        otu_id_iter=(240/(len(otu_ids)+1))
        counter=0    
        for i in otu_ids:
            taxons.append(taxon_colors % (i,str(counter),counter))
            counter=counter+otu_id_iter
        taxon_coloring.append(','.join(taxons))
    else:
        taxon_coloring.append(taxon_start % (str(1),str(1)))
        taxon_coloring.append(taxon_colors % ('Root;Bacteria',str(0),0))
    
    taxon_coloring.append("\n\t\t\t}\n\t\t}\n\t}")
    taxonomy_coloring_str=''.join(taxon_coloring)
    
    #Close and convert the sample_coloring dictionary to a string
    sample_coloring.append('\n\t},')
    sample_coloring_str=''.join(sample_coloring)
    
    #Close and convert the monte_carlo dictionary to a string
    monte_carlo1=',\n'.join(monte_carlo)
    monte_carlo_main.append(monte_carlo1)
    monte_carlo_main.append('\n\t},')
    monte_carlo_str=''.join(monte_carlo_main)
    
    #This iterates through the fields and creates the monte_carlo and fields
    #dictionary values
    for field in unique_dist_fields:
        field_dict.append(dist_fields % (field))
    
    #Close and convert the fields dictionary to a string
    field_dict1=',\n'.join(field_dict)
    field_dict_main.append(field_dict1)
    field_dict_main.append('\n\t],')
    field_dict_str=''.join(field_dict_main)
    
    #Add all the dictionary values to the prefs dictionary
    pref_lines.append(sample_coloring_str)
    pref_lines.append(monte_carlo_str)
    pref_lines.append(field_dict_str)
    pref_lines.append(taxonomy_coloring_str)
    pref_lines.append('\n}')

    return ''.join(pref_lines)

def format_map_file(headers, id_map, desc_key, sample_id_key, \
    description_map=None, run_description=None):
    """Generates string for formatted map file.
    
    Input:
        headers: list of strings corresponding to col headers
        id_map: dict of {id:{header:val}}
        description_map: dict of {id:description}
        run_description: either string, or list of strings
    """
    result = []
    if desc_key in headers:
        headers.remove(desc_key)
    if sample_id_key in headers:
        headers.remove(sample_id_key)
    header_line = '\t'.join([sample_id_key] + headers + [desc_key])
    if not header_line.startswith('#'):
        header_line = '#' + header_line
    result.append(header_line)
    if run_description:
        if not isinstance(run_description, str):
            run_description = '\n#'.join(run_description)
        if not run_description.startswith('#'):
            run_description = '#'+run_description
        result.append(run_description)
    for id_, fields in sorted(id_map.items()):
        curr_line = [id_]
        curr_line.extend([fields.get(h,'') for h in headers])
        curr_line.append(description_map.get(id_,''))
        result.append('\t'.join(map(str, curr_line)))
    return '\n'.join(result)
    
def format_histograms(pre_hist, post_hist, bin_edges):
    """Returns text-formatted histogram."""
    lines = []
    lines.append('Length\tBefore\tAfter')
    for edge, pre, post in zip(bin_edges, pre_hist, post_hist):
        lines.append('\t'.join(map(str, [edge, pre, post])))
    return '\n'.join(lines)

def format_unifrac_sample_mapping(sample_ids, otu_ids, otu_table_array):
    """Returns a unifrac sample mapping file from output of parse_otu_table
    """
    out = []
    for i, row in enumerate(otu_table_array):
        for j, val in enumerate(row):
            if val > 0:
                line = [otu_ids[i], sample_ids[j], str(val)]
                out.append('\t'.join(line))
    return out

def write_Fasta_from_name_seq_pairs(name_seqs, fh):
    """writes a list of (name,seqs) to filehandle.

    name_seqs: (name,seqs) pair such as from MinimalFASTAParser
    fh: an open filehandle
    """
    if fh==None:
        raise ValueError,"Need open file handle to write to." 

    for (name,seq) in name_seqs:
        fh.write("%s\n"% Sequence(name=name, seq = seq).toFasta())
    

