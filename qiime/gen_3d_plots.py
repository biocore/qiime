#!/usr/bin/env python
#file gen_3d_plots.py

__author__ = "Jesse Stombaugh and Micah Hamady"
__copyright__ = "Copyright 2009, the 454 Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Prototype"

"""
Author: Jesse Stombaugh (jesse.stombaugh@colorado.edu) and Micah Hamady
Status: Prototype

Requirements:
Python 2.5

Example 1: Create 3D plot from only the pca/pcoa data, where each ID is colored:
Usage: python gen_3d_plots.py -c raw_pca_data.txt

Example 2: Create a Kinemage with two coloring schemes (Day and Type):
Usage: python gen_3d_plots.py -c raw_pca_data.txt -m input_map.txt -b 'Day,Type'

Example 3: Create 3D plots for a combination of label headers from a mapping 
file:
Usage: python gen_3d_plots.py -c raw_pca_data.txt -m input_map.txt 
-b 'Type&&Day' -x ./test/

"""

from cogent.util.misc import flatten
from parse import parse_map, group_by_field, group_by_fields, parse_coords
from numpy import array
import os
from optparse import OptionParser
from random import choice, randrange

data_colors = ['blue','lime','red','aqua','fuchsia','yellow','green', \
              'maroon','teal','purple','olive','white','silver','gray']

class MissingFileError(IOError):
    pass

def make_3d_plots(coord_header, coords, pct_var, mapping, prefs):
    """Makes 3d plots given coords, mapping file, and prefs."""
    
    result = []
    #Iterate through prefs and color by given mapping labels
    for name, p in prefs.items():
        group_num=-1
        col_name = p['column']
        if 'colors' in p:
            colors = p['colors'].copy()    #copy so we can mutate
        else:
            colors={}
        labelname=prefs[name]['column']
        
        #Define groups and associate appropriate colors to each group
        groups = group_by_field(mapping, col_name)
        for g in groups:
            if g not in colors:
                group_num+=1
                if group_num==len(data_colors):
                    group_num=0
                colors[g] = data_colors[group_num]
        
        #Write to kinemage file using the groups, colors and coords 
        result.extend(make_mage_output(groups, colors, coord_header, coords, \
            pct_var, labelname,scaled=False))
        result.extend(make_mage_output(groups, colors, coord_header, coords,  \
            pct_var, labelname,scaled=True))

    return result

def scale_pc_data_matrix(coords, pct_var):
    """Scales pc data matrix by percent variation"""

    return coords * (pct_var / pct_var.max())

def auto_radius(coords,ratio=0.01):
    """Determine radius from coords"""
    dim1 = coords[:,0]
    range = dim1.max()-dim1.min()
    
    return ratio*range

def make_mage_output(groups, colors, coord_header, coords, pct_var, name='', \
    radius=None, alpha=.75, num_coords=10,scaled=False, coord_scale=1.05):
    """Convert groups, colors, coords and percent var into mage format"""
    result = []
    
    #Scale the coords and generate header labels
    if scaled:
        coords = scale_pc_data_matrix(coords, pct_var)
        header_suffix = '_scaled'
    else:
        header_suffix = '_unscaled'
        
    if radius is None:
        radius = auto_radius(coords)
        
    maxes = coords.max(0)[:num_coords]
    mins = coords.min(0)[:num_coords]
    pct_var = pct_var[:num_coords]    #scale from fraction
    
    #check that we didn't get fewer dimensions than we wanted
    if len(mins) < num_coords:
        num_coords = len(mins)
    min_maxes = flatten(zip(mins,maxes))
    
    #Write the header information
    result.append('@kinemage {%s}' % (name+header_suffix))
    result.append('@dimension '+' '.join(['{PC_'+str(i+1)+'}' \
        for i in range(num_coords)]))
    result.append('@dimminmax '+ ' '.join(map(str, min_maxes)))
    result.append('@master {points}')
    result.append('@master {labels}')
    result.append('@hsvcolor {aqua} 180 100 100')
    result.append('@hsvcolor {blue} 240 100 100')
    result.append('@hsvcolor {fuchsia} 300 100 100')
    result.append('@hsvcolor {gray} 300 0 50.2')
    result.append('@hsvcolor {green} 120 100 50.2')
    result.append('@hsvcolor {lime} 120 100 100')
    result.append('@hsvcolor {maroon} 0 100 50.2')
    result.append('@hsvcolor {olive} 60 100 50.2')
    result.append('@hsvcolor {purple} 300 100 50.2')
    result.append('@hsvcolor {red} 0 100 100')
    result.append('@hsvcolor {silver} 0 0 75.3')
    result.append('@hsvcolor {teal} 180 100 50.2')
    result.append('@hsvcolor {white} 180 0 100')
    result.append('@hsvcolor {yellow} 60 100 100')              
    
    #Write the groups, colors and coords
    coord_dict = dict(zip(coord_header, coords))
    for group_name, ids in sorted(groups.items()):
        result.append('@group {%s (n=%s)} collapsible' % (group_name, len(ids)))
        color = colors[group_name]
        coord_lines = []
        for id_ in sorted(ids):
            if id_ in coord_dict:
                coord_lines.append('{%s} %s' % \
                    (id_, ' '.join(map(str, coord_dict[id_][:num_coords]))))

        result.append('@balllist color=%s radius=%s alpha=%s dimension=%s \
master={points} nobutton' % (color, radius, alpha, num_coords))
        result.append('\n'.join(coord_lines))
        result.append('@labellist color=%s radius=%s alpha=%s dimension=%s \
master={labels} nobutton' % (color, radius, alpha, num_coords))
        result.append('\n'.join(coord_lines))

    #Write the axes on the bottom of the graph
    result.append('@group {axes} collapsible')
    state = 'on'
    axis_mins = mins*coord_scale
    axis_maxes = maxes*coord_scale
    for i, pct in enumerate(pct_var):
        if i == 3:
            state = 'off'
        result.append('@vectorlist {PC%s line} dimension=%s %s' % \
            (i+1, num_coords, state))
        result.append(' '.join(map(str, axis_mins)) + ' white')
        end = axis_mins.copy()
        end[i] = axis_maxes[i]
        result.append(' '.join(map(str, end)) + ' white')
        result.append('@labellist {PC%s (%0.2g%%)} dimension=%s %s' % \
            (i+1, pct, num_coords, state))
        end[i] *= coord_scale  #add scale factor to offset labels a little
        result.append( ('{PC%s (%0.2g%%)}' % (i+1, pct))  + \
            ' '.join(map(str, end)) + ' white')
            
    return result

def combine_map_label_cols(combinecolorby,mapping):
    """Merge two or more mapping columns into one column"""
    combinedmapdata=array(['']*len(mapping),dtype='a100')
    for p in range(len(combinecolorby)):                    
        for i in range(len(mapping[0])):
            if str(combinecolorby[p])==str(mapping[0][i]):
                for q in range(len(mapping)):
                    combinedmapdata[q]=combinedmapdata[q]+mapping[q][i]
    for i in range(len(combinedmapdata)):
        mapping[i].append(combinedmapdata[i])
    return mapping

def create_dir(dir_path):
    """Creates directory where data is stored.  If directory is not supplied in \
       the command line, a random folder is generated"""
       
    alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
    alphabet += alphabet.lower()
    alphabet += "01234567890"

    if dir_path==None or dir_path=='':
        random_dir_name=''.join([choice(alphabet) for i in range(10)])
        dir_path = './'+random_dir_name+'/'

    if dir_path:
        try:
            os.mkdir(dir_path)
        except OSError:
            pass

    return dir_path

def process_colorby(colorby,data):
    """Parses the colorby option from the command line"""
    prefs={}
    mapping=data['map']
    colorbydata = colorby.strip().strip("'").split(',')

    for j in range(len(colorbydata)):
        if '&&' in colorbydata[j]:
            #Create an array using multiple columns from mapping file
            combinecolorby=colorbydata[j].split('&&')
            data['map']=combine_map_label_cols(combinecolorby,mapping)
            prefs[str(j)]={}
            prefs[str(j)]['column']=''.join(combinecolorby)
        else:   
            #Color by only one column in mapping file      
            prefs[str(j)]={}
            prefs[str(j)]['column']=colorbydata[j]

    return prefs,data

def _make_path(paths):
    """Join together the paths (e.g. dir and subdir prefix), empty str default"""
    curr = ''
    for p in paths:
        if p:
            curr += p
            if curr and (not curr.endswith('/')):
                curr += '/'
    return curr

#The following functions were not unit_tested, however the parts within
#the functions are unit_tested
def get_map(options, data):
    """Opens and returns mapping data"""
    try:
        map_f = open(options.map_fname, 'U').readlines()
    except (TypeError, IOError):
        raise MissingFileError, 'Mapping file required for this analysis'
    data['map'] = parse_map(map_f)
    return data['map']

def get_coord(options, data):
    """Opens and returns coords data"""
    try:
        coord_f = open(options.coord_fname, 'U').readlines()
    except (TypeError, IOError):
        raise MissingFileError, 'Coord file required for this analysis'
    coord_header, coords, eigvals, pct_var = parse_coords(coord_f)
    data['coord'] = coord_header, coords, eigvals, pct_var
    return data['coord']

def _do_3d_plots(prefs, data, dir_path='', filename=None, \
                default_filename='out'):
    """Make 3d plots according to coloring options in prefs."""
    kinpath = _make_path([dir_path, filename])
    htmlpath = './'
    if kinpath:
        if filename:
            kinpath = kinpath[:-1]
        else:
            kinpath += default_filename
    if not kinpath:
        return
    coord_header, coords, eigvals, pct_var = data['coord']
    mapping=data['map']
    outf=kinpath+'.kin'
    
    res = make_3d_plots(coord_header, coords, pct_var,mapping,prefs)
    
    #Write kinemage file
    f = open(outf, 'w')
    f.write('\n'.join(res))
    f.close()
    
    #Write html page with the kinemage embedded
    f2 = open(htmlpath+filename+'_3D.html', 'w')
    f2.write("<html><head></head><body><applet code='king/Kinglet.class' \
archive='./king.jar' width=800 height=600> \
<param name='kinSource' value='%s'></body></html>" % (outf)) 
    f2.write('\n'.join(res))
    f2.close()
    
def _make_cmd_parser():
    """Returns the command-line options"""
    parser = OptionParser(usage="Usage: this_file.py -c <pca/pcoa output files>\
\nor\nUsage: this_file.py -c <pca/pcoa output files> -m <mapping output file>\
-b 'Mapping column to color by' -x <write to directory>")
    parser.add_option('-m', '--map', dest='map_fname', \
        help='name of mapping file')
    parser.add_option('-c', '--coord', dest='coord_fname', \
        help='name of coords file')
    parser.add_option('-x', '--dir-prefix', dest='dir_path',\
        help='directory prefix for all analyses')
    parser.add_option('-b', '--by', dest='colorby',\
        help='map header to color by')
    options, args = parser.parse_args()
    return options

def _process_prefs(options):
    """Opens files as necessary based on prefs"""
    data = {}
    dir_path = create_dir(options.dir_path)
    
    #Open and get coord data
    data['coord'] = get_coord(options, data)
    
    #Open and get mapping data, if none supplied create a pseudo mapping \
    #file
    if options.map_fname:
        data['map'] = get_map(options, data)
    else:
        data['map']=(([['#SampleID','Sample']]))
        for i in range(len(data['coord'][0])):
            data['map'].append([data['coord'][0][i],'Sample'])
    
    #Determine which mapping headers to color by, if none given, color by \
    #Sample ID's
    if options.colorby:
        prefs,data=process_colorby(options.colorby,data)
    else:
        prefs={}
        prefs['Sample']={}
        prefs['Sample']['column']='#SampleID'

    filepath=options.coord_fname
    filename=filepath.strip().split('/')[-1]

    action_str = '_do_3d_plots'
    try:
        action = eval(action_str)
    except NameError:
        action = None
    #Place this outside try/except so we don't mask NameError in action
    if action:
        action(prefs, data, dir_path,filename)

if __name__ == '__main__':
    from sys import argv, exit
    options = _make_cmd_parser()
    
    #Kept, just in case we allow for reading a prefs file
    #prefs = eval(open(options.pref_fname, 'U').read())

    _process_prefs(options)

