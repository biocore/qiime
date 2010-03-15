#!/usr/bin/env python
#file make_3d_plots.py

__author__ = "Jesse Stombaugh, Rob Knight, and Dan Knights"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Jesse Stombaugh", "Rob Knight", "Micah Hamady", "Dan Knights"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Pre-release"


from cogent.util.misc import flatten
from qiime.parse import parse_coords,parse_map,group_by_field
from qiime.colors import (natsort, data_color_order, data_colors,
        get_group_colors, color_groups, make_color_dict, combine_map_label_cols,
        process_colorby, linear_gradient,iter_color_groups,get_map)
from numpy import array, shape, apply_along_axis, dot, delete, vstack
import numpy as np
import os
from random import choice
import re
from time import strftime

'''
xdata_colors = {
        'aqua':     (180, 100, 100),
        'blue':     (240,100,100),
        'fuchsia':  (300,100,100),
        'gray':     (300,0,50.2),
        'green':    (120,100,50.2),
        'lime':     (120,100,100),
        'maroon':   (0,100,50.2),
        'olive':    (60,100,50.2),
        'purple':   (300,100,50.2),
        'red':      (0,100,100),
        'silver':   (0, 0, 75.3),
        'teal':     (180,100,50.2),
        'yellow':   (60,100,100),
}
'''
kinemage_colors = ['hotpink','blue', 'lime','gold', \
                       'red','sea','purple','green']

class MissingFileError(IOError):
    pass


def make_3d_plots(coord_header, coords, pct_var, mapping, prefs, \
                    background_color,label_color,custom_axes=None, \
                    data_colors=data_colors,data_color_order=data_color_order, \
                    edges=None):
    """Makes 3d plots given coords, mapping file, and prefs.
    
    Added quick-and-dirty hack for gradient coloring of columns, should
    replace with more general method. Current solution is to pass in any
    of the following:

    'colors':(('white', (0,100,100)),('red',(100,100,100)))

    makes gradient between white and red, applies to all samples

    'colors':{'RK':(('white',(0,0,100)),('red',(0,100,100))),
              'NF':(('white',(120,0,100)),('green',(120,100,100)))
             }
    pulls the combination samples starting with RK, colors with
    first gradient, then pulls the combination samples starting
    with NF, colors with the next gradient.
    """
    
    result = []
    #Iterate through prefs and color by given mapping labels
    #Sort by the column name first
    groups_and_colors=iter_color_groups(mapping,prefs,data_colors, \
                                        data_color_order)
    groups_and_colors=list(groups_and_colors)

    for i in range(len(groups_and_colors)):  
        #Write to kinemage file using the groups, colors and coords 
        labelname=groups_and_colors[i][0]
        groups=groups_and_colors[i][1]
        colors=groups_and_colors[i][2]
        data_colors=groups_and_colors[i][3]
        data_color_order=groups_and_colors[i][4]
        
        result.extend(make_mage_output(groups, colors, coord_header, coords, \
            pct_var,background_color,label_color,custom_axes,name=labelname, \
            scaled=False, data_colors=data_colors, edges=edges))
        result.extend(make_mage_output(groups, colors, coord_header, coords, \
            pct_var,background_color,label_color,custom_axes,name=labelname, \
            scaled=True, data_colors=data_colors, edges=edges))

    return result

def scale_pc_data_matrix(coords, pct_var):
    """Scales pc data matrix by percent variation"""

    return coords * (pct_var / pct_var.max())

def auto_radius(coords,ratio=0.01):
    """Determine radius from coords"""
    dim1 = coords[:,0]
    range = dim1.max()-dim1.min()
    
    return ratio*range

def make_mage_output(groups, colors, coord_header, coords, pct_var, \
                     background_color,label_color,custom_axes=None,name='', \
                     radius=None, alpha=.75, num_coords=10,scaled=False, \
                     coord_scale=1.05, data_colors=data_colors, edges=None):
    """Convert groups, colors, coords and percent var into mage format"""
    result = []
    
    #Scale the coords and generate header labels
    if scaled:
        scalars = pct_var
        if custom_axes:
            custom_scalars = scalars[0] * np.ones(len(custom_axes))
            scalars = np.append(custom_scalars,scalars)
        coords = scale_pc_data_matrix(coords, scalars)
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
    
    if custom_axes:
        axis_names = ['PC%s' %(i+1) for i in xrange(num_coords - len(custom_axes))]
        axis_names = custom_axes + axis_names
    else:
        axis_names = ['PC%s' %(i+1) for i in xrange(num_coords)]

    #Write the header information
    result.append('@kinemage {%s}' % (name+header_suffix))
    result.append('@dimension '+' '.join(['{%s}'%(name) for name in axis_names]))
    result.append('@dimminmax '+ ' '.join(map(str, min_maxes)))
    result.append('@master {points}')
    result.append('@master {labels}')
            
    for name, color in sorted(data_colors.items()):
        result.append(color.toMage())

    if background_color=='white':
        result.append('@whitebackground')
        result.append('@hsvcolor {black} 0.0 0.0 0.0')
    else:
        result.append('@hsvcolor {white} 180.0 0.0 100.0')
   
    #Write the groups, colors and coords
    coord_dict = dict(zip(coord_header, coords))
    for group_name in natsort(groups):
        ids = groups[group_name]
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

    #Write any custom axes first
    result.append('@group {axes} collapsible')
    state = 'on'
    axis_mins = mins*coord_scale
    axis_maxes = maxes*coord_scale

    #Write the axes on the bottom of the graph
    if not custom_axes:
        custom_axes = []
    for i in xrange(num_coords):
        if i == 3:
            state = 'off'            
        result.append('@vectorlist {%s line} dimension=%s %s' % \
            (axis_names[i], num_coords, state))
            
        result.append(' '.join(map(str, axis_mins)) + ' ' + label_color)
        end = axis_mins.copy()
        end[i] = axis_maxes[i]
        result.append(' '.join(map(str, end)) + ' ' + label_color)
        end[i] *= coord_scale  #add scale factor to offset labels a little
            
        if i < len(custom_axes):
            result.append('@labellist {%s} dimension=%s %s' % \
                              (axis_names[i], num_coords, state)) 
            result.append( ('{%s}' % (axis_names[i]))  + \
                               ' '.join(map(str, end)) + ' ' + label_color)
        else:
            pct = pct_var[i-len(custom_axes)]
            result.append('@labellist {%s (%0.2g%%)} dimension=%s %s' % \
                              (axis_names[i], pct, num_coords, state))
            result.append( ('{%s (%0.2g%%)}' % (axis_names[i], pct))  + \
                               ' '.join(map(str, end)) + ' ' + label_color)

    #Write edges if requested
    if edges:
        result += make_edges_output(coord_dict, edges, num_coords,label_color)
    return result

def make_edges_output(coord_dict, edges, num_coords,label_color):
    """creates make output to display edges (as a kinemage vectorlist)"""
    result = []
    result.append('@vectorlist {edges} dimension=%s on' % \
                      (num_coords))
    for edge in edges:
        id_fr, id_to = edge
        # get 'index' of the destination set from 'to' sampleID
        which_set = int(id_to[id_to.rindex('_')+1:]) - 1
        which_color = kinemage_colors[which_set % len(kinemage_colors)]
        # plot a color 'tip' on the line (10% of line length)
        pt_fr = coord_dict[id_fr][:num_coords]
        pt_to = coord_dict[id_to][:num_coords]
        diffs = (pt_to-pt_fr) * .66
        middles = pt_fr + diffs
        result.append('%s %s' % \
                          (' '.join(map(str, pt_fr)),label_color))
        result.append('%s %s P' % \
                          (' '.join(map(str, middles)),label_color))
        
        result.append('%s %s' % \
                          (' '.join(map(str, middles)), which_color))
        result.append('%s %s P' % \
                          (' '.join(map(str, pt_to)), which_color))            
    return result

def create_dir(dir_path,plot_type):
    """Creates directory where data is stored.  If directory is not supplied in\
       the command line, a random folder is generated"""
       
    alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
    alphabet += alphabet.lower()
    alphabet += "01234567890"

    
    if dir_path==None or dir_path=='':
        dir_path=''
        random_dir_name=''.join([choice(alphabet) for i in range(10)])
        dir_path ='./'+plot_type+strftime("%Y_%m_%d_%H_%M_%S")+random_dir_name+'/'

    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    return dir_path

def process_custom_axes(axis_names):
    """Parses the custom_axes option from the command line"""
    return axis_names.strip().strip("'").strip('"').split(',')

def process_coord_filenames(coord_filenames):
    """Parses the custom_axes option from the command line"""
    return coord_filenames.strip().strip("'").strip('"').split(',')

def get_custom_coords(axis_names,mapping, coords):
    """Gets custom axis coords from the mapping file."""
    for i, axis in enumerate(reversed(axis_names)):
        if not axis in mapping[0]:
            print 'Warning: could not find custom axis',axis,'in map headers:',mapping[0]
        else:
            col_idx = mapping[0].index(axis)
            col = zip(*mapping[1:])[col_idx]
            sample_IDs = zip(*mapping[1:])[0]
            new_coords = array([])
            for id in coords[0]:
                if id in sample_IDs:
                    row_idx = list(sample_IDs).index(id)
                    try:
                        as_float = float(col[row_idx])
                        new_coords = np.append(new_coords,as_float)
                    except ValueError:
                        new_coords = np.append(new_coords,np.nan)
            new_coords = np.transpose(np.column_stack(new_coords))
            coords[1] = np.hstack((new_coords,coords[1]))

def remove_nans(coords):
    """Deletes any samples with NANs in their coordinates"""
    s = np.apply_along_axis(sum,1,np.isnan(coords[1])) == 0
    coords[0] = (np.asarray(coords[0])[s]).tolist()
    coords[1] = coords[1][s,:]

def scale_custom_coords(custom_axes,coords):
    """Scales custom coordinates to match min/max of PC1"""

    to_mn = min(coords[1][:,len(custom_axes)])
    to_mx = max(coords[1][:,len(custom_axes)])

    for i in xrange(len(custom_axes)):
        from_mn = min(coords[1][:,i])
        from_mx = max(coords[1][:,i])
        coords[1][:,i] = (coords[1][:,i]  - from_mn) / (from_mx - from_mn)
        coords[1][:,i] = (coords[1][:,i]) * (to_mx-to_mn) + to_mn

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
def get_sample_ids(maptable):
    """Extracts list of sample IDs from mapping file."""
    return [line[0] for line in maptable[1:]]

def get_coord(coord_fname):
    """Opens and returns coords data"""
    try:
        coord_f = open(coord_fname, 'U').readlines()
    except (TypeError, IOError):
        raise MissingFileError, 'Coord file required for this analysis'
    coord_header, coords, eigvals, pct_var = parse_coords(coord_f)
    return [coord_header, coords, eigvals, pct_var]

def get_multiple_coords(coord_fnames):
    """Opens and returns coords data from multiple coords files"""
    # start with empty data structures
    coord_header = []
    coords = []
    edges = []

    for i,f in enumerate(coord_fnames):
        try:
            coord_f = open(coord_fnames[i], 'U').readlines()
        except (TypeError, IOError):
            raise MissingFileError, 'Coord file required for this analysis'
        coord_header_i, coords_i, eigvals_i, pct_var_i = parse_coords(coord_f)
        sampleIDs = coord_header_i
        coord_header_i = ['%s_%d' %(h,i) for h in coord_header_i]
        if i==0:
            eigvals = eigvals_i
            pct_var = pct_var_i
            coord_header = coord_header_i
            coords = coords_i
        else:
            coord_header.extend(coord_header_i)
            coords = vstack((coords,coords_i))
    # add edges from first file to others
    for id in sampleIDs:
        for i in xrange(1,len(coord_fnames)):
            edges += [('%s_%d' %(id,0), '%s_%d' %(id,i))]
    return edges, [coord_header, coords, eigvals, pct_var]

def remove_unmapped_samples(mapping,coords,edges=None):
    """Removes any samples not present in mapping file"""
    sample_IDs = zip(*mapping[1:])[0]
    # remove unmapped ids from headers and coords
    for i in xrange(len(coords[0])-1,-1,-1):
        if not coords[0][i] in sample_IDs:
            del(coords[0][i])
            coords[1] = np.delete(coords[1],i,0)
    # remove unmapped ids from edges
    if edges:
        for i in xrange(len(edges)-1,-1,-1):
            edge = edges[i]
            if not edge[0] in sample_IDs or not edge[1] in sample_IDs:
                del(edges[i])

def generate_3d_plots(prefs, data, custom_axes, background_color,label_color, \
                        dir_path='',data_file_path='',filename=None, \
                        default_filename='out'):
    """Make 3d plots according to coloring options in prefs."""
    kinpath = _make_path([(dir_path+data_file_path), filename])
    kinlink = './'+data_file_path+'/' + filename +'.kin'

    htmlpath = dir_path
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
    
    edges = None
    if data.has_key('edges'):
        edges = data['edges']

    res = make_3d_plots(coord_header, coords, pct_var,mapping,prefs, \
                        background_color,label_color,custom_axes,edges=edges)

    #Write kinemage file
    f = open(outf, 'w')
    f.write('\n'.join(res))
    f.close()
    
    #Write html page with the kinemage embedded
    f2 = open(htmlpath+filename+'_3D.html', 'w')
    f2.write("<html><head></head><body><applet code='king/Kinglet.class' \
archive='./jar/king.jar' width=800 height=600> \
<param name='kinSource' value='%s'></body></html>" % (kinlink)) 
    f2.write('\n'.join(res))
    f2.close()
    
