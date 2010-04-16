#!/usr/bin/env python
#file make_3d_plots.py

__author__ = "Jesse Stombaugh, Rob Knight, and Dan Knights"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Jesse Stombaugh", "Rob Knight", "Micah Hamady", "Dan Knights"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"


from cogent.util.misc import flatten
from qiime.parse import parse_coords,group_by_field,parse_mapping_file
from qiime.colors import natsort, get_group_colors, color_groups, make_color_dict, combine_map_label_cols, process_colorby, linear_gradient, iter_color_groups, get_map, kinemage_colors
from numpy import array, shape, apply_along_axis, dot, delete, vstack, sqrt
import numpy as np
import os
from random import choice
import re
from time import strftime
from biplots import make_mage_taxa
from qiime.util import load_pcoa_files, summarize_pcoas

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

data_colors={'blue':'#0000FF','lime':'#00FF00','red':'#FF0000', \
             'aqua':'#00FFFF','fuchsia':'#FF00FF','yellow':'#FFFF00', \
             'green':'#008000','maroon':'#800000','teal':'#008080', \
             'purple':'#800080','olive':'#808000', \
             'silver':'#C0C0C0','gray':'#808080'}
             


class MissingFileError(IOError):
    pass

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


def make_3d_plots(coord_header, coords, pct_var, mapping, prefs, \
                    background_color,label_color, \
                    taxa=None, custom_axes=None, \
                    edges=None, coords_low=None, coords_high=None, \
                    ellipsoid_prefs=None):
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
    groups_and_colors=iter_color_groups(mapping,prefs)
    groups_and_colors=list(groups_and_colors)

    for i in range(len(groups_and_colors)):  
        #Write to kinemage file using the groups, colors and coords 
        labelname=groups_and_colors[i][0]
        groups=groups_and_colors[i][1]
        colors=groups_and_colors[i][2]
        data_colors=groups_and_colors[i][3]
        data_color_order=groups_and_colors[i][4]
        
        result.extend(make_mage_output(groups, colors, coord_header, coords, \
            pct_var,background_color,label_color,data_colors, \
            taxa, custom_axes,name=labelname, \
            scaled=False, edges=edges,
            coords_low=coords_low, coords_high=coords_high, \
            ellipsoid_prefs=ellipsoid_prefs))
        result.extend(make_mage_output(groups, colors, coord_header, coords, \
            pct_var,background_color,label_color,data_colors, \
            taxa, custom_axes,name=labelname, \
            scaled=True, edges=edges, \
            coords_low=coords_low, coords_high=coords_high, \
            ellipsoid_prefs=ellipsoid_prefs))

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
                     background_color,label_color,data_colors, \
                     taxa=None, custom_axes=None,name='', \
                     radius=None, alpha=.75, num_coords=10,scaled=False, \
                     coord_scale=1.05, edges=None, coords_low=None, \
                     coords_high=None, ellipsoid_prefs=None):
    """Convert groups, colors, coords and percent var into mage format"""
    result = []
    
    #Scale the coords and generate header labels
    if scaled:
        scalars = pct_var
        if custom_axes:
            # create a dummy vector of ones to avoid scaling custom axes
            custom_scalars = scalars[0] * np.ones(len(custom_axes))
            scalars = np.append(custom_scalars,scalars)
        coords = scale_pc_data_matrix(coords, scalars)
        if not coords_low is None:
            coords_low = scale_pc_data_matrix(coords_low, scalars)
        if not coords_high is None:
            coords_high = scale_pc_data_matrix(coords_high, scalars)
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

    if not taxa is None:
        result.append('@master {taxa_points}')
        result.append('@master {taxa_labels}')

    for name, color in sorted(data_colors.items()):
        result.append(color.toMage())

    if background_color=='white':
        result.append('@whitebackground')
        result.append('@hsvcolor {black} 0.0 0.0 0.0')
    else:
        result.append('@hsvcolor {white} 180.0 0.0 100.0')
   
    #Write the groups, colors and coords
    coord_dict = dict(zip(coord_header, coords))
    if not coords_low is None:
        coord_low_dict = dict(zip(coord_header, coords_low))
    if not coords_high is None:
        coord_high_dict = dict(zip(coord_header, coords_high))
    for group_name in natsort(groups):
        ids = groups[group_name]
        result.append('@group {%s (n=%s)} collapsible' % (group_name, len(ids)))

        color = colors[group_name]
        coord_lines = []
        for id_ in sorted(ids):
            if id_ in coord_dict:
                coord_lines.append('{%s} %s' % \
                    (id_, ' '.join(map(str, coord_dict[id_][:num_coords]))))

        # create list of balls, one for each sample
        result.append('@balllist color=%s radius=%s alpha=%s dimension=%s \
master={points} nobutton' % (color, radius, alpha, num_coords))
        result.append('\n'.join(coord_lines))
        # make ellipsoids if low and high coord bounds were received
        if (not coords_low is None) and (not coords_high is None):
            # create one trianglelist for each sample to define ellipsoids
            result += make_mage_ellipsoids(ids, coord_dict, coord_low_dict,
                                           coord_high_dict, color, ellipsoid_prefs)

        # create list of labels 
        result.append('@labellist color=%s radius=%s alpha=%s dimension=%s \
master={labels} nobutton' % (color, radius, alpha, num_coords))
        result.append('\n'.join(coord_lines))

    if not taxa is None:
        result += make_mage_taxa(taxa, num_coords, pct_var,
                                 scaled=scaled, scalars=None, radius=radius)

    #Write the axes on the bottom of the graph
    result.append('@group {axes} collapsible')
    state = 'on'
    axis_mins = mins*coord_scale
    axis_maxes = maxes*coord_scale

    if not custom_axes:
        custom_axes = []
    # draw each axis
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
            
        # custom axes come first, no "percent variance" shown
        if i < len(custom_axes):
            result.append('@labellist {%s} dimension=%s %s' % \
                              (axis_names[i], num_coords, state)) 
            result.append( ('{%s}' % (axis_names[i]))  + \
                               ' '.join(map(str, end)) + ' ' + label_color)
        # if all custom axes have been drawn, draw normal PC axes
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

def make_edges_output(coord_dict, edges, num_coords,label_color,tip_fraction=0.4):
    """Creates make output to display edges (as a kinemage vectorlist).

       Params:
        coord_dict, a dict of (sampleID, coords), where coords is a numpy array
        edges, a list of pairs of sampleIDs (from, to)
        num_coords, the number of included dimensions in the PCoA plot
        label_color, the plain edge color.
        tip_fraction, the portion of each edge to be colored as the 'tip'

       Returns:
        result, a list of strings containing a kinemage vectorlist
    """
    result = []
    result.append('@vectorlist {edges} dimension=%s on' % \
                      (num_coords))
    for edge in edges:
        id_fr, id_to = edge
        # extract the coords of each vertex
        pt_fr = coord_dict[id_fr][:num_coords]
        pt_to = coord_dict[id_to][:num_coords]
        # get 'index' of the destination coords file from 'to' sampleID
        which_set = int(id_to[id_to.rindex('_')+1:]) - 1
        # different tip color for each destination coords file
        tip_color = kinemage_colors[which_set % len(kinemage_colors)]
        # plot a color 'tip' on the line (certain % of line length)
        # this gets the coords of the beginning of the 'tip'
        diffs = (pt_to-pt_fr) * (1-tip_fraction)
        middles = pt_fr + diffs
        # add a default-color line segment 
        result.append('%s %s' % \
                          (' '.join(map(str, pt_fr)),label_color))
        result.append('%s %s P' % \
                          (' '.join(map(str, middles)),label_color))
        # add the tip-colored line segment
        result.append('%s %s' % \
                          (' '.join(map(str, middles)), tip_color))
        result.append('%s %s P' % \
                          (' '.join(map(str, pt_to)), tip_color))            
    return result

def make_mage_ellipsoids(ids, coord_dict, coord_low_dict,
                         coord_high_dict, color, ellipsoid_prefs=\
                             {"smoothness":2,"alpha":.25}):
    """Makes ellipsoids with centers in coord_dict. 
       coord_low_dict and coord_high_dict are used to scale
       each axis of the ellipsoid.
    """
    alpha = ellipsoid_prefs['alpha']
    nsubdivs = ellipsoid_prefs['smoothness']
    result = []
    coord_lines = []
    for id_ in sorted(ids):
        if id_ in coord_dict:
            center = coord_dict[id_][:3]
            dims = coord_high_dict[id_][:3] - coord_low_dict[id_][:3]

            faces = make_ellipsoid_faces(center, dims, nsubdivs=nsubdivs)
            for face in faces:
                result.append("@trianglelist color=%s alpha=%f master={points} nobutton" %(color, alpha))
                for point in face:
                    result.append(' '.join(map(str,point)))
    return result

def make_ellipsoid_faces(center, dims, nsubdivs=2):
    """Returns a list of 3-tuples (triangles) of 3-tuples (points)
       defining an ellipsoid centered at center with axis 
       dimensions given in dims.

       nsubdivs determines the number of recursive divisions of
       the faces that will be made.
       nsubdivs value     No. faces
       0                  20
       1                 
    
    """
    t = (1+sqrt(5))/2.0
    s = sqrt(1+t**2)
    
    vertices = [(t/s,1/s,0), (-t/s,1/s,0), (t/s,-1/s,0),\
                (-t/s,-1/s,0), (1/s,0,t/s), (1/s,0,-t/s), (-1/s,0,t/s),(-1/s,0,-t/s),\
                    (0,t/s,1/s), (0,-t/s,1/s), (0,t/s,-1/s), (0,-t/s,-1/s)]

    v = vertices
    faces = [(v[0],v[8],v[4]),(v[1],v[10],v[7]),(v[2],v[9],v[11]),(v[7],v[3],v[1]),(v[0],v[5],v[10]),(v[3],v[9],v[6]),\
                 (v[3],v[11],v[9]),(v[8],v[6],v[4]),(v[2],v[4],v[9]),(v[3],v[7],v[11]),(v[4],v[2],v[0]),\
                 (v[9],v[4],v[6]),(v[2],v[11],v[5]),(v[0],v[10],v[8]),(v[5],v[0],v[2]),(v[10],v[5],v[7]),(v[1],v[6],v[8]),\
                 (v[1],v[8],v[10]),(v[6],v[1],v[3]),(v[11],v[7],v[5])]
    
    #subdivide each of the faces into 9 faces
    for i in xrange(nsubdivs):
        new_faces = []
        for face in faces:
            new_faces.extend(subdivide(face[0], face[1], face[2]))
        faces = new_faces
    faces = scale_faces(dims[0], dims[1], dims[2], faces)
    faces = translate_faces(center, faces)
    return faces

def subdivide(x,y,z):
    #look at x-y edge
    xy = [(x[0]*2/3.0+y[0]/3.0, x[1]*2/3.0+y[1]/3.0, x[2]*2/3.0+y[2]/3.0), (x[0]/3.0+y[0]*2/3.0, x[1]/3.0+y[1]*2/3.0, x[2]/3.0+y[2]*2/3.0)]
    #pull them to the surface of the sphere
    xy = [(i[0]/sqrt(i[0]**2+i[1]**2+i[2]**2), i[1]/sqrt(i[0]**2+i[1]**2+i[2]**2), i[2]/sqrt(i[0]**2+i[1]**2+i[2]**2)) for i in xy]
    
    #do the same for the other edges
    xz = [(x[0]*2/3.0+z[0]/3.0, x[1]*2/3.0+z[1]/3.0, x[2]*2/3.0+z[2]/3.0), (x[0]/3.0+z[0]*2/3.0, x[1]/3.0+z[1]*2/3.0, x[2]/3.0+z[2]*2/3.0)]
    xz = [(i[0]/sqrt(i[0]**2+i[1]**2+i[2]**2), i[1]/sqrt(i[0]**2+i[1]**2+i[2]**2), i[2]/sqrt(i[0]**2+i[1]**2+i[2]**2)) for i in xz]
    
    zy = [(z[0]*2/3.0+y[0]/3.0, z[1]*2/3.0+y[1]/3.0, z[2]*2/3.0+y[2]/3.0), (z[0]/3.0+y[0]*2/3.0, z[1]/3.0+y[1]*2/3.0, z[2]/3.0+y[2]*2/3.0)]
    zy = [(i[0]/sqrt(i[0]**2+i[1]**2+i[2]**2), i[1]/sqrt(i[0]**2+i[1]**2+i[2]**2), i[2]/sqrt(i[0]**2+i[1]**2+i[2]**2)) for i in zy]
    
    center = ((x[0]+y[0]+z[0])/3.0, (x[1]+y[1]+z[1])/3.0, (x[2]+y[2]+z[2])/3.0)
    center_len = sqrt(center[0]**2+center[1]**2+center[2]**2)
    center = (center[0]/center_len, center[1]/center_len, center[2]/center_len)

    #generate the new list of faces
    faces = [(x,xz[0],xy[0]), (xz[0],xy[0],center), (xz[0],xz[1],center), (xz[1],zy[0],center), (xz[1],z,zy[0]), (xy[0],center,xy[1]), 
             (xy[1],y,zy[1]), (xy[1],zy[1],center), (zy[1],zy[0],center)]
    
    return faces

def scale_faces(a,b,c,faces):
    for i,face in enumerate(faces):
        faces[i] = list(faces[i])
        for j,point in enumerate(face):
            faces[i][j] = (point[0]*a, point[1]*b, point[2]*c)
    return faces

def translate_faces(center, faces):
    for i,face in enumerate(faces):
        faces[i] = list(faces[i])
        for j,point in enumerate(face):
            faces[i][j] = (point[0]+center[0], point[1]+center[1], point[2]+center[2])
    return faces



def process_custom_axes(axis_names):
    """Parses the custom_axes option from the command line"""
    return axis_names.strip().strip("'").strip('"').split(',')

def process_coord_filenames(coord_filenames):
    """Parses the custom_axes option from the command line"""
    return coord_filenames.strip().strip("'").strip('"').split(',')

def get_custom_coords(axis_names,mapping, coords):
    """Gets custom axis coords from the mapping file.
       Appends custom as first column(s) of PCoA coords matrix.

       Params:
        axis_names, the names of headers of mapping file columns
        mapping, the mapping file object (with list of headers in element 0)
        coords, the PCoA coords object, with coords matrix in element 1
    """
    for i, axis in enumerate(reversed(axis_names)):
        if not axis in mapping[0]:
            print 'Warning: could not find custom axis',axis,'in map headers:',mapping[0]
        else:
            # get index of column in mapping file
            col_idx = mapping[0].index(axis)
            # extract column data
            col = zip(*mapping[1:])[col_idx]
            sample_IDs = zip(*mapping[1:])[0]
            new_coords = array([])
            # load custom coord for this axis for each sample ID 
            for id in coords[0]:
                if id in sample_IDs:
                    row_idx = list(sample_IDs).index(id)
                    try:
                        as_float = float(col[row_idx])
                        new_coords = np.append(new_coords,as_float)
                    except ValueError:
                        new_coords = np.append(new_coords,np.nan)
            new_coords = np.transpose(np.column_stack(new_coords))
            # append new coords to beginning column of coords matrix
            coords[1] = np.hstack((new_coords,coords[1]))

def remove_nans(coords):
    """Deletes any samples with NANs in their coordinates"""
    s = np.apply_along_axis(sum,1,np.isnan(coords[1])) == 0
    coords[0] = (np.asarray(coords[0])[s]).tolist()
    coords[1] = coords[1][s,:]

def scale_custom_coords(custom_axes,coords):
    """Scales custom coordinates to match min/max of PC1"""

    # the target min and max
    to_mn = min(coords[1][:,len(custom_axes)])
    to_mx = 2*max(coords[1][:,len(custom_axes)])

    # affine transformation for each custom axis
    for i in xrange(len(custom_axes)):
        from_mn = min(coords[1][:,i])
        from_mx = max(coords[1][:,i])
        coords[1][:,i] = (coords[1][:,i]  - from_mn) / (from_mx - from_mn)
        coords[1][:,i] = (coords[1][:,i]) * (to_mx-to_mn) + to_mn

#The following functions were not unit_tested, however the parts within
#the functions are unit_tested
def get_sample_ids(maptable):
    """Extracts list of sample IDs from mapping file."""
    return [line[0] for line in maptable[1:]]

def get_coord(coord_fname, method="IQR"):
    """Opens and returns coords location matrix and metadata.
       Also two spread matrices (+/-) if passed a dir of coord files.
       If only a single coord file, spread matrices are returned as None.
    """
    if not os.path.isdir(coord_fname):
        try:
            coord_f = open(coord_fname, 'U').readlines()
        except (TypeError, IOError):
            raise MissingFileError, 'Coord file required for this analysis'
        coord_header, coords, eigvals, pct_var = parse_coords(coord_f)
        return [coord_header, coords, eigvals, pct_var, None, None]
    else:
        master_pcoa, support_pcoas = load_pcoa_files(coord_fname)

        # get Summary statistics
        coords, coords_low, coords_high, eigval_average, coord_header = \
            summarize_pcoas(master_pcoa,support_pcoas, method=method)
        pct_var = master_pcoa[3] # should be getting this from an average

        # make_3d_plots expects coord_header to be a python list
        coord_header = list(master_pcoa[0])
        return [coord_header, coords, eigval_average, pct_var, coords_low, coords_high]

def get_multiple_coords(coord_fnames):
    """Opens and returns coords data and edges from multiple coords files.

       Params:
        coord_fnames, the names of the coordinate files

       Returns:
        edges, a list of pairs of sample IDs, (from, to)
        coords
            a list of [coord_header, coords, eigvals, pct_var]
            all coords are put in a single data matrix.
            Sample IDs from ith file have _i appended to them.
            eigvals, pct_var are taken from first coords file
    """
    # start with empty data structures
    coord_header = []
    coords = []
    edges = []

    # load all coords files into same data matrix
    for i,f in enumerate(coord_fnames):
        try:
            coord_f = open(coord_fnames[i], 'U').readlines()
        except (TypeError, IOError):
            raise MissingFileError, 'Coord file required for this analysis'
        coord_header_i, coords_i, eigvals_i, pct_var_i = parse_coords(coord_f)
        sampleIDs = coord_header_i
        # append _i to this file's sampleIDs
        coord_header_i = ['%s_%d' %(h,i) for h in coord_header_i]
        # get eigvals, pct_var from first coords file
        if i==0:
            eigvals = eigvals_i
            pct_var = pct_var_i
            coord_header = coord_header_i
            coords = coords_i
        # for second, third, etc coords files, just append to first file
        else:
            coord_header.extend(coord_header_i)
            coords = vstack((coords,coords_i))
    # add all edges
    for _id in sampleIDs:
        for i in xrange(1,len(coord_fnames)):
            # edges go from first file's points to other files' points
            edges += [('%s_%d' %(_id,0), '%s_%d' %(_id,i))]
    return edges, [coord_header, coords, eigvals, pct_var]

def get_taxa(taxa_fname, sample_ids):
    """Opens and returns coords data"""
    try:
        lines = open(taxa_fname, 'U').readlines()
    except (TypeError, IOError):
        raise MissingFileError, 'Taxa summary file required for this analysis'
    map = parse_mapping_file(lines)
    return map

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
                        default_filename='out', ellipsoid_prefs=None):
    """Make 3d plots according to coloring options in prefs."""

    if filename is None:
        filename = default_filename
    kinpath = os.path.join(data_file_path,filename) + ".kin"
    data_folder = os.path.split(data_file_path)[-1]
    kinlink = os.path.join('./',data_folder,filename) + ".kin"
    htmlpath = dir_path

    coord_header, coords, eigvals, pct_var, coords_low, coords_high = \
        data['coord']
    mapping=data['map']

    edges = None
    if data.has_key('edges'):
        edges = data['edges']

    taxa = None
    if data.has_key('taxa'):
        taxa = data['taxa']

    res = make_3d_plots(coord_header, coords, pct_var,mapping,prefs, \
                        background_color,label_color, \
                        taxa, custom_axes=custom_axes,edges=edges, \
                        coords_low=coords_low, coords_high=coords_high, \
                        ellipsoid_prefs=ellipsoid_prefs)

    #Write kinemage file
    f = open(kinpath, 'w')
    f.write('\n'.join(res))
    f.close()
    
    #Write html page with the kinemage embedded
    f2 = open(os.path.join(htmlpath,filename)+'_3D.html', 'w')
    f2.write("<html><head></head><body><applet code='king/Kinglet.class' \
archive='./jar/king.jar' width=800 height=600> \
<param name='kinSource' value='%s'></body></html>" % (kinlink)) 
    f2.write('\n'.join(res))
    f2.close()
    
