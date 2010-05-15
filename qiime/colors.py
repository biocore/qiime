#!/usr/bin/env python
#file colors.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Rob Knight","Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"

"""Code for coloring series based on prefs file.
"""
from colorsys import rgb_to_hsv, hsv_to_rgb
from parse import parse_mapping_file, group_by_field, parse_otu_table
from numpy import array
import os
import re

def string_to_rgb(s):
    """Converts hex string to RGB"""
    orig_s = s
    s = s.strip()
    if s.startswith('#'):
        s = s[1:]
    if not len(s) == 6:
        raise ValueError, "String %s doesn't look like a hex string" % orig_s
    return int(s[:2], 16), int(s[2:4],16), int(s[4:], 16)

def rgb_tuple_to_hsv(rgb):
    """Converts rgb tuple to hsv on Mage's scale"""
    rgb_0_to_1 = array(rgb)/255.0
    hsv = rgb_to_hsv(*tuple(rgb_0_to_1))
    return hsv[0]*360, hsv[1]*100, hsv[2]*100

def mage_hsv_tuple_to_rgb(hsv):
    """Converts hsv tuple on Mage scale to rgb on 0-255 scale"""
    hsv_0_to_1 = hsv[0]/360.0, hsv[1]/100.0, hsv[2]/100.0
    rgb = hsv_to_rgb(*tuple(hsv_0_to_1))
    return int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)

class Color(object):
    """Stores a color object: name, HSV, ability to write as HTML or Mage.
    
    Note: the reason we store as HSV, not RGB, is that you frequently want
    to do gradient colors by hue going from e.g. white to blue, white to red,
    etc. Unfortunately, in RGB, you can't specify _which_ white you have
    in e.g. #FFFFFF, whereas to get the right gradient you need to be able
    to specify that you want (0,0,100) or (180,0,100) or whatever. Hence
    the colorspace gymnastics.
    """

    def __init__(self, name, coords, colorspace='rgb'):
        """Returns new Color object. Init with name and coords as (R,G,B).
        
        Can also initialize with coords as (H,S,V) or #aabbcc format.
        """
        self.Name = name

        if isinstance(coords, str): #assume is hex format
            self.Coords = rgb_tuple_to_hsv(string_to_rgb(coords))
        elif colorspace == 'rgb':
            self.Coords = rgb_tuple_to_hsv(tuple(coords))
        elif colorspace == 'hsv':
            self.Coords = tuple(coords)
        else:
            raise ValueError, \
                "Unknown colorspace %s: valid values are rgb, hsv" % colorspace

    def toRGB(self):
        """Returns self as r, g, b tuple."""
        return mage_hsv_tuple_to_rgb(self.Coords)

    def toMage(self):
        """Returns self as Mage/KiNG-format string"""
        h, s, v = self.Coords
        return '@hsvcolor {%s} %3.1f %3.1f %3.1f' % (self.Name, h,s,v)

    def toHex(self):
        """Returns self as hex string."""
        rgb = self.toRGB()
        return ('#%02s%02s%02s' % (hex(rgb[0])[2:], hex(rgb[1])[2:], \
            hex(rgb[2])[2:])).replace(' ','0')

    def __str__(self):
        """Return string representation of self"""
        return str(self.Name) + ':' + self.toHex()

def _natsort_key(item):
    """Provides normalized version of item for sorting with digits.

    From: 
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    chunks = re.split('(\d+(?:\.\d+)?)', item)
    for ii in range(len(chunks)):
        if chunks[ii] and chunks[ii][0] in '0123456789':
            if '.' in chunks[ii]: numtype = float
            else: numtype = int
            # wrap in tuple with '0' to explicitly specify numbers come first
            chunks[ii] = (0, numtype(chunks[ii]))
        else:
            chunks[ii] = (1, chunks[ii])
    return (chunks, item)

def natsort(seq):
    """Sort a sequence of text strings in a reasonable order.

    From: 
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    alist = list(seq)
    alist.sort(key=_natsort_key)
    return alist

data_color_order = ['blue','lime','red','aqua','fuchsia','yellow','green', \
              'maroon','teal','purple','olive','silver','gray']

def color_dict_to_objects(d, colorspace='hsv'):
    """Converts color dict to dict of Color objects"""
    result = {}
    for k, v in d.items():
        result[k] = Color(k, v, colorspace)
    return result


#Note: these are all in Mage HSV colorspace
data_color_hsv = {
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
        'yellow':   (60,100,100)
}

data_colors = color_dict_to_objects(data_color_hsv)

kinemage_colors = ['hotpink','blue', 'lime','gold','red','sea','purple','green']

def iter_color_groups(mapping, prefs):
    """Iterates over color groups for each category given mapping file/prefs.

    See get_group_colors for details of algorithm.
    """
    #Iterate through prefs and color by given mapping labels
    for key in natsort(prefs.keys()):
        col_name = prefs[key]['column']
        if 'colors' in prefs[key]:
            if isinstance(prefs[key]['colors'], dict):
                colors = prefs[key]['colors'].copy()    #copy so we can mutate
            else:
                colors = prefs[key]['colors'][:]
        else:
            colors={}
        labelname=prefs[key]['column']
        
        #Define groups and associate appropriate colors to each group
        groups = group_by_field(mapping, col_name)
        colors, data_colors, data_color_order = \
            get_group_colors(groups, colors)

        yield labelname, groups, colors, data_colors, data_color_order

def get_group_colors(groups, colors, data_colors=data_colors, data_color_order=data_color_order):
    """Figures out group colors for a specific series based on prefs.

    Algorithm is as follows:

    - For each name, color pair we know about:
        - Check if the name is one of the groups (exact match)
        - If it isn't, assume it's a prefix and pull out all the matching groups
        - If the color is just a string, set everything to the color with that
          name
        - Otherwise, assume that either it's a new color we're adding, or that
          it's a range for gradient coloring.
        - If it's a new color, create it and add it to added_data_colors.
        - If it's a gradient, make up all the new colors and add them to
          added_data_colors

    The current method for gradient coloring of columns (should perhaps
    replace with more general method) is to pass in any of the following:

    'colors':(('white', (0,0,100)),('red',(0,100,100)))

    makes gradient between white and red, applies to all samples
    
    'colors':{'RK':(('white',(0,0,100)),('red',(0,100,100))),
              'NF':(('white',(120,0,100)),('green',(120,100,100)))
             }
    pulls the combination samples starting with RK, colors with
    first gradient, then pulls the combination samples starting
    with NF, colors with the next gradient.

    Return values are:
    - colors: dict of {group_value:color_name}
    - data_colors: dict of {color_name:color_object}
    - data_color_order: order in which the data colors are used/written.
    """
    
    added_data_colors = {}
    if isinstance(colors, dict):
        #assume we're getting some of the colors out of a dict
        if colors.items() <> []:
            for k, v in sorted(colors.items()):
                if k not in groups: #assume is prefix
                    k_matches = [g for g in groups if g.startswith(k)]
                    if isinstance(v, str):  #just set everything to this color
                        for m in k_matches:
                            colors[m] = v
                    else:   #assume is new color or range
                        first, second = v
                        if isinstance(first, str): #new named color?
                            if first not in data_colors:
                                added_data_colors[first] = Color(first, second)
                            for m in k_matches:
                                colors[m] = first
                        else:   #new color range?
                            start_color, end_color = map(get_color,
                                                            [first,second])
                            num_colors = len(k_matches)
                            curr_data_colors = color_dict_to_objects(
                                make_color_dict(start_color,
                                start_hsv,end_color,end_hsv,num_colors))
                            curr_colors = {}
                            color_groups(k_matches, curr_colors,
                                natsort(curr_data_colors))
                            colors.update(curr_colors)
                            added_data_colors.update(curr_data_colors)
                    del colors[k]
                elif not isinstance(v, str):    #assume val is new color
                    color = get_color(v)
                    if color.Name not in data_colors:
                        added_data_colors[color.Name] = color
                    colors[k] = color.Name
            #handle any leftover groups
            color_groups(groups, colors, data_color_order)
            #add new colors
            data_colors.update(added_data_colors)
            data_color_order.append(natsort(added_data_colors))
        else:
            #handle case where no prefs is used
            color_groups(groups, colors, data_color_order)
    else:
        #handle the case where colors is a tuple for gradients
        start_color, end_color = map(get_color, colors)
        start_hsv=start_color.Coords
        end_hsv=end_color.Coords
        num_colors = len(groups)
        data_colors = color_dict_to_objects(
            make_color_dict(start_color, start_hsv, end_color, 
            end_hsv, num_colors))
        data_color_order = list(natsort(data_colors.keys()))
        colors = {}
        color_groups(groups, colors, data_color_order)

    return colors, data_colors, data_color_order

def get_color(color, data_colors=data_colors):
    """Gets a color by looking up its name or initializing with name+data"""
    if isinstance(color, str):
        if color in data_colors:
            
            return data_colors[color]
        else:
            raise ValueError, "Color name %s in prefs not recognized" % color
    else:
        name, coords = color
        if isinstance(coords, str):
            colorspace = 'rgb'
        else:
            colorspace = 'hsv'
        return Color(name, coords, colorspace)

def color_groups(groups, colors, data_color_order):
    """Colors a set of groups in data_color_order, handling special colors.
    
    Modifies colors in-place.

    Cycles through data colors (i.e. wraps around when last color is reached).
    """
    group_num=-1
    for g in natsort(groups):
        if g not in colors:
            group_num+=1
            if group_num==len(data_color_order):
                group_num=0
            colors[g] = data_color_order[group_num]
 
def make_color_dict(start_name, start_hsv, end_name, end_hsv,n):
    """Makes dict of color gradient"""
    colors = linear_gradient(start_hsv,end_hsv,n)
    names = ['%sto%s%s_%s' % (start_name,end_name,n,i) for i in range(n)]
    return dict(zip(names,colors))

def combine_map_label_cols(combinecolorby,mapping):
    """Merge two or more mapping columns into one column"""
    combinedmapdata=array(['']*len(mapping),dtype='a100')
    title=[]
    match=False
    for p in range(len(combinecolorby)):                    
        for i in range(len(mapping[0])):
            if str(combinecolorby[p])==str(mapping[0][i]):
                match=True
                for q in range(len(mapping)):
                    combinedmapdata[q]=combinedmapdata[q]+mapping[q][i]
                break
            else:
                match=False
        if not match:
            raise ValueError, 'One of the columns you tried to combine does not exist!'
        title.append(combinecolorby[p])
    combinedmapdata[0]='&&'.join(title)
    for i in range(len(combinedmapdata)):
        mapping[i].append(combinedmapdata[i])
    
    return mapping

def process_colorby(colorby,data,color_prefs=None):
    """Parses the colorby option from the command line.

    color_prefs is required if colorby is not passed.
    """
    match=False
    prefs = {}
    mapping=data['map']
    colorbydata=[]
    if colorby==None and color_prefs==None:
        #if coloby option are prefs file not given, color by all categories
        #in mapping file
        colorbydata = mapping[0]
    elif colorby and color_prefs:
        #if both the colorby option and prefs file are given, use the categories
        #from the colorby option with their appropriate colors in the prefs file
        prefs_colorby = [color_prefs[i]['column'] for i in color_prefs]
        cmd_colorby=colorby.strip().strip("'").split(',')
        for i in range(len(cmd_colorby)):
            for j in range(len(prefs_colorby)):
                if cmd_colorby[i]==prefs_colorby[j]:
                    colorbydata.append(prefs_colorby[j])
                    match=True
                    break
                else:
                    match=False
            if not match:
                colorbydata.append(cmd_colorby[i])
        names = list(colorbydata)
    elif colorby:
        #if only the colorby option is passed
        colorbydata = colorby.strip().strip("'").split(',')
    else:
        #if only the prefs file is passed
        colorbydata = [color_prefs[i]['column'] for i in color_prefs]
        names = list(color_prefs)
        
    match=False
    for j, col in enumerate(colorbydata):
        key = str(col)
        #transfer over old color data if it was present
        if '&&' in col:
            #Create an array using multiple columns from mapping file
            combinecolorby=col.split('&&')
            data['map']=combine_map_label_cols(combinecolorby,mapping)
            prefs[key]={}
            prefs[key]['column']='&&'.join(combinecolorby)
        else:
            #Color by only one column in mapping file      
            prefs[key]={}
            prefs[key]['column']=col
            
        if color_prefs:
            for p in color_prefs:
                if 'column' in color_prefs[p] and color_prefs[p]['column']==col:
                    if 'colors' in color_prefs[p]:
                        prefs[key]['colors'] = color_prefs[p]['colors']
                    else:
                        prefs[key]['colors'] = (('white', (0,0,100)),('red',(0,100,100)))
                    match=True
                    break
                else:
                    match=False
            if not match:
                prefs[key]={}
                prefs[key]['column']=col
                prefs[key]['colors'] = (('white', (0,0,100)),('red',(0,100,100)))
                
    return prefs,data


def linear_gradient(start, end, nbins, eps=1e-10):
    """Makes linear color gradient from start to end, using nbins.
    
    Returns list of (x, y, z) tuples in current colorspace.
    eps is used to prevent the case where start and end are the same.
    """
    start = array(start)
    end = array(end)
    result = []
    n_minus_1 = max(float(nbins-1), eps)
    for i in range(nbins):
        result.append(list((start*(n_minus_1-i)/n_minus_1)+(end*(i/n_minus_1))))
    return result

#The following functions were not unit_tested, however the parts within
#the functions are unit_tested
def get_map(options, data):
    """Opens and returns mapping data"""
    try:
        map_f = open(options.map_fname, 'U').readlines()
    except (TypeError, IOError):
        raise MissingFileError, 'Mapping file required for this analysis'
    data['map'] = parse_mapping_file(map_f)
    return data['map']

def map_from_coords(coords):
    """Makes pseudo mapping file from coords.
    
    set data['map'] to result of this if coords file supplied but not map.

    TODO: write equivalent function for other inputs, e.g. for rarefaction --
    basic principle is that you need data structure that you can extract list
    of sample ids from.
    """
    result=(([['SampleID','Sample']]))
    for i in range(len(data['coord'][0])):
            data['map'].append([data['coord'][0][i],'Sample'])

def sample_color_prefs_and_map_data_from_options(options):
    """Returns color prefs and mapping data based on options.
    
    Note: opens files as needed. Only returns the info related to metadata
    coloring and category maps. If you need additional info, it is necessary
    to get that info explicitly (e.g. coord files, rarefaction files, etc.).

    For example, you might modify the data dict afterwards to add coords,
    rarefaction info, etc. depending on the application.
    """
    data = {}

    #Open and get mapping data, if none supplied create a pseudo mapping \
    #file

    mapping,headers,comments = get_map(options, data)
    new_mapping=[]
    new_mapping.append(headers)
    for i in range(len(mapping)):
        new_mapping.append(mapping[i])
    data['map']=new_mapping
    #need to set some other way from sample ids
    #Determine which mapping headers to color by, if none given, color by \
    #Sample ID's
    
    try:
        colorby = options.colorby
    except AttributeError:
        colorby = None
    
    if options.prefs_path:
        prefs = eval(open(options.prefs_path, 'U').read())
        color_prefs, data=process_colorby(colorby, data, \
                                                prefs['sample_coloring'])
        
        if prefs.has_key('background_color'):
            background_color= prefs['background_color']
        else:
            background_color='black'
    else:
        background_color='black'
        color_prefs, data=process_colorby(colorby, data, None)
    
    if options.prefs_path and options.background_color:
        background_color=options.background_color
    elif options.background_color:
        background_color=options.background_color

    if background_color=='black':
        label_color='white'
    else:
        label_color='black'

    return color_prefs,data,background_color,label_color

def taxonomy_color_prefs_and_map_data_from_options(options):
    """Returns color prefs and counts data based on options.
       counts data is any file in a format that can be parsed by parse_otu_table  
    """
    data = {}
    data['counts'] = {}
    taxonomy_levels = []

    #need to set some other way from sample ids
    #Determine which mapping headers to color by, if none given, color by \
    #Sample ID's
    taxonomy_count_files = options.counts_fname.strip().strip("'").split(',')
    for f in taxonomy_count_files:
        try:
            counts_f = open(f, 'U').readlines()
        except (TypeError, IOError):
            raise MissingFileError, 'Counts file required for this analysis'
        sample_ids, otu_ids, otu_table, lineages = \
                       parse_otu_table(counts_f,count_map_f=float)

        data['counts'][f] = (sample_ids, otu_ids, otu_table, lineages)
        level = max([len(t.split(';')) - 1 for t in otu_ids])

        taxonomy_levels.append(str(level))

    if options.prefs_path:
        prefs = eval(open(options.prefs_path, 'U').read())
        color_prefs=taxonomy_process_prefs(taxonomy_levels, \
                                                prefs['taxonomy_coloring'])

        if prefs.has_key('background_color'):
            background_color= prefs['background_color']
        else:
            background_color='black'
    else:
        background_color='black'
        color_prefs=taxonomy_process_prefs(taxonomy_levels, None)

    if options.prefs_path and options.background_color:
        background_color=options.background_color
    elif options.background_color:
        background_color=options.background_color

    if background_color=='black':
        label_color='white'
    else:
        label_color='black'
    return color_prefs,data,background_color,label_color


def taxonomy_process_prefs(taxonomy_levels, color_prefs = None):
    """Creates taxonomy prefs dict given specific taxonomy levels.

            color_prefs is not required
            taxonomy_levels is a list of the level number i.e. Phylum is 2
            prefs will include a 'colors' dictionary for each given level
               if there is a cooresponding level in color_prefs that is the
               dictionary for the level otherwise it adds and empty dict
    """
    prefs = {}
    for j, col in enumerate(taxonomy_levels):
        key = str(col)
        col = str(col)
        #Color by only one level    
        prefs[key]={}
        prefs[key]['column']=col

        if color_prefs:
            for p in color_prefs:
                if 'column' in color_prefs[p] and str(color_prefs[p]['column'])==col:
                    if 'colors' in color_prefs[p]:
                        prefs[key]['colors'] = color_prefs[p]['colors'].copy()
                    else:
                        prefs[key]['colors'] = {}
                    match=True
                    break
                else:
                    match=False
            if not match:
                prefs[key]={}
                prefs[key]['column']=col
                prefs[key]['colors'] = {}

    return prefs
