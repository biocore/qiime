#!/usr/bin/env python
from __future__ import division

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jeremy Widmann","Rob Knight","Jesse Stombaugh",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Jeremy Widmann"
__email__ = "jeremy.widmann@colorado.edu"

from matplotlib import use
use('Agg',warn=False)
from qiime.parse import parse_mapping_file, parse_distmat, group_by_field,\
    group_by_fields
from qiime.pycogent_backports.test import t_two_sample
from numpy import array, mean, average, arange, concatenate
from collections import defaultdict
from string import strip
from matplotlib.pylab import savefig, clf, gca, gcf,close
from cogent.draw.util import hist
from matplotlib.patches import Ellipse, Polygon
from random import choice
from numpy.random import permutation
from qiime.colors import data_colors, Color, rgb_tuple_to_hsv, \
    mage_hsv_tuple_to_rgb
from math import ceil
from os import mkdir, path

def matplotlib_rgb_color(rgb_color):
    """Returns RGB color in matplotlib format.
       ex: (255,0,255) will return (1.0,0.0,1.0)
    """
    return tuple([i/255. for i in rgb_color])

def average_colors(color1, color2):
    """Returns average of two RGB colors. 
       - color1 and color2 are RGB tuples.
    """
    avg_list = []
    for i,j in zip(color1,color2):
        avg_list.append((i+j)/2.)
    return tuple(avg_list)

def average_all_colors(paired_field_names, field_to_color_prefs):
    """Returns dict of paired_field_names matched to color.
    
        - paired_field_names: list of names of fields as follows:
            FieldName_FieldData1_to_FieldData2.
            ex:
            Hand_Right_to_Left
        - field_to_color_prefs:
            dict mapping field name to groups, colors, data_colors, 
            data_color_order from qiime.colors.iter_color_groups
    """
    paired_field_to_color = {}
    for name in paired_field_names:
        field, data = name.split('###FIELDDATA###')
        first,second = data.split('_to_')
        field_prefs = field_to_color_prefs[field]
        color1 = field_prefs[1][first]
        color2 = field_prefs[1][second]
        color1_rgb = field_prefs[2][color1].toRGB()
        color2_rgb = field_prefs[2][color2].toRGB()
        avg_color = average_colors(color1_rgb,color2_rgb)
        paired_field_to_color[name]=matplotlib_rgb_color(avg_color)
    return paired_field_to_color

def assign_unassigned_colors(unassigned_colors, color_dict=data_colors):
    """Assigns unassigned_colors in order.
    """
    label_to_color = {}
    num_colors = len(data_colors)
    color_list = data_colors.values()
    for i, label in enumerate(unassigned_colors):
        color_index = i % num_colors
        color = matplotlib_rgb_color(color_list[color_index].toRGB())
        label_to_color[label]=color
    return label_to_color

def assign_mapped_colors(mapped_labels, field_to_color_prefs):
    """Returns map of label to color using color prefs.
    """
    label_to_color = {}
    
    for name in mapped_labels:
        field, data = name.split('_Within_')
        label = data.split('_Distances')[0]
        field_prefs = field_to_color_prefs[field]
        color = field_prefs[1][label]

        color_rgb = field_prefs[2][color].toRGB()

        label_to_color[name]=matplotlib_rgb_color(color_rgb)
    return label_to_color

def between_sample_distances(dmat):
    """Returns all upper triangle distances from dmat.

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.
    """
    distances = []
    dmat_len = len(dmat)
    for i in range(dmat_len):
        for j in range(i+1,dmat_len):
            distances.append(dmat[i][j])
    return {'All_Between_Sample_Distances':distances}

def within_category_distances_grouped(single_field, \
    label_suffix='_All_Within_Category_Distances'):
    """Returns all within category distances grouped for every field.
    
        - single_field is from calling group_distances and taking single_field
            result.
    """
    distances = defaultdict(list)
    for field, groups in single_field.items():
        for data in groups:
            if data[0] == data[1]:
                all = array(data[2])
                distances[field+label_suffix].extend(all.flat)
    return distances

def between_category_distances_grouped(single_field,\
    label_suffix='_All_Between_Category_Distances'):
    """Returns all between category distances grouped for every field.
    
        - single_field is from calling group_distances and taking single_field
            result.
    """
    distances = defaultdict(list)
    for field, groups in single_field.items():
        for data in groups:
            if data[0] != data[1]:
                all = array(data[2])
                distances[field+label_suffix].extend(all.flat)
        
    return distances

def within_category_distances(single_field):
    """Returns all within category distances, broken down by category.
    
        - single_field is from calling group_distances and taking single_field
            result.
    """
    distances = defaultdict(list)
    for field, groups in single_field.items():
        for data in groups:
            if data[0] == data[1]:
                all = array(data[2])
                distances[field+'_Within_'+data[0]+\
                    '_Distances'].extend(all.flat)
    return distances

def within_and_between_fields(paired_field):
    """Returns all within field and all between field comparisons.
    
        - paired_field is from calling group_distances and taking paired_field
            result.
    """
    distances = defaultdict(list)
    for field, groups in paired_field.items():
        first, second = field.split('_to_')
        if first == second:
            for data in groups:
                if data[0] != data[1]:
                    continue
                all = array(data[2])
                distances['Within_All_Fields'].extend(all.flat)
        else:
            for data in groups:
                if data[0] == data[1]:
                    continue
                all = array(data[2])
                distances['Between_All_Fields'].extend(all.flat)
    return distances

def all_category_distances(single_field):
    """Returns all within category distances, broken down by category.
    
        - single_field is from calling group_distances and taking single_field
            result.
    """
    distances = defaultdict(list)
    for field, groups in single_field.items():
        for data in groups:
            all = array(data[2])
            distances[field + '###FIELDDATA###' + data[0] + '_to_' + \
                      data[1]].extend(all.flat)
    return distances

def draw_all_histograms(single_field, paired_field, dmat, histogram_dir,\
    field_to_color_prefs,background_color):
    """Draws all combinations of histograms.
    """
    #make dict of label to histogram filename
    label_to_histogram_filename = {}
    
    #list of different distances
    distances_dict = {}
    
    #Unassigned color list.  These must be manually assigned.
    unassigned_colors = []
    
    #all category colors list.  These must be averaged for display.
    colors_to_average = []
    
    #list of colors to map from prefs.
    colors_to_map = []
    
    #Get all between sample distances
    all_between = between_sample_distances(dmat)
    distances_dict['All_Between_Sample_Distances']=all_between
    
    #add to unassigned colors list.
    unassigned_colors.extend(all_between.keys())
    
    label_to_histogram_filename.update(\
        _make_histogram_filenames(all_between,histogram_dir))
    
    #Get all within category distances grouped together
    all_within_category_grouped = \
        within_category_distances_grouped(single_field)
        
    distances_dict['All_Within_Category_Grouped']=all_within_category_grouped
    
    #add to unassigned colors list.
    unassigned_colors.extend(all_within_category_grouped.keys())

    label_to_histogram_filename.update(\
        _make_histogram_filenames(all_within_category_grouped,histogram_dir))

    #Get all between category distances grouped together
    all_between_category_grouped = \
        between_category_distances_grouped(single_field)
    distances_dict['All_Between_Category_Grouped']=all_between_category_grouped
    
    #add to unassigned colors list.
    unassigned_colors.extend(all_between_category_grouped.keys())
    
    label_to_histogram_filename.update(\
        _make_histogram_filenames(all_between_category_grouped,histogram_dir))

    
    #Get all within category distances by category
    all_within_category_individual = \
        within_category_distances(single_field)
    distances_dict['All_Within_Categories']=all_within_category_individual
    
    colors_to_map.extend(all_within_category_individual)
    
    label_to_histogram_filename.update(\
        _make_histogram_filenames(all_within_category_individual,histogram_dir))
    
    #Get all within and between field distances
    all_within_and_between_fields = within_and_between_fields(paired_field)
    distances_dict['All_Within_And_Between_Fields']=\
        all_within_and_between_fields
    
    #add to unassigned colors list.
    unassigned_colors.extend(all_within_and_between_fields.keys())

    label_to_histogram_filename.update(\
        _make_histogram_filenames(all_within_and_between_fields,histogram_dir))

    #Get all category distances
    all_categories = all_category_distances(single_field)
    distances_dict['All_Category_Pairs']=all_categories

    #add colors to be averaged
    colors_to_average.extend(all_categories.keys())
    
    label_to_histogram_filename.update(\
        _make_histogram_filenames(all_categories,histogram_dir))
    
    #assign all colors
    label_to_color = {}
    label_to_color.update(average_all_colors(colors_to_average,\
        field_to_color_prefs))
    label_to_color.update(assign_unassigned_colors(unassigned_colors))
    label_to_color.update(assign_mapped_colors(colors_to_map,\
        field_to_color_prefs))
        
    max_val = max(dmat.flat)
    bin_size = max_val/20.
    BINS=arange(0,ceil(max_val)+.01,bin_size)
    xscale, yscale = get_histogram_scale(distances_dict,nbins=BINS)
    #draw histograms

    for d_name, d_dict in distances_dict.items():
        for field, data in d_dict.items():
            #If there are no distances, remove from distances_dict and 
            # label_to_histogram_filename
            if len(data) < 1:
                label_to_histogram_filename.pop(field)
                continue
            
            color = label_to_color[field]
            outfile_name = label_to_histogram_filename[field]
            histogram = draw_histogram(distances=data,color=color,nbins=BINS, \
                outfile_name=outfile_name,xscale=xscale,yscale=yscale,\
                background_color=background_color)
    
    return distances_dict, label_to_histogram_filename

def get_histogram_scale(distances_dict, nbins):
    """Draws histogram to outfile_name.
    """
    scale_dict = defaultdict(list)
    #draw histograms
    for d_dict in distances_dict.values():
        for i, (field, data) in enumerate(d_dict.items()):
            if len(data) < 1:
                continue
            histogram = hist(data,bins=nbins)
            
            fig  = gcf()
            axis = fig.gca()

            #get height scale: y/x
            ymin,ymax = axis.get_ylim()
        
            xmin,xmax = axis.get_xlim()
            scale_dict['ymin'].append(ymin)
            scale_dict['ymax'].append(ymax)
            scale_dict['xmin'].append(xmin)
            scale_dict['xmax'].append(xmax)

            clf()
    
    yscale = (min(scale_dict['ymin']),max(scale_dict['ymax']))
    xscale = (min(scale_dict['xmin']),max(scale_dict['xmax']))
    
    return xscale,yscale

def draw_histogram(distances, color, nbins, outfile_name,\
    xscale=None, yscale=None, background_color='white',\
    title='', **kwargs):
    """Draws histogram to outfile_name.
    """    
    average = mean(distances)
    maximum = max(distances)
    histogram = hist(distances,bins=nbins,facecolor=color, \
        normed=True,**kwargs)
    
    fig  = gcf()
    axis = fig.gca()
    
    #set labels
    axis.set_xlabel('Distance')
    axis.set_ylabel('Normalized Counts of Pairs in Group')
    axis.set_title(title)
    #get figure scale: width/height
    fig_scale = fig.get_figwidth()/float(fig.get_figheight())
    
    if xscale is not None:
        axis.set_xlim(xscale)
    if yscale is not None:
        axis.set_ylim(yscale)
    
    #get height scale: y/x
    ylim = axis.get_ylim()
    ylen = ylim[1]-ylim[0]
    xlim = axis.get_xlim()
    xlen = xlim[1]-xlim[0]
    height_scale = ylen/float(xlen)
    
    #set width
    width = xlen/20.
    height = width*height_scale*fig_scale
    
    #draw circle at average distance
    ellipse = Ellipse([average,0.0],width=width, \
        height=height, \
        edgecolor=color, fill=False)
    axis.add_artist(ellipse)        
    #draw line at center of circle
    y1 = -height/2.
    y2 = height/2.

    line = Polygon([[average, y1] ,[average, y2]], edgecolor=color)
    axis.add_artist(line)
    
    transparent=True
    if background_color != "white":
        axis.set_axis_bgcolor(background_color)
        transparent=False

    savefig(outfile_name,format='png',dpi=72, transparent=transparent)

    close()
    return histogram

NAV_HTML_TR = '''<span class="smnorm"><input type="checkbox" id="%s" %s onclick="visibilityAndOpacity(this, %s)" />
<a onmouseover="mouseoverVisible(%s)"; onmouseout="mouseoverHidden(%s)">%s</a></span><br />'''

NAV_HTML_TR_BREAK = '''<span class="normal">%s</span><br />'''

NAV_HTML_FRAME_START = '''<td>
    <div style="overflow:scroll;white-space:nowrap;width:300px;height: 400px;">
    <p>
'''
NAV_HTML_FRAME_END = '''
    </p>
    </div>
    </td>
'''

FULL_HTML_TITLE_FRAME = """
<html><head> <title>
%s
</title>
<script type="text/javascript" src="./js/histograms.js"></script>
 <style type="text/css">
.smnorm {color: blue; font-family:Arial,Verdana; font-size:10; font-weight: bold;}
.normal {color: black; font-family:Arial,Verdana; font-size:11; font-weight: bold;}
"""
FULL_HTML_JS_FRAME = """
</style>
</head>
<body> 
<div id="overDiv" style="position:absolute; visibility:hidden; z-index:1000;"></div>
"""

FULL_HTML_MAIN_IMAGE = '''<div>
    <table>
        <tr>
            %s
        </tr>
    </table>
</div>

'''

SINGLE_IMAGE_BLOCK = '''
<td style="position:absolute;left:0;">
<img style="z-index:1; opacity:1.0;filter:alpha(opacity=100); visibility:%s;" id="%s" name="%s" src="%s" border="0"></td>'''

FULL_HTML_NAV_FRAME = """
<table style="position:absolute;left:600">
    <tr>
    %s
    </tr>
</table>
"""

def make_nav_html(distances_dict, label_to_histogram_filename, \
    default='All_Between_Sample_Distances'):
    """Returns HTML string with mouseover table linked to filenames.
    """
    html_list = []
    for main_block, distances in distances_dict.items():
        html_list.append(NAV_HTML_TR_BREAK%(main_block))
        for sub_label in sorted(distances.keys()):
            if sub_label not in label_to_histogram_filename:
                continue
            hist_filename = label_to_histogram_filename[sub_label].strip('./')
            hist_filename = "'"+hist_filename+"'"
            #sub_label=sub_label.replace('###FIELDDATA###','_')
            checked = ''
            if main_block == default:
                checked = 'checked'
            sub_label_quoted = "'"+sub_label+"'"
            html_list.append(NAV_HTML_TR%\
                ('check_'+sub_label, checked, \
                    sub_label_quoted, sub_label_quoted,\
                    sub_label_quoted, sub_label.replace('###FIELDDATA###','_')))
    nav_html = '\n'.join(html_list)
    return NAV_HTML_FRAME_START + nav_html + NAV_HTML_FRAME_END

def make_main_html(distances_dict, label_to_histogram_filename, root_outdir, \
    outfile_name, title='QIIME - Distance Histograms', \
    default='All_Between_Sample_Distances'):
    """Returns full HTML string to show distance histograms.
    """
    main_html_list = []
    #Add title HTML
    main_html_list.append(FULL_HTML_TITLE_FRAME%(title))
    
    #Add javascript in html
    main_html_list.append(FULL_HTML_JS_FRAME)
    
    #Add default image HTML
    all_images = []
    all_images.append(SINGLE_IMAGE_BLOCK%('visible',default,'visible',\
        label_to_histogram_filename[default]))
        
    for i,(name,src) in enumerate(label_to_histogram_filename.items()):
        all_images.append(SINGLE_IMAGE_BLOCK%('hidden',\
            name,'hidden',src))
        
    main_html_list.append(FULL_HTML_MAIN_IMAGE%('\n'.join(all_images)))
    
    #Add nav html
    nav_html = make_nav_html(distances_dict, label_to_histogram_filename)
    main_html_list.append(FULL_HTML_NAV_FRAME%(nav_html))
    
    main_html_out = open(root_outdir+'/'+outfile_name,'w')
    main_html_out.write(''.join(main_html_list))
    main_html_out.close

def get_valid_indices(input_items, wanted_items):
    """returns indices of wanted_items in input_items if present."""
    try:
        return map(input_items.index, wanted_items)
    except ValueError:  #missing ids?
        return [input_items.index(i) for i in wanted_items\
            if i in input_items]

def distances_by_groups(distance_header, distance_matrix, groups):
    """Splits distances by group membership, returns vals for each pair.
    
    Omits the zeros along the diagonal.

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.
    """
    result = []
    group_items = groups.items()

    for i, (row_group, row_ids) in enumerate(group_items):
        row_indices = get_valid_indices(distance_header, row_ids)
        #handle the case where indices are separate: just return blocks
        for j in range(i+1, len(groups)):
            col_group, col_ids = group_items[j]
            col_indices = get_valid_indices(distance_header, col_ids)
            vals = distance_matrix[row_indices][:,col_indices]
            result.append([row_group, col_group, vals])
            
        #handle the case where indices are the same so need to omit diag
        block = distance_matrix[row_indices][:,row_indices]
        size = len(row_indices)
        indices = []
        for i in range(size):
            for j in range(i,size):
                if i != j:
                    indices.append(block[i][j])
        result.append([row_group, row_group, array(indices)])
    return result

def write_distance_files(group_distance_dict,dir_prefix = '', \
    subdir_prefix='distances'):
    """writes distance files for each col of mapping file."""
    path_prefix = path.join(dir_prefix,subdir_prefix)
    try:
        mkdir(path_prefix)
    except OSError:     #raised if dir exists
        pass
        
    for field, data in group_distance_dict.items(): #skip sample id field
        fname = path.join(path_prefix, 'dist_' + field + '.txt')
        outfile = open(fname, 'w')
        for d in data:
            if subdir_prefix.endswith('pairs'):
                outfile.write('\t'.join([':'.join(d[0]) + '_to_' + \
                              ':'.join(d[1])] + map(str, d[2].flat)))
            else:
                outfile.write('\t'.join([str(d[0]) + '_to_' + \
                              str(d[1])] + map(str, d[2].flat)))
            outfile.write('\n')
        outfile.close()


def group_distances(mapping_file,dmatrix_file,fields,dir_prefix='',\
    subdir_prefix='group_distances'):
    """Calculate all lists of distance groups.
    
    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.
    """
    distance_groups = {}
    mapping, header, comments = parse_mapping_file(open(mapping_file,'U'))
    header = [header]
    header.extend(mapping)
    mapping=header

    distance_header, distance_matrix = \
        parse_distmat(open(dmatrix_file,'U'))

    if fields == []:
        raise ValueError, 'Since no fields were defined and the values within your fields are either all the same or all unique, a field was not chosen for analysis. Please define a field to analyse.'
        
    single_field = defaultdict(dict)
    for i in range(len(fields)):
        field = fields[i]
        groups = group_by_field(mapping, field)
        data = distances_by_groups(distance_header, distance_matrix, groups)
        #Need to remove pound signs from field name.
        field_name = field.replace('#','')
        single_field[field_name]=data

    write_distance_files(group_distance_dict=single_field,\
        dir_prefix=dir_prefix,subdir_prefix=subdir_prefix+'_single')
        
    paired_field = defaultdict(dict)
    paired_field_for_writing = defaultdict(dict)
    for i in range(len(fields)):
        for j in range(i,len(fields)):
            fieldi = fields[i]
            fieldj = fields[j]
            groups = group_by_fields(mapping, [fieldi,fieldj])
            data = distances_by_groups(distance_header, distance_matrix, groups)
            paired_field[fieldi+'_to_'+fieldj]=data
            paired_field_for_writing[fieldi+'_to_'+field]=data
    
    write_distance_files(group_distance_dict=paired_field_for_writing,\
        dir_prefix=dir_prefix,subdir_prefix=subdir_prefix+'_pairs')
    
    return single_field, paired_field, distance_matrix

def build_monte_carlo_prefs(fields,default_iters):
    """Builds prefs dict for monte_carlo_group_distances when not provided.
    """
    field_to_iters = dict([(f,default_iters) for f in fields])
    prefs = {'MONTE_CARLO_GROUP_DISTANCES':field_to_iters}
    return prefs

def monte_carlo_group_distances(mapping_file, dmatrix_file, prefs, \
    dir_prefix = '', subdir_prefix='monte_carlo_group_distances',\
    default_iters=10, fields=None):
    """Calculate Monte Carlo stats for specified group distances.
    
    Specifically:
    - find the groups for each specified col (or combination of cols)
    - do t test between each pair of groups
    - randomize matrix n times and find empirical value of t for each pair
    - compare the actual value of t to the randomized values

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.
    """
    mapping, header, comments = parse_mapping_file(open(mapping_file,'U'))
    header = [header]
    header.extend(mapping)
    mapping=header

    distance_header, distance_matrix = \
        parse_distmat(open(dmatrix_file,'U'))

    orig_distance_matrix = distance_matrix.copy()

    path_prefix = path.join(dir_prefix,subdir_prefix)
    
    #if dir doesn't exist
    if not path.isdir(path_prefix):
        # make directory
        mkdir(path_prefix)
    
    if fields is None:
        fields = [mapping[0][0]]
        
    if prefs is None:
        prefs = {}
 
    if 'MONTE_CARLO_GROUP_DISTANCES' not in prefs:
        prefs = build_monte_carlo_prefs(fields,default_iters)
            
    for field, num_iters in prefs['MONTE_CARLO_GROUP_DISTANCES'].items():
        if '&&' in field:
            groups = group_by_fields(mapping, field.split('&&'))
        else:
            groups = group_by_field(mapping, field)
        outfile = open(path.join(path_prefix,
                                 'group_distances_'+field+'.txt'), 'w')
        outfile.write('\t'.join(['Category_1a','Category_1b','Avg',\
            'Category_2a','Category_2b','Avg','t','p',\
            'p_greater','p_less','Iterations\n']))
        real_dists = distances_by_groups(distance_header, distance_matrix,\
            groups)
 
        #iterate over the groups
        for i, (first_g1, second_g1, distances_g1) in \
            enumerate(real_dists[:-1]):

            real_dist_1 = average(distances_g1)

            #then for each other pair (not including same group)
            for j in range(i+1,len(real_dists)):
                first_g2, second_g2, distances_g2 = real_dists[j]

                real_dist_2 = average(distances_g2)

                # permute distances just within these groups!
                rand_dists_1, rand_dists_2 = \
                        permute_between_groups(distances_g1, 
                                               distances_g2,
                                               num_iters)

                ttests = [t_two_sample(rand_dists_1[n].flatten(),rand_dists_2[n].flatten())[0] \
                    for n in range(num_iters)]
                real_ttest = t_two_sample(distances_g1.flatten(), distances_g2.flatten())
                curr_line = [first_g1, second_g1, real_dist_1, \
                    first_g2, second_g2, real_dist_2]
                curr_line.extend([real_ttest[0], real_ttest[1],\
                    (array(ttests)>real_ttest[0]).sum()/float(num_iters), \
                    (array(ttests)<real_ttest[0]).sum()/float(num_iters), \
                    num_iters])
                outfile.write('\t'.join(map(str, curr_line)))
                outfile.write('\n')


def monte_carlo_group_distances_within_between(single_field, \
    paired_field, dmat, dir_prefix = '', \
    subdir_prefix='monte_carlo_group_distances',\
    num_iters=10):
    """Calculate Monte Carlo stats within and between fields.
    
    Specifically:
    - find the groups for each specified col (or combination of cols)
    - do t test between each pair of groups
    - randomize matrix n times and find empirical value of t for each pair
    - compare the actual value of t to the randomized values

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.
    """

    path_prefix = path.join(dir_prefix,subdir_prefix)
    #if dir doesn't exist
    if not path.isdir(path_prefix):
        # make directory
        mkdir(path_prefix)
    
    real_dists = []
    within_category_distances = \
        within_category_distances_grouped(single_field,label_suffix='')
    real_dists.extend([['Within',field,distances] for field,\
        distances in within_category_distances.items()])
        
    between_category_distances = \
        between_category_distances_grouped(single_field,label_suffix='')
    real_dists.extend([['Between',field,distances] for field,\
        distances in between_category_distances.items()])
    
    within_and_between = \
        within_and_between_fields(paired_field)
    
    real_dists.extend([[field.split('_',1)[0],\
        field.split('_',1)[1],distances] for \
        field, distances in within_and_between.items()])
    
    outfile = open(path.join(path_prefix,
                            'group_distances_within_and_between.txt'), 'w')
    outfile.write('\t'.join(['Comparison','Category_1','Avg',\
        'Comparison','Category_2','Avg','t','p',\
        'p_greater','p_less','Iterations\n']))

    rand_distances = get_random_dists(real_dists, dmat, num_iters)
    
    #iterate over the groups
    for i, (first_g1, second_g1, distances_g1) in \
        enumerate(real_dists[:-1]):
        real_dist_1 = average(distances_g1)
        rand_dists_1 = [rand_distances[n][i][-1] for n in range(num_iters)]
        #then for each other pair (not including same group)
        for j in range(i+1,len(real_dists)):
            first_g2, second_g2, distances_g2 = real_dists[j]
            real_dist_2 = average(distances_g2)
            rand_dists_2 = [rand_distances[n][j][-1] \
                for n in range(num_iters)]
            ttests = [t_two_sample(rand_dists_1[n],rand_dists_2[n])[0] \
                for n in range(num_iters)]
            real_ttest = t_two_sample(distances_g1, distances_g2)
            curr_line = [first_g1, second_g1, real_dist_1, \
                first_g2, second_g2, real_dist_2]
            curr_line.extend([real_ttest[0], real_ttest[1],\
                (array(ttests)>real_ttest[0]).sum()/float(num_iters), \
                (array(ttests)<real_ttest[0]).sum()/float(num_iters), \
                num_iters])
            outfile.write('\t'.join(map(str, curr_line)))
            outfile.write('\n')
    
def permute_between_groups(g1, g2, num_iters, permute_f=permutation):
    """Returns num_iters permuted versions of g1, g2.
       Values are permuted between g1 and g2.
    
       g1, g2 are numpy arrays.
    """
    g1 = g1.flatten()
    g2 = g2.flatten()
    n1 = len(g1)
    n2 = len(g2)
    n = n1 + n2
    combined = concatenate([g1, g2])
    # generate a list of all permutations
    perms = [permute_f(n) for i in xrange(num_iters)]
    # use random permutations to split into groups
    rand_g1 = [combined[perm[:n1]] for perm in perms]
    rand_g2 = [combined[perm[n1:n]] for perm in perms]
    return rand_g1, rand_g2

def permute_for_monte_carlo(dist_matrix):
    """Returns permuted copy of distance matrix for Monte Carlo tests."""
    size = len(dist_matrix)
    p = permutation(size)
    return dist_matrix[p][:,p]

def get_random_dists(real_dists, dmat, num_iters):
    """Returns random distances same size as real_dists.
    
        - real_dists: list of distances:
            [[first_group, second_group, distancdes],...]
        - dmat: full distance matrix
        - num_iters: integer number of random dmats to make.

    WARNING: Only symmetric, hollow distance matrices may be used as input.
    Asymmetric distance matrices, such as those obtained by the UniFrac Gain
    metric (i.e. beta_diversity.py -m unifrac_g), should not be used as input.
    """
    rand_dists = []
    upper_triangle = between_sample_distances(dmat).values()[0]
    for i in range(num_iters):
        curr_rand_dists = []
        for first,second,real_dist in real_dists:    
            curr_rand_dists.append([first,second,[choice(upper_triangle) for j \
                in range(len(real_dist))]])
        rand_dists.append(curr_rand_dists)
    return rand_dists
        

def _make_histogram_filenames(distances,histogram_dir):
    """From distances dict, returns dict of keys to histogram filenames.
    
        - distances: dict of label to distances.
        
        returns dict of label to filename: {label: label_randomchars.png}
    """
    filename_dict = {}
    
    for label in distances.keys():
        filename_dict[label]=_make_random_filename(base_dir=histogram_dir, \
            suffix='.png')
    
    return filename_dict

def _make_relative_paths(label_to_path_dict, prefix):
    """Returns relative path from full path where prefix is replaced with ./
    """
    label_to_path_dict_relative = {}
    for k,v in label_to_path_dict.items():
        label_to_path_dict_relative[k] = v.replace(prefix,'./',1)
    return label_to_path_dict_relative

def _make_random_filename(base_dir='',suffix='',num_chars=20):
    """Returns filename with random characters between prefix and suffix.
    """
    all = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    rand_region = ''.join([choice(all) for i in range(num_chars)])
    return path.join(base_dir,rand_region+suffix)


