#!/usr/bin/env python
#file make_2d_plots.py

__author__ = "Jesse Stombaugh and Micah Hamady"
__copyright__ = "Copyright 2009, the 454 Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Prototype"

"""
Author: Jesse Stombaugh (jesse.stombaugh@colorado.edu) and Micah Hamady
Status: Prototype

Requirements:
MatPlotLib 0.98.5.2
Python 2.5

Example 1: Create 2D plots from only the pca/pcoa data:
Usage: python make_2d_plots.py -i raw_pca_data.txt

Example 2: Create two separate html files, one for Day and one for Type:
Usage: python make_2d_plots.py -i raw_pca_data.txt -m input_map.txt -b 'Day,Type'

Example 3: Create 2D plots for a combination of label headers from a mapping 
file:
Usage: python make_2d_plots.py -i raw_pca_data.txt -m input_map.txt 
-b 'Type&&Day' -o ./test/

"""

import matplotlib,re

matplotlib.use('Agg')
from matplotlib.pylab import *
from commands import getoutput
from string import strip
from parse import parse_map,parse_coords,group_by_field,group_by_fields
from numpy import array
from optparse import OptionParser
from time import strftime
import shutil
import os

from random import choice, randrange
from make_3d_plots import combine_map_label_cols,get_map,get_coord,\
                         process_colorby,create_dir

matplotlib_version = re.split("[^\d]", matplotlib.__version__)
matplotlib_version_info = tuple([int(i) for i in matplotlib_version if \
                          i.isdigit()])
                          
if matplotlib_version_info != (0,98,5,3):
     print "This code was only tested with Matplotlib-0.98.5.3"

TABLE_HTML = """<table cellpadding=0 cellspacing=0 border=0>
<tr>
<td class=normal align=center>%s</td>
<td class=normal align=center>%s</td>
<td class=normal align=center>%s</td>
</tr>
</table>""" 

PAGE_HTML = """
<html>
<head>
<style type="text/css">
.normal { color: black; font-family:Arial,Verdana; font-size:12; 
font-weight:normal;}
</style>
<script type="text/javascript" src="js/overlib.js"></script>
<title>%s</title>
</head>
<body>
<div id="overDiv" style="position:absolute; visibility:hidden; z-index:1000;">\
</div>
%s
</body>
</html>
"""

IMG_SRC = """<img src="%s" border=0 />"""
DOWNLOAD_LINK = """<a href="%s" >%s</a>"""

AREA_SRC = """<AREA shape="circle" coords="%d,%d,5" href="#%s"  \
onmouseover="return overlib('%s');" onmouseout="return nd();">\n"""
IMG_MAP_SRC = """<img src="%s" border="0" ismap usemap="#points%s" width="%d" \
height="%d" />\n"""

MAP_SRC = """
<MAP name="points%s">
%s
</MAP>
"""

shape = [ 
    's', #: square
    'o', # : circle
    '^', # : triangle up
    '>', # : triangle right
    'v', # : triangle down
    '<', # : triangle left
    'd', # : diamond
    'p', # : pentagon
    'h', # : hexagon
]

data_colors={'blue':'#0000FF','lime':'#00FF00','red':'#FF0000', \
             'aqua':'#00FFFF','fuchsia':'#FF00FF','yellow':'#FFFF00', \
             'green':'#008000','maroon':'#800000','teal':'#008080', \
             'purple':'#800080','olive':'#808000', \
             'silver':'#C0C0C0','gray':'#808080'}

default_colors=['blue','lime','red','aqua','fuchsia','yellow','green', \
               'maroon','teal','purple','olive','silver','gray']

def make_interactive_scatter(plot_label,dir_path,data_file_link,xy_coords, 
                                props, x_len=8, y_len=4, size=10,
                                draw_axes=False, generate_eps=True):
    """Write interactive plot  

    xy_coords: a dict of form {series_label:([x data], [y data], \
    [xy point label],[color])}
    """
    alpha=1.0
    my_axis=None    
    rc('font', size='8')
    rc('patch', linewidth=0)
    rc('axes', linewidth=.5,edgecolor='black')
    rc('axes', labelsize=8)
    rc('xtick', labelsize=8,color='black')
    rc('ytick', labelsize=8,color='black')

    sc_plot=draw_scatterplot(props,xy_coords,x_len,y_len,alpha,size)
    
    mtitle = props.get("title","Groups")
    x_label = props.get("xlabel","X")
    y_label = props.get("ylabel","Y")
        
    title('%s' % mtitle, fontsize='10',color='black')
    xlabel(x_label, fontsize='8',color='black')
    ylabel(y_label, fontsize='8',color='black')

    if draw_axes:
        axvline(linewidth=.5, x=0, color='black')
        axhline(linewidth=.5, y=0, color='black')
    if my_axis is not None:
        axis(my_axis)
    
    ax=sc_plot.get_axes()
    for line in ax.yaxis.get_ticklines():
        # Color the tick lines in the 2D
        line.set_color('white')

    img_name = x_label[0:2]+'_vs_'+y_label[0:2]+'_plot.png'
    savefig(dir_path + img_name, dpi=80,facecolor='white')
    
    #Create zipped eps files
    eps_link = ""
    if generate_eps:
        eps_img_name = str(x_label[0:2]+'vs'+y_label[0:2]+'plot.eps')
        savefig(dir_path + eps_img_name,format='eps')
        out = getoutput("gzip -f " + dir_path + eps_img_name)
        eps_link =  DOWNLOAD_LINK % ((data_file_link + eps_img_name + ".gz"), \
                                     "Download Figure")

    all_cids,all_xcoords,all_ycoords=transform_xy_coords(xy_coords,sc_plot)

    xmap,img_height,img_width=generate_xmap(x_len,y_len,all_cids,all_xcoords,\
                                            all_ycoords)


    points_id = plot_label+x_label[1:2]+y_label[1:2]

    return IMG_MAP_SRC % (data_file_link + img_name, points_id, img_width, \
                          img_height), MAP_SRC % (points_id, ''.join(xmap)), \
                          eps_link 

def generate_xmap(x_len,y_len,all_cids,all_xcoords,all_ycoords):
    """Generates the html interactive image map"""
    #Determine figure height and width"""
    img_height = x_len * 80
    img_width = y_len * 80

    #Write html script which allows for mouseover of labels
    xmap = []
    for cid, x, y in zip(all_cids, all_xcoords, all_ycoords):
        xmap.append(AREA_SRC % (x, img_height-y, cid, cid))
    
    return xmap,img_height,img_width

def draw_scatterplot(props,xy_coords,x_len,y_len,alpha,size):
    """Create scatterplot figure"""
    
    fig = figure(figsize=(x_len,y_len))
    
    sorted_keys = xy_coords.keys()
    scatters = {} 
    size_ct =  shape_ct = 0 
    #Iterate through coords and add points to the scatterplot
    for s_label in sorted_keys:
        s_data = xy_coords[s_label]
        if s_data[0]==[]:
            pass
        else:
            c = s_data[3]
            m =shape[shape_ct % len(shape)]
            ax = fig.add_subplot(111,axisbg='white')
            sc_plot = ax.scatter(s_data[0], s_data[1], c=c, marker=m, \
                                 alpha=alpha,s=size, linewidth=1,edgecolor=c) 
            size_ct += 1 
            shape_ct += 1
            scatters[s_label] = sc_plot
        
    return sc_plot
                          
def transform_xy_coords(xy_coords,sc_plot):
    """Transform the coords from the scatterplot into coords that can be \
       referenced in the html page"""
    sorted_keys = xy_coords.keys()
    all_cids = [] 
    all_xcoords = []
    all_ycoords = []
    sc_plot.set_transform(sc_plot.axes.transData)
    trans=sc_plot.get_transform()
    
    for s_label in sorted_keys: 
        s_data = xy_coords[s_label]
        if s_data[0]==[]:
            pass
        else:
            icoords=trans.transform(zip(s_data[0],s_data[1]))
            xcoords, ycoords = zip(*icoords)
            all_cids.extend(s_data[2])
            all_xcoords.extend(xcoords)
            all_ycoords.extend(ycoords)
    
    return all_cids,all_xcoords,all_ycoords
    
def draw_pca_graph(plot_label,dir_path,data_file_link,coord_1,coord_2,data,prefs,groups,colors,\
                   generate_eps=True):
    """Draw PCA graphs"""
    coords,pct_var=convert_coord_data_to_dict(data)
    mapping = data['map']

    if coord_1 not in coords:
        raise ValueError, "Principal coordinate: %s not available." % coord_1

    if coord_2 not in coords:
        raise ValueError, "Principal coordinate: %s not available." % coord_2

    #Handle matplotlib scale bug when all coords are 0.0
    if not len([x for x in map(float, coords[coord_2]) if x != 0.0]):
        for ix in range(len(coords[coord_2])):
            coords[coord_2][ix] = '1e-255'
    if not len([x for x in map(float, coords[coord_1]) if x != 0.0]):
         for ix in range(len(coords[coord_1])):
            coords[coord_1][ix] = '1e-255'

    #Write figure labels
    props = {}
    props["title"] = "PCA - P%s vs P%s" % (coord_1, coord_2)
    props["ylabel"] = "P%s - Percent variation explained %.2f%%" \
                        % (coord_2, float(pct_var[coord_2]))
    props["xlabel"] = "P%s - Percent variation explained %.2f%%" \
                        % (coord_1, float(pct_var[coord_1])) 

    labels = coords['pc vector number']

    p1 = map(float, coords[coord_2])
    p2 = map(float, coords[coord_1])

    if len(p1) != len(p2):
        raise ValueError, "Principal coordinate vectors unequal length."

    p1d = dict(zip(labels, p1))
    p2d = dict(zip(labels, p2))

    xy_coords=extract_and_color_xy_coords(p1d,p2d,colors,groups,coords)
 
    img_src, img_map, eps_link =  make_interactive_scatter(plot_label,dir_path,data_file_link,
                                     xy_coords=xy_coords,props=props,x_len=4.5, 
                                     y_len=4.5,size=20,draw_axes=True,
                                     generate_eps=generate_eps)

    return img_src + img_map, eps_link
    
def extract_and_color_xy_coords(p1d,p2d,colors,groups,coords):
    """Extract coords from appropriate columns and attach their \
       corresponding colors based on the group"""

    xy_coords = {}
    for group_name, ids in sorted(groups.items()):
        x=0
        color = colors[group_name]
        for id_ in sorted(ids):
            cur_labs = []
            cur_x = []
            cur_y = []
            cur_color = []
            if id_ in coords['pc vector number']:
                cur_labs.append(id_+': '+group_name)
                cur_x.append(p2d[id_])
                cur_y.append(p1d[id_])
                cur_color.append(colors[group_name])
 
            xy_coords["%s" % id_] = (cur_x, cur_y, cur_labs,cur_color)
    
    return xy_coords

def create_html_filename(coord_filename,name_ending):
    """Generate html filename using the given coord filename"""
    outpath = coord_filename.split('/')[-1] + name_ending 
    return outpath

def convert_coord_data_to_dict(data):
    """Convert the coord data into a dictionary"""
    coord_header=data['coord'][0]
    coords=data['coord'][1]
    pct_var=data['coord'][3]
    coords_dict={}
    pct_var_dict={}
    coords_dict['pc vector number']=coord_header
    for x in range(len(coords)):
        coords_dict[str(x+1)]=coords[0:,x]
        pct_var_dict[str(x+1)]=pct_var[x]

    return coords_dict,pct_var_dict 

def write_html_file(out_table,outpath):  
    """Write 2D plots into an html file"""
    page_out = PAGE_HTML % (outpath, out_table)
    out = open(outpath, "w+") 
    out.write(page_out)
    out.close()

#The following functions are not unit_tested
def _make_cmd_parser():
    """Returns the command-line options"""
    parser = OptionParser(usage="Usage: this_file.py -i <pca/pcoa output files>\
\nor\nUsage: this_file.py -i <pca/pcoa output files> -m <mapping output file>\
-b 'Mapping column to color by' -o <write to directory>")
    parser.add_option('-i', '--coord_fname', dest='coord_fname', \
        help='name of coords file [REQUIRED]')
    parser.add_option('-m', '--map_fname', dest='map_fname', \
        help='name of mapping file [default=%default]')
    parser.add_option('-b', '--colorby', dest='colorby',\
        help='map header to color by [default=%default]')
    parser.add_option('-o', '--dir_path', dest='dir_path',\
        help='directory prefix for all analyses [default=%default]',default='')
    options, args = parser.parse_args()
    return options

def _process_prefs(options):
    """opens files as necessary based on prefs"""
    data = {}
        
    #Open and get coord data
    data['coord'] = get_coord(options, data)

    #Open and get mapping data, if none supplied create a pseudo mapping
    #file
    if options.map_fname:
        mapping = get_map(options, data)
    else:
        data['map']=(([['#SampleID','Sample']]))
        for i in range(len(data['coord'][0])):
            data['map'].append([data['coord'][0][i],'Sample'])

    #Determine which mapping headers to color by, if none given, color by
    #Sample ID's
    if options.colorby:
        prefs,data=process_colorby(options.colorby,data)
    else:
        prefs={}
        prefs['Sample']={}
        prefs['Sample']['column']='#SampleID'
    
    filepath=options.coord_fname
    filename=filepath.strip().split('/')[-1]

    file_path=__file__.split('/')
    if len(file_path)==1:
        qiime_dir='./'
    else:
        qiime_dir='';
        for i in range(len(file_path)-1):
            qiime_dir+=file_path[i]+'/'
    
    
    dir_path=options.dir_path
    if not dir_path.endswith("/"):
        dir_path=dir_path+"/"
    
    dir_path=create_dir(dir_path,'2d_plots_')
    
    js_dir_path = dir_path+'/js/'
    try:
        os.mkdir(js_dir_path)
    except OSError:
        pass
    shutil.copyfile(qiime_dir+'/js/overlib.js', js_dir_path+'overlib.js')

    action_str = '_do_2d_plots'
    try:
        action = eval(action_str)
    except NameError:
        action = None
    #Place this outside try/except so we don't mask NameError in action
    if action:
        action(prefs, data, dir_path,filename)

def _do_2d_plots(prefs,data,dir_path,filename):
    """Generate interactive 2D scatterplots"""
    coord_tups = [("1", "2"), ("3", "2"), ("1", "3")]
    mapping=data['map']
    out_table=''
    #Iterate through prefs and generate html files for each colorby option
    for name, p in prefs.items():
        group_num=-1
        col_name = p['column']
        new_col_name=col_name.strip('#')
        
        alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
        alphabet += alphabet.lower()
        alphabet += "01234567890"

        data_file_path=''.join([choice(alphabet) for i in range(10)])
        data_file_path='2d_plots_'+strftime("%Y_%m_%d_%H_%M_%S")+data_file_path+new_col_name+'/'
        data_file_dir_path = dir_path+data_file_path
        data_file_link='./'+data_file_path

        data_file_dir_path=create_dir(data_file_dir_path,'')
        colors={}
        groups = group_by_field(mapping, col_name)
        
        #Iterate through groups and associate colors to each group
        for g in groups:
            if g not in colors:
                group_num+=1
                if group_num==len(data_colors):
                    group_num=0
                colors[g] = data_colors[default_colors[group_num]]
        
        img_data = {}
        plot_label=p['column']
        for coord_tup in coord_tups: 
            coord_1, coord_2 = coord_tup
            img_data[coord_tup] = draw_pca_graph(plot_label,data_file_dir_path,data_file_link,coord_1,coord_2,data,\
                                                 prefs,groups,colors,\
                                                 generate_eps=True)    

        out_table += TABLE_HTML % ("<br>".join(img_data[("1", "2")]),
                                  "<br>".join(img_data[("3", "2")]),
                                   "<br>".join(img_data[("1", "3")]))
    outfile = create_html_filename(filename+new_col_name,'_pca_2D.html')
    outfile=dir_path+outfile
        
    write_html_file(out_table,outfile)

if __name__ == "__main__":
    from sys import argv, exit
    options = _make_cmd_parser()
    
    #Kept, just in case we allow for reading a prefs file
    #prefs = eval(open(options.pref_fname, 'U').read()) 

    _process_prefs(options)

        