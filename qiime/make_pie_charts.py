#!/usr/bin/env python
#file make_pie_charts.py

__author__ = "Julia Goodrich and Micah Hamady"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Julia Goodrich"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ ="Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__status__ = "Development"


"""
Author: Julia Goodrich (julia.goodrich@colorado.edu) and Micah Hamady
Status: Prototype

Requirements:
MatPlotLib 0.98.5.2
Python 2.5

Example: Create 2D pie charts using taxonomy counts:
Usage: python make_pie_charts.py -c phylum.txt,class.txt,genus.txt -x ./webfiles

"""

import matplotlib,re
matplotlib.use('Agg',warn=False)
from pylab import rc, axis, title, axes, pie, figlegend, clf, savefig, figure\
     ,close
from commands import getoutput
from string import strip
from numpy import array
from optparse import OptionParser
from collections import defaultdict
from time import strftime
from random import choice, randrange
from matplotlib.font_manager import FontProperties
import os
import shutil
from qiime.util import get_qiime_project_dir
from qiime.colors import natsort, data_color_order, data_colors, \
                            get_group_colors,iter_color_groups
from qiime.parse import group_by_field
ALPHABET = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
ALPHABET += ALPHABET.lower()
ALPHABET += "01234567890"

TITLE = "%s: %s (of %s) from %s categories displayed"
TITLE_include = """%s: %s (of %s) from %s categories displayed, including %s \
from %s categories ('All Other Categories')"""
TITLE_exclude = """%s: %s (of %s) from %s categories displayed, excluding %s \
from %s categories ('All Other Categories')"""
TABLE_graph = """<table cellpadding=2 cellspacing=2 border=0 width=800>
<tr><td class="normal" colspan=2>&nbsp;</td></tr>
<tr><td class="header" colspan=2>Taxonomy Summary. Current Level: %s</td></tr>
<tr><td class="ntitle"><a href="#" target="_blank" onmouseover="return \
overlib('%s', STICKY, MOUSEOFF, FGCOLOR, '#FFFFFF', ABOVE, LEFT, CAPTION, \
'Pie Chart', CAPCOLOR, '#000000', BGCOLOR, '#DDDDDD', CLOSECOLOR, 'blue');" \
onmouseout="return nd();" >%s</a><br></td>
<td class="ntitle">&nbsp;</td>
</tr>
<tr><td class="ntitle"><a href="#" target="_blank" onmouseover="return \
overlib('%s', STICKY, MOUSEOFF, FGCOLOR, '#FFFFFF', ABOVE, LEFT, CAPTION, \
'Pie Chart', CAPCOLOR, '#000000', BGCOLOR, '#DDDDDD', CLOSECOLOR, 'blue');" \
onmouseout="return nd();" >%s</a><br></td>
</tr>
<tr><td class="normal" colspan=2><table cellpadding=2 cellspacing=2 border=0 \
width=600><tr class=ntitle> <td valign=bottom class=header>Count</td> \
<td valign=bottom  class=header nowrap>Pct</td> <td valign=bottom class=header \
nowrap>Taxonomy</td></tr>
%s
</table>
</table>
"""

PAGE_HTML = """
<html>
<head>
<link rel="stylesheet" href="./css/qiime_style.css" type="text/css">
<script type="text/javascript" src="./js/overlib.js"></script>

<script type="text/javascript">
<!-- Begin
function set_target(new_target)
{
    sf = document.getElementById("search_form");
    sf.target = new_target;
}
function gg(targetq)
{
        window.open("http://www.google.com/search?q=" + targetq, 'searchwin');
}


//  End -->
</script>
<title>%s</title>
</head>
<body>
<div id="overDiv" style="position:absolute; visibility:hidden; z-index:1000;">\
</div>
%s
</body>
</html>
"""

DATA_HTML = """<tr class=normal><td>%s</td> <td nowrap>%.2f%%</td>\
<td class="normal" ><a onmouseover="return overlib('<b>Taxonomy:</b><br>%s<br>\
<a href=javascript:gg(\\'%s\\');>%s</a> ',STICKY,MOUSEOFF,RIGHT);" \
onmouseout="return nd();">%s</a></td></tr>"""

IMG_SRC = """<img src=\\'%s\\' border=1 />"""

IMG_SRC_2 = """<img src='%s' height=200 width=400 border=1 />"""

DOWNLOAD_LINK = """<a href=\\'%s\\' >%s</a>"""


def strip_eps_font(filename):
    """ Rewrite file """
    inf = open(filename)
    filecache = []
    in_ttf = False
    for line in inf:
        if "Bitstream" in line:
            line = line.replace("BitstreamVeraSans-Roman", "Arial")
        if line.startswith("""%%BeginFont"""):
            in_ttf = True
        if line.startswith("""%%EndFont"""):
            in_ttf = False
            continue 
        if in_ttf:
            continue
        else:
            filecache.append(line)
    inf.close()
    ouf = open(filename, "w+")
    ouf.write(''.join(filecache))
    ouf.close()
    
def make_pie_chart(data, dir_path,level,color_data,prefs,background_color,label_color, file_prefix = None,props={},
    y_len=6.5, dpi=80, generate_eps=False, generate_pdf=True, 
    others_key = "All Other Categories",
    others_color = "#eeeeee", should_capitalize=True):
    """
    Write interactive piechart 

    data: [fraction:label,...] 

    trunc_len: truncates labels after this many chars

    """
    if not data:
        raise ValueError, "No data available for pie chart."

    all_fracs = []
    all_labels = []
    colors = []
    for key in prefs.keys():
        if prefs[key]['column'] != str(level):
            continue
        col_name = 'Taxon'
        mapping = [['Taxon']]
        mapping.extend([[m] for m in color_data[1]])
        if 'colors' in prefs[key]:
            if isinstance(prefs[key]['colors'], dict):
                pref_colors = prefs[key]['colors'].copy()    #copy so we can mutate
            else:
                pref_colors = prefs[key]['colors'][:]
        else:
            pref_colors={}
        labelname=prefs[key]['column']

        #Define groups and associate appropriate colors to each group
        groups = group_by_field(mapping, col_name)
        pref_colors, data_colors, data_color_order = \
            get_group_colors(groups, pref_colors)

    # set up labels and colors for pie chart 
    for color_ix, (c_label, c_frac) in enumerate(data):
        #commented out the following line, since the key becomes invalid when
        #replacing part of the string.
        #c_label = c_label.replace("_", " ")
        # we also want to color others category same every time
        if c_label == others_key:
            colors.append(others_color)
        else: 
            colors.append(data_colors[pref_colors[c_label]].toHex())
        all_fracs.append(c_frac)
        if should_capitalize:
            capital = "%s (%.2f%%)"%(c_label.capitalize(),(c_frac*100.0))
            all_labels.append(capital)
        else:
            all_labels.append("%s (%.2f%%)" % (c_label, (c_frac * 100.0)))
    rc('font', size='10')
    rc('text', color=label_color)
    rc('patch', linewidth=.1)
    rc('axes', linewidth=.5,edgecolor=label_color)
    rc('text', usetex=False)
    fig = figure(randrange(10000), figsize=(2*y_len,y_len))

    fp = FontProperties()
    fp.set_size('8')
    if len(data) > 30:
        loc = 4
    else:
        loc = 5
    mtitle = "Pie Chart"
    if "title" in props:
        mtitle = props["title"]
    axis('off')
    title(mtitle, fontsize='10',color=label_color)
    ax = axes([0.0, 0.0, .5, 1])
    p1 = pie(all_fracs,  shadow=False, colors=colors)
    flg = figlegend(p1[0],labels = all_labels, loc = loc, borderpad=0.3, \
                 labelspacing=0.3, prop = fp) 
    flg.legendPatch.set_alpha(0.0)
    #write out
    if file_prefix is None:
        img_name = make_img_name()
    else:
        img_name = file_prefix
    img_abs =  os.path.join(dir_path,'pie_charts', img_name)
    savefig(img_abs, dpi=dpi,facecolor=background_color)
    eps_link = ""
    eps_abs = ""

    if generate_pdf:
        if file_prefix is None:
            eps_img_name = make_img_name(file_ext=".pdf")
        else:
            eps_img_name = file_prefix + ".pdf"
        savefig(os.path.join(dir_path,'pie_charts', eps_img_name),facecolor=background_color)
        eps_abs = os.path.join('pie_charts', eps_img_name)
        eps_link = DOWNLOAD_LINK % ((os.path.join('pie_charts',\
                eps_img_name)),\
        IMG_SRC % (os.path.join('pie_charts',img_name))) 
    if generate_eps:
        if file_prefix is None:
            eps_img_name = make_img_name(file_ext=".eps")
        else:
            eps_img_name = file_prefix + ".eps"
        savefig(os.path.join(dir_path,'pie_charts', eps_img_name),facecolor=background_color)
        strip_eps_font(os.path.join(dir_path,'pie_charts', eps_img_name))
        out = getoutput("gzip " + os.path.join(dir_path,'pie_charts', eps_img_name))
        eps_abs = os.path.join(dir_path,'pie_charts',eps_img_name) + ".gz"
        eps_link=DOWNLOAD_LINK % ((os.path.join('pie_charts', eps_img_name)+".gz"),\
        IMG_SRC % (os.path.join('pie_charts',img_name))) 
    close(fig)
    clf()
    return eps_link, IMG_SRC_2 % (os.path.join('pie_charts',img_name))

def make_img_name(file_ext='.png'):
    """ Generate a random file name """
    fn = []
    # format seqs and write out to temp file
    for i in range(0,30):
        fn.append(choice(ALPHABET))
    return ''.join(fn) + file_ext

def write_html_file(out_table,outpath):
    """Write pie charts into an html file"""
    page_out = PAGE_HTML % (outpath, out_table)
    out = open(outpath, "w+")
    out.write(page_out)
    out.close()


def get_fracs(counts, num_categories, total):
    """"Returns the fractions to be used in matplotlib piechart"""
    fracs_labels_other = []
    fracs_labels = []
    all_counts = []
    other_cat = 0
    other_frac = 0
    red = 0
    counts.sort()
    counts.reverse()
    for j,(n, t, s) in enumerate(counts):
        frac = float(n)/total
        if j < num_categories-1:
            red += n
            fracs_labels_other.append((t,frac))
        tax = s.strip().split("<br>")[-1]
        tax = tax.replace('"','')
        for_overlib = s.strip().rpartition("<br>")[0]
        for_overlib = for_overlib.replace('"','')
        all_counts.append(DATA_HTML % (n,frac*100,for_overlib,tax, tax,t))
    if len(counts) > num_categories:
        other_cat = len(counts) - (num_categories-1)
        new_counts = counts[0:num_categories-1]
        other = sum([c_over[0] for c_over in counts[num_categories-1:]])
        other_frac = float(other)/total
        fracs_labels = [(t,float(n)/red) for n,t,s in new_counts]

    fracs_labels_other.sort()
    fracs_labels.sort()

    return fracs_labels_other,fracs_labels,all_counts, other_cat, red,other_frac

def make_HTML_table(l,other_frac,total,red,other_cat,
                    fracs_labels_other,fracs_labels,dir_path,all_counts,level,color_data, color_prefs,background_color, label_color):
    """Makes the HTML table for one set of pie charts """
    img_data = []
    if other_cat > 0:
        fracs_labels_other.append(("All Other Categories",other_frac))
        title = TITLE_include % (l,total,total,\
                        len(fracs_labels_other), total-red, other_cat)
        all_taxons = [l]
        pie = make_pie_chart(fracs_labels_other,dir_path,level,color_data, color_prefs,background_color,label_color,\
                            props = {'title':title})
        all_taxons.extend(pie)
        title = TITLE_exclude % (l,red,total,len(fracs_labels), \
                total-red, other_cat)
        pie  = make_pie_chart(fracs_labels,
                              dir_path,level,color_data, color_prefs,background_color,label_color,props = {'title':title})
        all_taxons.extend(pie)
        all_taxons.append('\n'.join(all_counts))
        img_data.append(TABLE_graph % tuple(all_taxons))
    else:
        title = TITLE % (l,total,total,len(fracs_labels_other))
        all_taxons = [l]
        
        pie = make_pie_chart(fracs_labels_other,dir_path,level,color_data, color_prefs,background_color,label_color,\
                                props = {'title':title})
        all_taxons.extend(pie)
        all_taxons.extend(("",""))
        all_taxons.append('\n'.join(all_counts))
        img_data.append(TABLE_graph % tuple(all_taxons))


    return img_data
                
def get_counts(label,colorby, num_categories, dir_path,level,color_data, color_prefs,background_color,label_color):
    """gets all the counts for one input file"""
    img_data = []
    labels = []
    level_counts = []
    
    sample_ids, otu_ids, otu_table, lineages = color_data
    labels = sample_ids

    for idx, counts in enumerate(otu_table):
        taxonomy = otu_ids[idx]
        split_label = [i for i in taxonomy.strip().split(";")]
        taxonomy = ';'.join(split_label)
        level_counts.append((sum(map(float,counts)), taxonomy,
                            '<br>'.join(split_label)))
    all_sum = sum([c_over[0] for c_over in level_counts])
    fracs_labels_other,fracs_labels,all_counts, other_cat, red, other_frac= \
                                get_fracs(level_counts, num_categories, all_sum)

    img_data.extend(make_HTML_table(label,other_frac,all_sum,red,other_cat,
                    fracs_labels_other,fracs_labels,dir_path,all_counts,level,color_data, color_prefs,background_color,label_color))

    if colorby is not None:
        for i, l in enumerate(colorby):
            total = 0
            sample_counts = []
            for idx, counts in enumerate(otu_table):
                taxonomy =otu_ids[idx]
                split_label = [j for j in taxonomy.strip().split(";")]
                taxonomy = ';'.join(split_label)

                c = float(counts[i])
                if c > 0:
                    total += c
                    sample_counts.append((c,taxonomy,'<br>'.join(split_label)))
            fracs_labels_other,fracs_labels,all_counts,other_cat,red,other_frac= \
                                get_fracs(sample_counts, num_categories, total)
            img_data.extend(make_HTML_table('_'.join([label,l.strip()]),
                            other_frac,total,red,other_cat,fracs_labels_other,
                            fracs_labels,dir_path,all_counts,level,color_data, color_prefs,background_color,label_color))

    return img_data


def make_all_pie_charts(data, dir_path, filename,num_categories, colorby,args,color_data, color_prefs,background_color,label_color):
    """Generate interactive pie charts in one HTML file"""
    img_data = []
    
    for label,f_name in data:
        f = color_data['counts'][f_name]
        level = max([len(t.split(';')) - 1 for t in f[1]])
        img_data.extend(get_counts(label.strip(),colorby,num_categories, 
                 dir_path,level,f, color_prefs,background_color,label_color))
        
    outpath = os.path.join(dir_path,'taxonomy_summary_pie_chart.html')
    out_table = '\n'.join(img_data)
    write_html_file(out_table,outpath)

