#!/usr/bin/env python
#file make_pie_charts.py

__author__ = "Jesse Stobmaugh"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Jesse Stobmaugh","Julia Goodrich",\
               "Micah Hamady"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"


"""
Requirements:
MatPlotLib 0.98.5.2
Python 2.5

Example: Create 2D charts using taxonomy counts:
Usage: python plot_taxa_summary.py -c phylum.txt,class.txt,genus.txt -x ./webfiles

"""

import matplotlib,re
matplotlib.use('Agg',warn=False)
from matplotlib.font_manager import FontProperties
from pylab import rc, axis, title, axes, pie, figlegend, clf, savefig, figure\
     ,close
from commands import getoutput
from string import strip
from numpy import array
import numpy as numpy
from optparse import OptionParser
from collections import defaultdict
from time import strftime
from random import choice, randrange
import os
import shutil
from qiime.util import get_qiime_project_dir
from qiime.colors import natsort, data_color_order, data_colors, \
                            get_group_colors,iter_color_groups
from qiime.parse import group_by_field

ALPHABET = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
ALPHABET += ALPHABET.lower()
ALPHABET += "01234567890"

'''define some html elements that we can populate with data'''
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
'Chart', CAPCOLOR, '#000000', BGCOLOR, '#DDDDDD', CLOSECOLOR, 'blue');" \
onmouseout="return nd();" >%s</a><br></td>
<td class="ntitle">&nbsp;</td>
</tr>
<tr><td class="ntitle"><a href="#" target="_blank" onmouseover="return \
overlib('%s', STICKY, MOUSEOFF, FGCOLOR, '#FFFFFF', ABOVE, LEFT, CAPTION, \
'Chart', CAPCOLOR, '#000000', BGCOLOR, '#DDDDDD', CLOSECOLOR, 'blue');" \
onmouseout="return nd();" >%s</a><br></td>
</tr>
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

DATA_TABLE_HTML="""
<table cellpadding=2 cellspacing=2 border=0 \
width=600><tr class=ntitle> <td valign=bottom class=header>Count</td> \
<td valign=bottom  class=header nowrap>Pct</td> <td valign=bottom class=header \
nowrap>Taxonomy</td></tr>
%s
</table>"""

DATA_HTML = """<tr class=normal><td>%s</td> <td nowrap>%.2f%%</td>\
<td class="normal" ><a onmouseover="return overlib('<b>Taxonomy:</b><br>%s<br>\
<a href=javascript:gg(\\'%s\\');>%s</a> ',STICKY,MOUSEOFF,RIGHT);" \
onmouseout="return nd();">%s</a></td></tr>"""

IMG_SRC = """<img src=\\'%s\\' border=1 />"""

IMG_SRC_2 = """<img src='%s' height=200 width=400 border=1 />"""

DOWNLOAD_LINK = """<a href=\\'%s\\' >%s</a>"""


def strip_eps_font(filename):
    """ Rewrite file by changing font to handle eps bug"""
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
    
def make_pie_chart(data, dir_path,level,prefs,pref_colors,background_color,\
                  label_color, file_prefix = None,props={}, y_len=6.5, dpi=80,\
                  generate_eps=False, generate_pdf=True,\
                  others_key = "All Other Categories",\
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

    # set up labels and colors for pie chart 
    for color_ix, (c_label, c_frac) in enumerate(data):
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
    
    #define figure parameters
    rc('font', size='10')
    rc('text', color=label_color)
    rc('patch', linewidth=.1)
    rc('axes', linewidth=.5,edgecolor=label_color)
    rc('text', usetex=False)
    
    #generate figure object
    fig = figure(randrange(10000), figsize=(2*y_len,y_len))
    
    fp = FontProperties()
    fp.set_size('8')
    if len(data) > 30:
        loc = 4
    else:
        loc = 5
    
    #define title and legend
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
    img_abs =  os.path.join(dir_path,'charts', img_name)
    savefig(img_abs, dpi=dpi,facecolor=background_color)
    eps_link = ""
    eps_abs = ""

    #creates the plot as pdf
    if generate_pdf:
        if file_prefix is None:
            eps_img_name = make_img_name(file_ext=".pdf")
        else:
            eps_img_name = file_prefix + ".pdf"
        savefig(os.path.join(dir_path,'charts', eps_img_name),\
                             facecolor=background_color)
        eps_abs = os.path.join('charts', eps_img_name)
        eps_link = DOWNLOAD_LINK % ((os.path.join('charts',\
                eps_img_name)),\
        IMG_SRC % (os.path.join('charts',img_name))) 

    #creates the image as an eps
    if generate_eps:
        if file_prefix is None:
            eps_img_name = make_img_name(file_ext=".eps")
        else:
            eps_img_name = file_prefix + ".eps"
        savefig(os.path.join(dir_path,'charts', eps_img_name),\
                facecolor=background_color)
        strip_eps_font(os.path.join(dir_path,'charts', eps_img_name))
        out = getoutput("gzip " + os.path.join(dir_path,'charts', eps_img_name))
        eps_abs = os.path.join(dir_path,'charts',eps_img_name) + ".gz"
        eps_link=DOWNLOAD_LINK % ((os.path.join('charts', eps_img_name)+".gz"),\
        IMG_SRC % (os.path.join('charts',img_name))) 
    close(fig)
    clf()
    
    return eps_link, IMG_SRC_2 % (os.path.join('charts',img_name))


def make_area_chart(sample_ids,taxa_percents,taxa,dir_path,level,prefs,\
                    pref_colors,\
                    background_color,label_color,\
                    file_prefix = None,props={},y_len=6.5,\
                    dpi=80,generate_eps=False, generate_pdf=True,\
                    others_key = "All Other Categories",\
                    others_color = "#eeeeee", should_capitalize=True):
    """
    Write interactive area chart 
    data: [fraction:label,...] 
    trunc_len: truncates labels after this many chars
    
    This function is very similar to the make_pie_charts function
    """
    
    if not taxa_percents:
        raise ValueError, "No data available for area chart."

    all_fracs = []
    all_labels = []
    colors = []
    
    #define figure parameters
    rc('font', size='10')
    rc('text', color=label_color)
    rc('patch', linewidth=.1)
    rc('axes', linewidth=.5,edgecolor=label_color)
    rc('text', usetex=False)
    rc('xtick', labelsize=8,color=label_color)
    fig = figure(randrange(10000), figsize=(2*y_len,y_len))
    ax1 = fig.add_subplot(111)
    
    fp = FontProperties()
    fp.set_size('8')
    if len(sample_ids) > 30:
        loc = 4
    else:
        loc = 5
    
    x = numpy.arange(0, len(sample_ids))
    y_data = numpy.row_stack((zip(*taxa_percents)))

    y_data_stacked = numpy.cumsum(y_data, axis=0)

    ax1.fill_between(x, 0, y_data_stacked[0,:],\
                     facecolor=data_colors[pref_colors[taxa[0]]].toHex(),\
                     alpha=1)

    for i,j in enumerate(y_data_stacked):
        if i < len(y_data_stacked)-1:
            next=i+1
            ax1.fill_between(x, y_data_stacked[i,:], y_data_stacked[next,:],\
                        facecolor=data_colors[pref_colors[taxa[i+1]]].toHex(),\
                        alpha=1)
        else:
            ax1.fill_between(x, y_data_stacked[i,:], 1,\
                        facecolor=data_colors[pref_colors[taxa[i+1]]].toHex(),\
                             alpha=1)   

    if "title" in props:
        mtitle = props["title"]

    ax1.xaxis.set_ticks(x)
    ax1.set_xticklabels(sample_ids,rotation='vertical')
    ax1.set_yticks([])
   
    #title(mtitle, fontsize='10',color=label_color)

    if file_prefix is None:
        img_name = make_img_name()
    else:
        img_name = file_prefix
    img_abs =  os.path.join(dir_path,'charts', img_name)
    savefig(img_abs, dpi=dpi,facecolor=background_color)
    eps_link = ""
    eps_abs = ""
    
    if generate_pdf:
        if file_prefix is None:
            eps_img_name = make_img_name(file_ext=".pdf")
        else:
            eps_img_name = file_prefix + ".pdf"
        savefig(os.path.join(dir_path,'charts', eps_img_name),\
                             facecolor=background_color)
        eps_abs = os.path.join('charts', eps_img_name)
        eps_link = DOWNLOAD_LINK % ((os.path.join('charts',\
                eps_img_name)),\
        IMG_SRC % (os.path.join('charts',img_name))) 

    if generate_eps:
        if file_prefix is None:
            eps_img_name = make_img_name(file_ext=".eps")
        else:
            eps_img_name = file_prefix + ".eps"
        savefig(os.path.join(dir_path,'charts', eps_img_name),\
                facecolor=background_color)
        strip_eps_font(os.path.join(dir_path,'charts', eps_img_name))
        out = getoutput("gzip " + os.path.join(dir_path,'charts', eps_img_name))
        eps_abs = os.path.join(dir_path,'charts',eps_img_name) + ".gz"
        eps_link=DOWNLOAD_LINK % ((os.path.join('charts', eps_img_name)+".gz"),\
        IMG_SRC % (os.path.join('charts',img_name))) 
    close(fig)
    clf()
    
    return eps_link, IMG_SRC_2 % (os.path.join('charts',img_name))

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

def get_fracs(counts, num_categories, total, chart_type,sort_data=True):
    """"Returns the fractions to be used in matplotlib piechart"""
    fracs_labels_other = []
    fracs_labels = []
    all_counts = []
    other_cat = 0
    other_frac = 0
    red = 0
    
    # added in the ability to turn off sorting, since we want the data to be
    # unsorted for the area charts
    if sort_data:
        counts.sort()
        counts.reverse()
    
    area_table_out=[]
    
    # this loop iterates over the OTU table and generates html code for the 
    # data table
    for j,(n, t, s) in enumerate(counts):
        frac = float(n)/total
        if j < num_categories-1:
            red += n
            fracs_labels_other.append((t,frac))
        tax = s.strip().split("<br>")[-1]
        tax = tax.replace('"','')
        for_overlib = s.strip().rpartition("<br>")[0]
        for_overlib = for_overlib.replace('"','')

        ###Added this code because the data table is being presented
        ###differently for the area charts
        if chart_type=='pie':
            all_counts.append(DATA_HTML % (n,frac*100,for_overlib,tax, tax,t))
        elif chart_type=='area':
            area_table_out.append(str(n))

    #returning a dictionary for the case of area charts, which is different
    #than the array passed by the pie charts
    if chart_type=='area':
        all_counts=area_table_out
        
    if len(counts) > num_categories:
        other_cat = len(counts) - (num_categories-1)
        new_counts = counts[0:num_categories-1]
        other = sum([c_over[0] for c_over in counts[num_categories-1:]])
        other_frac = float(other)/total
        fracs_labels = [(t,float(n)/red) for n,t,s in new_counts]

    # added in the ability to turn off sorting, since we want the data to be
    # unsorted for the area charts
    if sort_data:
        fracs_labels_other.sort()
        fracs_labels.sort()

    return fracs_labels_other,fracs_labels,all_counts,other_cat,red,other_frac

def make_HTML_table(l,other_frac,total,red,other_cat,fracs_labels_other,\
                    fracs_labels,dir_path,all_counts,level,\
                    prefs,pref_colors,background_color, label_color,chart_type,\
                    label):
    """Makes the HTML table for one set of pie charts """
    img_data = []

    if chart_type=='pie':
        if other_cat > 0:
        
            fracs_labels_other.append(("All Other Categories",other_frac))
            title = TITLE_include % (l,total,total,\
                            len(fracs_labels_other), total-red, other_cat)
            all_taxons = [l]
        
            pie = make_pie_chart(fracs_labels_other,dir_path,level,\
                                prefs,pref_colors,background_color,label_color,\
                                props = {'title':title})
            all_taxons.extend(pie)
        
            title = TITLE_exclude % (l,red,total,len(fracs_labels), \
                    total-red, other_cat)
            pie  = make_pie_chart(fracs_labels,dir_path,level,\
                                prefs,pref_colors,background_color,label_color,\
                                props = {'title':title})
            all_taxons.extend(pie)

            img_data.append(TABLE_graph % tuple(all_taxons))
            img_data.append(DATA_TABLE_HTML % '\n'.join(all_counts))

        else:
            title = TITLE % (l,total,total,len(fracs_labels_other))
            all_taxons = [l]
        
            pie = make_pie_chart(fracs_labels_other,dir_path,level,\
                            prefs,pref_colors,background_color,label_color,\
                            props = {'title':title})

            all_taxons.extend(pie)
            all_taxons.extend(("",""))
            img_data.append(TABLE_graph % tuple(all_taxons))
            img_data.append(DATA_TABLE_HTML % '\n'.join(all_counts))

    elif chart_type=='area':
        taxa_percents=fracs_labels_other
        sample_ids=l
        taxa=other_cat
        all_categories=[]
        title = TITLE % (label,total,total,len(fracs_labels_other))
        all_taxons = [label]
        area = make_area_chart(sample_ids,taxa_percents,taxa,dir_path,\
                               level,prefs,pref_colors,\
                               background_color,label_color,\
                               props = {'title':title})
        all_taxons.extend(area)
        all_taxons.extend(("",""))
        img_data.append(TABLE_graph % tuple(all_taxons))
            
    return img_data
                
def get_counts(label,colorby,num_categories,dir_path,level,color_data,\
               prefs,pref_colors,background_color,label_color,chart_type):
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
                                get_fracs(level_counts,num_categories,all_sum,\
                                chart_type,True)
    if chart_type=='pie':
        img_data.extend(make_HTML_table(label,other_frac,all_sum,red,other_cat,\
                                fracs_labels_other,fracs_labels,dir_path,\
                                all_counts,level,prefs,pref_colors,\
                                background_color,label_color,chart_type,\
                                label))

        if colorby is not None:
            for i, l in enumerate(sample_ids):
                if l not in colorby:
                    continue
                total = 0
                sample_counts = []
                for idx, counts in enumerate(otu_table):
                    taxonomy =otu_ids[idx]
                    split_label = [j for j in taxonomy.strip().split(";")]
                    taxonomy = ';'.join(split_label)

                    c = float(counts[i])
                    if c > 0:
                        total += c
                        sample_counts.append((c,taxonomy,\
                                                '<br>'.join(split_label)))
            
                fracs_labels_other,fracs_labels,all_counts,\
                other_cat,red,other_frac= get_fracs(sample_counts,\
                                                    num_categories,\
                                                    total,chart_type,True)

                img_data.extend(make_HTML_table('_'.join([label,l.strip()]),
                            other_frac,total,red,other_cat,fracs_labels_other,
                            fracs_labels,dir_path,all_counts,level,\
                            prefs,pref_colors,background_color,label_color,\
                            chart_type,l.strip()))
    
    elif chart_type=='area':
        area_plot_arr=[]
        area_plot_sample_ids=[]
        area_plot_taxa_arr=[]
        taxa_html=[]
        total_area_table_out=[]
        if colorby is not None:
            for i, l in enumerate(sample_ids):
                if l not in colorby:
                    continue
                total = 0
                area_plot_sample_ids.append(l)
                sample_counts = []
                for idx, counts in enumerate(otu_table):
                    taxonomy =otu_ids[idx]
                    split_label = [j for j in taxonomy.strip().split(";")]
                    taxonomy = ';'.join(split_label)

                    c = float(counts[i])
                    total += c
                    sample_counts.append((c,taxonomy,'<br>'.join(split_label)))
                    
                fracs_labels_other,fracs_labels,all_counts,\
                other_cat,red,other_frac=get_fracs(sample_counts,\
                                                len(sample_counts),total,\
                                                chart_type,False)
                
                total_area_table_out.append(all_counts)
            
                area_plot_per=[]
                area_plot_taxa=[]
                for i in fracs_labels_other:
                    area_plot_per.append(i[1])
                    area_plot_taxa.append(i[0])
                area_plot_arr.append(area_plot_per)
                area_plot_taxa_arr.append(area_plot_taxa)
        
        taxa_html.append('<tr><th>'+l.strip()+\
                         '</th></tr>'+''.join(all_counts)+'')
        
        data_table=zip(*total_area_table_out)
        data_html_str='<table cellpadding=2 cellspacing=2 border=0 ' + \
                      'width=600><tr class=ntitle><td valign=bottom ' + \
                      'class=header>Taxonomy</td>\n'
        
        for i in area_plot_sample_ids:
            data_html_str+='<td valign=bottom class=header>%s</td>\n' % (i)
        data_html_str+='</tr>'
        for ct,dat in enumerate(otu_ids):
            tax=dat
            split_label = [i for i in tax.strip().split(";")]
            joined_label='<br>'.join(split_label)
            data_html_str+="<tr><td bgcolor=\"%s\" class=\"normal\" ><a onmouseover=\"return  overlib(\'<b>Taxonomy:</b><br>%s<br><a href=javascript:gg(\\'%s\\');>%s</a> \',STICKY,MOUSEOFF,RIGHT);\" onmouseout=\"return nd();\">%s</a></td>" % \
            (data_colors[pref_colors[tax]].toHex(),\
             joined_label.replace('"',''),split_label[-1].replace('"',''),\
             split_label[-1].replace('"',''),tax.replace('"',''))
            
            
            for per_tax in data_table[ct]:
                if float(per_tax)>0:
                    data_html_str+='<td bgcolor=\"%s\">%5.2f</td>\n' % \
                    (data_colors[pref_colors[tax]].toHex(),(float(per_tax)*100))
                else:
                    data_html_str+='<td>%5.2f</td>\n' % \
                    ((float(per_tax)*100))
            data_html_str+='</tr>\n'
        
        data_html_str+='</table>'
        #print area_plot_taxa_arr
        for i in range(len(area_plot_taxa_arr)-1):
            if area_plot_taxa_arr[i]!=area_plot_taxa_arr[i+1]:
                raise ValueError, 'The taxonomies are out of order!'

        img_data.extend(make_HTML_table(area_plot_sample_ids,
                    other_frac,all_sum,red,otu_ids,area_plot_arr,
                    fracs_labels,dir_path,[' '.join(taxa_html)],level,\
                    prefs,pref_colors,background_color,label_color,chart_type,\
                    label))
        img_data.append(data_html_str)
        
    return img_data

def make_all_charts(data,dir_path,filename,num_categories,colorby,args,\
                        color_data, prefs,background_color,label_color,
                        chart_type):
    """Generate interactive pie charts in one HTML file"""
    #iterate over the preferences and assign colors according to taxonomy
    
    img_data = []
    for label,f_name in data:
        f = color_data['counts'][f_name]
        level = max([len(t.split(';')) - 1 for t in f[1]])
        
        prefs=prefs

        for key in prefs.keys():
            if prefs[key]['column'] != str(level):
                continue
            col_name = 'Taxon'
            mapping = [['Taxon']]
            mapping.extend([[m] for m in f[1]])
            if 'colors' in prefs[key]:
                if isinstance(prefs[key]['colors'], dict):
                    pref_colors = prefs[key]['colors'].copy() 
                    #copy so we can mutate
                else:
                    pref_colors = prefs[key]['colors'][:]
            else:
                pref_colors={}
            labelname=prefs[key]['column']

            #Define groups and associate appropriate colors to each group
            groups = group_by_field(mapping, col_name)
            pref_colors, data_colors, data_color_order = \
                get_group_colors(groups, pref_colors)
        
        img_data.extend(get_counts(label.strip(),colorby,num_categories,\
                        dir_path,level,f,prefs,pref_colors,background_color,\
                        label_color,chart_type))

    outpath = os.path.join(dir_path,'taxonomy_%s_summary_chart.html' % \
                                                                    chart_type)
    out_table = '\n'.join(img_data)
    write_html_file(out_table,outpath)
