#!/usr/bin/env python
#file gen_pie_charts.py

__author__ = "Julia Goodrich and Michah Hamady"
__copyright__ = "Copyright 2009, Qiime" #consider project name
__credits__ = ["Julia Goodrich"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "julia.goodrich@colorado.edu"
__status__ = "Prototype"


"""
Author: Julia Goodrich (julia.goodrich@colorado.edu) and Micah Hamady
Status: Prototype

Requirements:
MatPlotLib 0.98.5.2
Python 2.5

Example: Create 2D pie charts using taxonomy counts:
Usage: python gen_pie_charts.py -c phylum.txt,class.txt,genus.txt -x ./webfiles

"""

import matplotlib,re
matplotlib.use('Agg')
from pylab import rc, axis, title, axes, pie, figlegend, clf, savefig, figure\
     ,close
from commands import getoutput
from string import strip
#from parse import parse_map,parse_coords,group_by_field,group_by_fields
from numpy import array
from optparse import OptionParser
from collections import defaultdict
from time import strftime
from random import choice, randrange
from matplotlib.font_manager import FontProperties
import os
matplotlib_version = re.split("[^\d]", matplotlib.__version__)
matplotlib_version_info = tuple([int(i) for i in matplotlib_version if \
                          i.isdigit()])

if matplotlib_version_info != (0,98,5,3):
    print "This code was only tested with Matplotlib-0.98.5.3"

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
<link rel="stylesheet" href="css/qiime_style.css" type="text/css">
<script type="text/javascript" src="js/overlib.js"></script>

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




COLORS =  ["#FF0000", "#FF6600", "#FF9900", "#FFCC00", "#FFFF99", "#FFFF00",
    "#CCFF00", "#99FF00","#00FF00", "#66CC00", "#009900", "#006666", "#66CC99",
    "#6699CC","#6666FF", "#0033FF", "#330099", "#663399", "#660066", "#9999CC",
    "#000066","#0000FF","#00CCFF","#660066","#6600CC","#9966CC","#CC0066",
    "#CC6666","#eeeeee"]


ALT_COLORS = ["#8B0000","#B22222","#FF0000", "#FF4500", "#FF8C00", "#FFA500",
              "#FFD700","#FFFF00","#F0E68C","#FFFFE0", "#006400","#008000",
              "#32CD32","#00FF00", "#7FFF00","#98FB98","#00008B", "#0000CD",
              "#0000FF","#1E90FF", "#00BFFF", "#87CEFA","#008080","#20B2AA",
              "#48D1CC","#7FFFD4", "#4B0082", "#800080","#9932CC","#9966CC",
              "#BA55D3","#FF00FF", "#C71585","#FF1493","#FF69B4","#DB7093",
              "#DDA0DD","#8B4513","#D2691E","#eeeeee"]

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
    
def make_pie_chart(data, dir_path, file_prefix = None,props={},
	y_len=6.5, colors=None, dpi=80,
    generate_eps=False, generate_pdf=True, others_key = "All Other Categories",
    others_color = "#eeeeee", should_capitalize=True):
    """
"    Write interactive piechart 

    data: [fraction:label,...] 

    trunc_len: truncates labels after this many chars

    """
    if not colors:
        colors = ALT_COLORS

    if not data:
        raise ValueError, "No data available for pie chart."

    all_fracs = []
    all_labels = []

    # set up labels and colors for pie chart 
    for color_ix, (c_frac, c_label) in enumerate(data):
        c_label = c_label.replace("_", " ")
        # we also want to color others category same every time
        if c_label == others_key:
            colors[color_ix] = others_color

        all_fracs.append(c_frac)
        if should_capitalize:
            capital = "%s (%.2f%%)"%(c_label.capitalize(),(c_frac*100.0))
            all_labels.append(capital)
        else:
            all_labels.append("%s (%.2f%%)" % (c_label, (c_frac * 100.0)))
    rc('font', size='12')
    rc('patch', linewidth=.1)
    rc('axes', linewidth=.5)
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
    title(mtitle, fontsize='14')
    ax = axes([0.0, 0.0, .5, 1])
    p1 = pie(all_fracs,  shadow=False, colors=colors)
    figlegend(p1[0],labels = all_labels, loc = loc, borderpad=0.3, \
				    labelspacing=0.3, prop = fp) 

    #write out
    if file_prefix is None:
        img_name = make_img_name()
    else:
        img_name = file_prefix
    img_abs =  os.path.join(dir_path, img_name)
    savefig(img_abs, dpi=dpi)
    eps_link = ""
    eps_abs = ""

    if generate_pdf:
        if file_prefix is None:
            eps_img_name = make_img_name(file_ext=".pdf")
        else:
            eps_img_name = file_prefix + ".pdf"
        savefig(os.path.join(dir_path, eps_img_name))
        eps_abs = os.path.join(dir_path, eps_img_name)
        eps_link = DOWNLOAD_LINK % ((os.path.join(dir_path,eps_img_name)),\
			    IMG_SRC % (os.path.join(dir_path,img_name))) 
    if generate_eps:
        if file_prefix is None:
            eps_img_name = make_img_name(file_ext=".eps")
        else:
            eps_img_name = file_prefix + ".eps"
        savefig(os.path.join(dir_path, eps_img_name))
        strip_eps_font(os.path.join(dir_path, eps_img_name))
        out = getoutput("gzip " + os.path.join(dir_path, eps_img_name))
        eps_abs = os.path.join(dir_path, eps_img_name) + ".gz"
        eps_link=DOWNLOAD_LINK % ((os.path.join(dir_path, eps_img_name)+".gz"),\
				IMG_SRC % (os.path.join(dir_path,img_name))) 
    close(fig)
    clf()
    return eps_link, IMG_SRC_2 % (os.path.join(dir_path,img_name))

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
            fracs_labels_other.append((frac,t))
        tax = s.strip().split("<br>")[-1]
        for_overlib = s.strip().rpartition("<br>")[0]
        all_counts.append(DATA_HTML % (n,frac*100,for_overlib,tax, tax,t))
    if len(counts) > num_categories:
        other_cat = len(counts) - (num_categories-1)
        new_counts = counts[0:num_categories-1]
        other = sum([c_over[0] for c_over in counts[num_categories-1:]])
        other_frac = float(other)/total
        fracs_labels = [(float(n)/red,t) for n,t,s in new_counts]

    return fracs_labels_other,fracs_labels,all_counts, other_cat, red,other_frac

def make_HTML_table(l,other_frac,total,red,other_cat,
                    fracs_labels_other,fracs_labels,dir_path,all_counts):
    """Makes the HTML table for one set of pie charts """
    img_data = []
    if other_cat > 0:
        fracs_labels_other.append((other_frac,"All Other Categories"))
        title = TITLE_include % (l,total,total,\
                        len(fracs_labels_other), total-red, other_cat)
        all_taxons = [l]
        pie = make_pie_chart(fracs_labels_other,dir_path,\
                            props = {'title':title})
        all_taxons.extend(pie)
        title = TITLE_exclude % (l,red,total,len(fracs_labels), \
                total-red, other_cat)
        pie  = make_pie_chart(fracs_labels,
                              dir_path,props = {'title':title})
        all_taxons.extend(pie)
        all_taxons.append('\n'.join(all_counts))
        img_data.append(TABLE_graph % tuple(all_taxons))
    else:
        title = TITLE % (l,total,total,len(fracs_labels_other))
        all_taxons = [l]
        
        pie = make_pie_chart(fracs_labels_other,dir_path,\
                                props = {'title':title})
        all_taxons.extend(pie)
        all_taxons.extend(("",""))
        all_taxons.append('\n'.join(all_counts))
        img_data.append(TABLE_graph % tuple(all_taxons))


    return img_data
                
def get_counts(lines, label,do_sample, num_categories, dir_path):
    """gets all the counts for one input file"""
    labels = lines[1].strip().split("Consensus Lineage")[0].split()[1:-1]
    img_data = []

    level_counts = []
    for line in lines[2:]:
        counts = line.strip().split()
        taxonomy = counts[0]
        split_label = taxonomy.strip().split(";")
        level_counts.append((sum(map(int,counts[1:])), taxonomy,
                            '<br>'.join(split_label)))

    all_sum = sum([c_over[0] for c_over in level_counts])
    fracs_labels_other,fracs_labels,all_counts, other_cat, red, other_frac= \
                                get_fracs(level_counts, num_categories, all_sum)

    img_data.extend(make_HTML_table(label,other_frac,all_sum,red,other_cat,
                    fracs_labels_other,fracs_labels,dir_path,all_counts))

    
    if do_sample:
        for i, l in enumerate(labels):
            total = 0
            sample_counts = []
            for line in lines[2:]:
                counts = line.strip().split()
                taxonomy = counts[0]
                split_label = taxonomy.strip().split(";")
                c = int(counts[i+1])
                if c > 0:
                    total += c
                    sample_counts.append((c,taxonomy,'<br>'.join(split_label)))
            fracs_labels_other,fracs_labels,all_counts,other_cat,red,other_frac= \
                                get_fracs(sample_counts, num_categories, total)
            img_data.extend(make_HTML_table('_'.join([label,l.strip()]),
                            other_frac,total,red,other_cat,fracs_labels_other,
                            fracs_labels,dir_path,all_counts))

    return img_data

usage_str = """usage: %prog [options] {-i LIST_OF_TAXONOMY_FILES> \
-l LIST_OF_LABELS}

[] indicates optional input (order unimportant) 
{} indicates required input (order unimportant) 

Requirements:
MatPlotLib 0.98.5.2
Python 2.5

Example: Create pie charts using taxonomy counts for combined samples:
python gen_pie_charts.py -i phylum.txt,class.txt,genus.txt  \
-l phylum,class,genus -o ./webfiles

Create pie charts using taxonomy counts for combined samples and individual \
samples:
python gen_pie_charts.py -i phylum.txt,class.txt,genus.txt  \
-l phylum,class,genus -o ./webfiles -s

"""

def _make_cmd_parser():
    """Returns the command-line options"""
    usage =usage_str
    version = 'Version: %prog ' +  __version__
    
    parser = OptionParser(usage=usage_str, version=version)
    
    parser.add_option('-i', '--input_files', dest='counts_fname',
                help='list of files with sample counts by taxonomy [REQUIRED]')
    parser.add_option('-l', '--labels', dest='labels',
            help='list of labels for pie chart(i.e. Phylum,Class)[REQUIRED]')
    parser.add_option('-s', '--sample_flag', dest='do_sample',
     help='if True pie charts will be created for each sample',default=False,
                      action = 'store_true')
    parser.add_option('-n', '--num', dest='num_categories', \
                help='name of file containing metadata [default: %default]', \
                      default='20')
    parser.add_option('-o', '--dir-prefix', dest='dir_path',\
	           help='directory prefix for all analyses')
    opts, args = parser.parse_args()

    if not opts.counts_fname:
        parser.error("A list of input files must be specified")

    if not opts.labels:
        parser.error("A list of label names cooresponding to files must be \
                    specified")
        
    return opts

def create_dir(dir_path,plot_type):
    """Creates directory where data is stored.  If directory is not supplied in\
       the command line, a random folder is generated"""
       
    alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
    alphabet += alphabet.lower()
    alphabet += "01234567890"

    
    if dir_path==None or dir_path=='':
        dir_path=''
        random_dir_name=''.join([choice(alphabet) for i in range(10)])
        dir_path = os.path.join(os.getcwd(),
                    plot_type+strftime("%Y_%m_%d_%H_%M_%S")+random_dir_name)

    if dir_path:
        try:
            os.mkdir(dir_path)
        except OSError:
            pass


    return dir_path


def _process_prefs(options):
    """opens files as necessary based on prefs"""
    data = []
    
    dir_path = create_dir(options.dir_path, "webfiles")

    do_sample = options.do_sample
    counts_fname = options.counts_fname
    labels = options.labels
    data = [(label,f.strip()) \
            for f,label in zip(counts_fname.split(","),labels.split(","))]
    
    filepath=data[0][1]
    filename=filepath.strip().rpartition('/')[0]
    num_categories = int(options.num_categories)
    action_str = '_do_pie_charts'
    try:
        action = eval(action_str)
    except NameError:
        action = None
    #Place this outside try/except so we don't mask NameError in action
    if action:
        action(data,dir_path,filename,num_categories, do_sample)


def _do_pie_charts(data, dir_path, filename,num_categories, do_sample):
    """Generate interactive pie charts in one HTML file"""
    img_data = []
    for label,f_name in data:
        f = open(f_name)
        lines = f.readlines()
        f.close()
        img_data.extend(get_counts(lines,label.strip(),do_sample,
                                   num_categories, dir_path))
        
    outpath = os.path.join(filename,'taxonomy_summary_pie_chart.html')
    out_table = '\n'.join(img_data)
    write_html_file(out_table,outpath)

if __name__ == "__main__":
    from sys import argv, exit
    options = _make_cmd_parser()
    _process_prefs(options)

