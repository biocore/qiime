from __future__ import division
#!/usr/bin/env python
# file plot_taxa_summary.py

__author__ = "Jesse Stobmaugh"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stobmaugh", "Julia Goodrich",
               "Micah Hamady"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

"""
Requirements:
MatPlotLib 1.1.0
Python 2.7
"""

import matplotlib

matplotlib.use('Agg', warn=False)
from matplotlib.font_manager import FontProperties
from matplotlib.pyplot import (rc, axis, title, axes, pie, clf,
                               savefig, figure, close)
import matplotlib.patches as mpatches
from itertools import cycle
import numpy as numpy
from random import choice, randrange
import os
import shutil
from qiime.colors import (data_colors, get_group_colors)
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

'''the following is for the area and bar charts'''
TABLE_graph = """<table cellpadding=2 cellspacing=2 border=0>
<tr><td class="normal" colspan=2>&nbsp;</td></tr>
<tr><td class="header" colspan=2>Taxonomy Summary. Current Level: %s</td></tr>
<tr><td class="ntitle">&nbsp;&nbsp;%s&nbsp;&nbsp;%s<br></td>
<td class="ntitle">&nbsp;</td>
</tr>
<tr><td class="ntitle">%s<br></td>
</tr>
</table>
%s
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

DATA_TABLE_HTML = """
<table cellpadding=2 cellspacing=2 border=0><tr class=ntitle> \
<td valign=bottom  class=header nowrap>Pct</td> <td valign=bottom class=header \
nowrap>Taxonomy</td></tr>
%s
</table>"""

DATA_HTML = """<tr class=normal><td nowrap>%.2f%%</td>\
<td class="normal" ><a onmouseover="return overlib('<b>Taxonomy:</b><br>%s<br>\
<a href=javascript:gg(\\'%s\\');>%s</a> ',STICKY,MOUSEOFF,RIGHT);" \
onmouseout="return nd();">%s</a></td></tr>"""

AREA_SRC = """<AREA shape="rect" coords="%.4f%%,%.4f%%,%.4f%%,%.4f%%" href="#%s"  \
onmouseover="return overlib('%s&nbsp;(%.2f%%)');" onmouseout="return nd();">\n"""

IMG_MAP_SRC = """<img src="%s" border="0" ismap usemap="#points%s" width="%d" \
height="%d" />\n"""

MAP_SRC = """
<MAP name="points%s">
%s
</MAP>
"""

IMG_SRC = """<img src=\\'%s\\' border=1 />"""

IMG_SRC_2 = """<img src='%s' border=1 ismap usemap="#points%s" />&nbsp;<img src='%s' border=1/>"""
IMG_SRC_minus_legend = """<img src='%s' border=1 ismap usemap="#points%s" />"""
DOWNLOAD_LINK = """<a href=\\'%s\\' target="_blank">%s</a>"""

PDF_LINK = """<a href=\'%s\' target="_blank">%s</a>"""

LEGEND_LINK = """<a href=\'%s\' target="_blank">%s</a>"""


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


def make_legend(data_ids, colors, plot_width, plot_height, label_color,
                background_color, img_abs, generate_image_type, dpi):
    '''This is a simple function which generates a legend figure for plots '''

    img_path = os.path.splitext(img_abs)[0]

    out_fpath = img_path + '_legend.' + generate_image_type
    fname = os.path.split(out_fpath)[-1]

    # determine the number of text values and the max length of those values
    num_ids = len(data_ids)
    max_id_len = 0
    for i in data_ids:
        if len(i) > max_id_len:
            max_id_len = len(i)

    # numbers multiplied by were tweaked by hand
    figlegend = figure(figsize=(max_id_len * 0.15, num_ids * 0.22))

    # set some of the legend parameters
    fsize = 6
    rc('font', size=fsize)
    rc('text', color=label_color)
    rc('patch', linewidth=0)
    rc('axes', linewidth=0, edgecolor=label_color)
    rc('text', usetex=False)

    proxies = [mpatches.Patch(color=c, label=lab)
               for c, lab in zip(cycle(colors), data_ids)]
    l = figlegend.legend(handles=proxies, labels=data_ids,
                         loc='center left',
                         shadow=False, fancybox=False)
    l.legendPatch.set_alpha(0)

    figlegend.savefig(out_fpath, dpi=dpi, facecolor=background_color)
    close(figlegend)

    return fname


def make_pie_chart(data, dir_path, level, prefs, pref_colors, background_color,
                   label_color, generate_image_type,
                   plot_width, plot_height, bar_width, dpi, include_html_legend,
                   file_prefix=None, props={},
                   others_key="All Other Categories",
                   others_color="#eeeeee", should_capitalize=True):
    """
    Write interactive piechart
    data: [fraction:label,...]
    trunc_len: truncates labels after this many chars
    """

    # Raise value error if no data is supplied
    if not data:
        raise ValueError("No data available for pie chart.")

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
            capital = "%s (%.2f%%)" % (c_label.capitalize(), (c_frac * 100.0))
            all_labels.append(capital)
        else:
            all_labels.append("%s (%.2f%%)" % (c_label, (c_frac * 100.0)))

    # define figure parameters
    rc('font', size='10')
    rc('text', color=label_color)
    rc('patch', linewidth=.1)
    rc('axes', linewidth=.5, edgecolor=label_color)
    rc('text', usetex=False)

    # generate figure object
    fig = figure(randrange(10000), figsize=(plot_width, plot_height))

    fp = FontProperties()
    fp.set_size('8')
    if len(data) > 30:
        loc = 4
    else:
        loc = 5

    # define title and legend
    mtitle = "Pie Chart"
    if "title" in props:
        mtitle = props["title"]
    axis('off')

    # create pie object and legend
    title(mtitle, fontsize='10', color=label_color)
    ax = axes([0.0, 0.0, .5, 1])
    p1 = pie(all_fracs, shadow=False, colors=colors)

    # Remove legend from plotting area
    # flg = figlegend(p1[0],labels = all_labels, loc = loc, borderpad=0.3, \
    #             labelspacing=0.3, prop = fp)
    # flg.legendPatch.set_alpha(0.0)

    # write out
    if file_prefix is None:
        img_name = make_img_name(file_ext='.png')
    else:
        img_name = file_prefix

    # define filepath
    img_abs = os.path.join(dir_path, 'charts', img_name)
    savefig(img_abs, dpi=dpi, facecolor=background_color)
    eps_link = ""
    eps_abs = ""

    # generate the image as a pdf
    if file_prefix is None:
        eps_img_name = make_img_name(file_ext=".%s" % generate_image_type)
    else:
        eps_img_name = file_prefix + ".%s" % generate_image_type

    savefig(os.path.join(dir_path, 'charts', eps_img_name),
            facecolor=background_color)

    # generate the image as an eps
    if generate_image_type == 'eps':
        strip_eps_font(os.path.join(dir_path, 'charts', eps_img_name))

    eps_abs = os.path.join(dir_path, 'charts', eps_img_name)
    eps_link = PDF_LINK % (os.path.join('charts',
                                        eps_img_name), 'View Figure (.%s)' % generate_image_type)

    close(fig)
    clf()

    # this converts the sample_ids to a sample_id array and a corresponding
    # color array, so we can set the order based on array
    updated_taxa = []
    updated_colors = []
    for i in data:
        if i[0] != others_key:
            updated_taxa.append(i[0].replace('"', ''))
            updated_colors.append(data_colors[pref_colors[i[0]]].toHex())
        else:
            updated_taxa.append(others_key)
            updated_colors.append(others_color)

    if include_html_legend:
        # first make a low-res png legend for display on webpage
        legend_fname_png = make_legend(updated_taxa, updated_colors,
                                       plot_width, plot_height, label_color,
                                       background_color, img_abs, 'png', 80)
        legend_fpath_png = (os.path.join('charts', legend_fname_png))

    # make hi-res legend
    legend_fname = make_legend(updated_taxa, updated_colors,
                               plot_width, plot_height, label_color,
                               background_color, img_abs, generate_image_type, dpi)

    legend_fpath = (os.path.join('charts', legend_fname))
    legend_link = LEGEND_LINK % (legend_fpath, 'View Legend (.%s)' %
                                 generate_image_type)

    # the points id is generic along with xmap_html since they are used for
    # area and bar charts
    points_id = ''
    xmap_html = ''

    if not include_html_legend:
        IMG_TEXT = IMG_SRC_minus_legend % (os.path.join('charts', img_name),
                                           points_id)
    else:
        IMG_TEXT = IMG_SRC_2 % (os.path.join('charts', img_name),
                                points_id, legend_fpath_png)

    return eps_link, legend_link, IMG_TEXT, xmap_html


def transform_and_generate_xmap(
        ax1, bar_y_data, bar_width, taxa, x, plot_height,
        dpi, taxa_percents, sample_ids, chart_type):
    '''This function takes the bar graph data and generate html coordinates
       which can be used with an area map'''

    # transform the data into coordinates
    ax1.set_transform(ax1.axes.transData)
    trans = ax1.get_transform()
    # get transformed edges for the start/end of the plot
    starty = trans.transform((0, 0))[1]
    endy = trans.transform((0, 1))[1]
    startx = trans.transform((0, 0))[0]
    endx = trans.transform((len(x) - 1, 0))[0]

    # define the distance between x-axis points, so we can estimate width of
    # rectangle to use
    iterx = endx - startx

    # iterate over the bar graph data and put them into x/y coordinate arrays
    all_cids = []
    all_xcoords = []
    all_ycoords = []
    for s, t in enumerate(bar_y_data):
        y_start = t + numpy.sum(bar_y_data[:s], axis=0)
        icoords = trans.transform(zip(x, y_start))
        xcoords, ycoords = zip(*icoords)
        all_cids.append(taxa[s])
        all_xcoords.append(xcoords)
        all_ycoords.append(ycoords)

    # define xmap array and height of image, which is based on plot_height
    xmap = []
    img_height = plot_height * 80

    # determine width of rectangle to use
    if chart_type == 'area':
        half_iterx = iterx / 2 * bar_width
        # half_iterx=iterx*bar_width/17
    elif chart_type == 'bar':
        if len(all_xcoords[0]) > 1:
            half_iterx = (
                (all_xcoords[0][1] - all_xcoords[0][0]) / 2) * bar_width
        else:
            half_iterx = all_xcoords[0][0] / 2

    # iterate over y-coordinates and define area_map
    for i, j in enumerate(all_ycoords):
        # iterate over the x-coordinates
        for r, s in enumerate(x):
            # when in the first row, starting points are different, so need to
            # handle this case
            if i != 0:
                # we take the height of the image minus the distance up to the
                # current point to get the min of the field, then take the
                # height minus the current point to get the total height of
                # rectangle
                # use +/- half distance between x(0)->x(1) to make the width
                # of rectangle
                prev = i - 1

                # also we don't want to plot points that have an area of 0
                if all_ycoords[i][r] != all_ycoords[prev][r]:
                    xmap.append(AREA_SRC % (all_xcoords[i][r] - half_iterx,
                                            img_height - all_ycoords[prev][r],
                                            all_xcoords[i][r] + half_iterx,
                                            img_height - all_ycoords[i][r],
                                            taxa[i], sample_ids[r] + ':' + taxa[i], taxa_percents[r][i] * 100))
            else:
                # if at the beginning of the array and the value is not 0
                if all_ycoords[i][r] != starty:
                    xmap.append(AREA_SRC % (all_xcoords[i][r] - half_iterx,
                                            img_height -
                                            starty, all_xcoords[i][r] +
                                            half_iterx,
                                            img_height - all_ycoords[i][r],
                                            taxa[i], sample_ids[r] +
                                            ':' + taxa[i],
                                            taxa_percents[r][i] * 100))

    return xmap


def make_area_bar_chart(
        sample_ids, taxa_percents, taxa, dir_path, level, prefs,
        pref_colors,
        background_color, label_color, chart_type,
        generate_image_type,
        plot_width, plot_height, bar_width, dpi, resize_nth_label,
        label_type, include_html_legend,
        file_prefix=None, props={},
        others_key="All Other Categories",
        others_color="#eeeeee", should_capitalize=True):
    """
    Write interactive area chart
    data: [fraction:label,...]
    trunc_len: truncates labels after this many chars

    This function is very similar to the make_pie_charts function
    """
    # verify there is data in the file
    if not taxa_percents:
        raise ValueError("No data available for area chart.")

    all_fracs = []
    all_labels = []
    colors = []

    # set font-size based on the number of samples
    fsize = 0
    for i in range(7):
        fsize = 11 - i
        if len(sample_ids) <= (i * 10):
            break

    # define figure parameters
    rc('font', size=fsize)
    rc('text', color=label_color)
    rc('patch', linewidth=.1)
    rc('axes', linewidth=0, edgecolor=background_color)
    rc('text', usetex=False)
    rc('xtick', labelsize=fsize, color=label_color)

    # define figure
    fig = figure(figsize=(plot_width, plot_height))
    ax1 = fig.add_subplot(111, axisbg=background_color)

    # change the tick colors and width
    for tick in ax1.xaxis.get_ticklines():
        tick.set_color(label_color)
        tick.set_markersize(0)
        tick.set_markeredgewidth(.5)

    # create an iterative array for length of sample_ids
    if label_type == 'categorical':
        x = numpy.arange(0, len(sample_ids))
    elif label_type == 'numeric':
        x = sorted(map(lambda x: float(x), sample_ids))

        # numerical numbers must be sorted or else it gets screwed up in the
        # plot
    else:
        raise ValueError('Label type is not valid!')

    # get the raw data into a form, we can use for plotting areas and bars
    y_data = numpy.row_stack((zip(*taxa_percents)))
    bar_y_data = zip(*taxa_percents)
    y_data_stacked = numpy.cumsum(y_data, axis=0)

    # if area chart we use fill_between
    if chart_type == 'area':
        # bar_width is for mouseovers, and since area charts are more polygonal
        # we use a small width, so user can at least mouseover on the x-axis
        # positions
        bar_width = 0.005
        # fill the first taxa
        ax1.fill_between(x, 0, y_data_stacked[0, :], linewidth=0,
                         facecolor=data_colors[pref_colors[taxa[0]]].toHex(),
                         alpha=1)

        # fill all taxa up to the last one
        for i, j in enumerate(y_data_stacked):
            if i < len(y_data_stacked) - 1:
                next = i + 1
                ax1.fill_between(x, y_data_stacked[i, :],
                                 y_data_stacked[next, :], linewidth=0,
                                 facecolor=data_colors[
                                     pref_colors[taxa[i + 1]]].toHex(),
                                 alpha=1)
            # fill the last taxa to the total height of 1/
            else:
                ax1.fill_between(x, y_data_stacked[i, :], 1, linewidth=0,
                                 facecolor=data_colors[
                                     pref_colors[taxa[i]]].toHex(),
                                 alpha=1)

        # this cleans up the whitespace around the subplot
        # ax1.set_xlim((0,len(x)))
        ax1.set_ylim((0, 1))

    # if area chart we use bar
    elif chart_type == 'bar':

        # iterate over the data and make stacked bars
        for i, j in enumerate(bar_y_data):
            # if we are not in the first row of array, append more taxa
            if i > 0:
                ax1.bar(x, bar_y_data[i], width=bar_width, linewidth=0,
                        color=data_colors[pref_colors[taxa[i]]].toHex(),
                        bottom=numpy.sum(bar_y_data[:i], axis=0), align='center')
            # make the bars for the first row of array
            else:
                ax1.bar(x, bar_y_data[i], width=bar_width, linewidth=0,
                        color=data_colors[pref_colors[taxa[i]]].toHex(),
                        align='center')
        # this cleans up the whitespace around the subplot
        # ax1.set_xlim((-0.5,len(sample_ids)-0.5))
        ax1.set_ylim((0, 1))

    # transform bar_data into an area map for html mouseovers
    xmap = transform_and_generate_xmap(ax1, bar_y_data, bar_width, taxa, x,
                                       plot_height, dpi, taxa_percents, sample_ids,
                                       chart_type)

    # rename each area map based on the level passed in.
    points_id = 'rect%s' % (level)

    # append the area map html
    map_html = MAP_SRC % (points_id, ''.join(xmap))

    # set the values for the x-ticks
    # ax1.xaxis.set_ticks(x)
    # create an iterative array for length of sample_ids
    if label_type == 'numeric':
        x_axis_labels = numpy.arange(min(x), max(x))
        output_labels = []
        if resize_nth_label > 0:
            # xlabels=ax1.get_xticklabels()
            start = x_axis_labels[0]
            iterator_size = 0
            for i, l in enumerate(x_axis_labels):
                if l == start:
                    output_labels.append(l)
                else:
                    if (iterator_size == (resize_nth_label - 1)):
                        output_labels.append(l)
                        iterator_size = 0
                    else:
                        iterator_size = iterator_size + 1

            ax1.xaxis.set_ticks(output_labels)
            ax1.set_xticklabels(output_labels, rotation='vertical')
        else:
            ax1.xaxis.set_ticks(x)
            ax1.set_xticklabels(x, rotation='vertical')

    else:
        x_axis_labels = numpy.arange(0, len(sample_ids))
        # if the user would like to create larger labels for every nth label
        # this iterates over the labels and adds the value 4 to the font-size,
        # thereby making the fontsize larger
        x_tick_locations = []
        output_labels = []
        if resize_nth_label > 0:
            xlabels = ax1.get_xticklabels()
            iterator_size = 0
            for i, l in enumerate(x_axis_labels):
                if i == 0:
                    output_labels.append(sample_ids[i])
                    x_tick_locations.append(i)
                else:
                    if iterator_size == (resize_nth_label - 1):
                        output_labels.append(sample_ids[i])
                        x_tick_locations.append(i)
                        iterator_size = 0
                    else:
                        iterator_size = iterator_size + 1
            ax1.xaxis.set_ticks(x_tick_locations)
            ax1.set_xticklabels(output_labels, rotation='vertical')
        else:
            ax1.xaxis.set_ticks(x)
            ax1.set_xticklabels(sample_ids, rotation='vertical')

    ax1.set_yticks([])

    # write out
    if file_prefix is None:
        img_name = make_img_name(file_ext='.png')
    else:
        img_name = file_prefix

    # define filepath
    img_abs = os.path.join(dir_path, 'charts', img_name)
    savefig(img_abs, dpi=80, facecolor=background_color)
    eps_link = ""
    eps_abs = ""

    # generate the image as a pdf
    if file_prefix is None:
        eps_img_name = make_img_name(file_ext=".%s" % generate_image_type)
    else:
        eps_img_name = file_prefix + ".%s" % generate_image_type

    savefig(os.path.join(dir_path, 'charts', eps_img_name),
            facecolor=background_color)

    # generate the image as an eps
    if generate_image_type == 'eps':
        strip_eps_font(os.path.join(dir_path, 'charts', eps_img_name))

    eps_abs = os.path.join(dir_path, 'charts', eps_img_name)
    eps_link = PDF_LINK % (os.path.join('charts',
                                        eps_img_name), 'View Figure (.%s)' % generate_image_type)

    close(fig)
    clf()

    # this converts the sample_ids to a sample_id array and a corresponding
    # color array, so we can set the order based on array
    updated_taxa = []
    updated_colors = []
    for i in taxa:
        if i != others_key:
            updated_taxa.append(i.replace('"', ''))
            updated_colors.append(data_colors[pref_colors[i]].toHex())
        else:
            updated_taxa.append(others_key)
            updated_colors.append(others_color)

    if include_html_legend:
        # first make a low-res png legend for display on webpage
        legend_fname_png = make_legend(updated_taxa, updated_colors,
                                       plot_width, plot_height, label_color,
                                       background_color, img_abs, 'png', 80)

        legend_fpath_png = (os.path.join('charts', legend_fname_png))

    # make high-res legend
    legend_fname = make_legend(updated_taxa, updated_colors,
                               plot_width, plot_height, label_color,
                               background_color, img_abs, generate_image_type, dpi)

    legend_fpath = (os.path.join('charts', legend_fname))
    legend_link = LEGEND_LINK % (legend_fpath, 'View Legend (.%s)' %
                                (generate_image_type))

    if not include_html_legend:
        IMG_TEXT = IMG_SRC_minus_legend % (os.path.join('charts', img_name),
                                           points_id)
    else:
        IMG_TEXT = IMG_SRC_2 % (os.path.join('charts', img_name),
                                points_id, legend_fpath_png)

    return eps_link, legend_link, IMG_TEXT, map_html


def make_img_name(file_ext='.png'):
    """ Generate a random file name """
    fn = []
    # format seqs and write out to temp file
    for i in range(0, 30):
        fn.append(choice(ALPHABET))
    return ''.join(fn) + file_ext


def write_html_file(out_table, outpath):
    """Write pie charts into an html file"""
    page_out = PAGE_HTML % ('Taxa Summaries', out_table)
    out = open(outpath, "w+")
    out.write(page_out)
    out.close()


def get_fracs(counts, num_categories, total, chart_type, sort_data=True):
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

    area_table_out = []

    # this loop iterates over the OTU table and generates html code for the
    # data table
    for j, (n, t, s) in enumerate(counts):
        frac = float(n) / total
        if chart_type == 'pie':
            if j < num_categories - 1:
                red += n
                fracs_labels_other.append((t, frac))
        elif chart_type == 'area' or chart_type == 'bar':
            if j < num_categories:
                red += n
                fracs_labels_other.append((t, frac))

        tax = s.strip().split("<br>")[-1]
        tax = tax.replace('"', '')
        for_overlib = s.strip().rpartition("<br>")[0]
        for_overlib = for_overlib.replace('"', '')

        # Added this code because the data table is being presented
        # differently for the area charts
        if chart_type == 'pie':
            all_counts.append(
                DATA_HTML % (
                    frac * 100,
                    for_overlib,
                    tax,
                    tax,
                    t))
        elif chart_type == 'area' or chart_type == 'bar':
            area_table_out.append(str(n))

    # returning a dictionary for the case of area charts, which is different
    # than the array passed by the pie charts
    if chart_type == 'area' or chart_type == 'bar':
        all_counts = area_table_out

    if len(counts) > num_categories:
        other_cat = len(counts) - (num_categories - 1)
        new_counts = counts[0:num_categories - 1]
        other = sum([c_over[0] for c_over in counts[num_categories - 1:]])
        other_frac = float(other) / total
        fracs_labels = [(t, float(n) / red) for n, t, s in new_counts]

    # added in the ability to turn off sorting, since we want the data to be
    # unsorted for the area charts
    if sort_data:
        fracs_labels_other.sort()
        fracs_labels.sort()

    return (
        fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac
    )


def make_HTML_table(l, other_frac, total, red, other_cat, fracs_labels_other,
                    fracs_labels, dir_path, all_counts, level,
                    prefs, pref_colors, background_color, label_color, chart_type,
                    label, generate_image_type,
                    plot_width, plot_height, bar_width, dpi, resize_nth_label,
                    label_type, include_html_legend):
    """Makes the HTML table for one set of charts """
    img_data = []

    # generate html for pie charts
    if chart_type == 'pie':
        # in the case the user wants to trim down the number of taxa
        if other_cat > 0:
            # first generate the pie charts containing an other group for all
            # taxa below the cutoff.
            fracs_labels_other.append(("All Other Categories", other_frac))
            title = TITLE_include % (l, total, total,
                                     len(fracs_labels_other), total - red, other_cat)
            all_taxons = [l]
            pie_charts_placement = []

            # make pie chart image
            pie = make_pie_chart(fracs_labels_other, dir_path, level,
                                 prefs, pref_colors, background_color, label_color,
                                 generate_image_type,
                                 plot_width, plot_height, bar_width, dpi,
                                 include_html_legend,
                                 props={'title': title})

            pie_charts_placement.append(pie[0] + '&nbsp;&nbsp;' + pie[1] +
                                        '</td></tr><tr><td>' + pie[2] +
                                        '</td></tr><tr><td class="ntitle">')

            # second generate the pie charts where the other category is removed
            # and percents are recalculated
            title = TITLE_exclude % (l, red, total, len(fracs_labels),
                                     total - red, other_cat)

            # make pie chart image
            pie = make_pie_chart(fracs_labels, dir_path, level,
                                 prefs, pref_colors, background_color, label_color,
                                 generate_image_type,
                                 plot_width, plot_height, bar_width, dpi,
                                 include_html_legend,
                                 props={'title': title})

            pie_charts_placement.append(pie[0] + '&nbsp;&nbsp;' + pie[1] +
                                        '</td></tr><tr><td class="ntitle">' +
                                        pie[2])

            all_taxons.extend(pie_charts_placement)
            all_taxons.extend((" ", " "))

            # put the charts into the html image data
            img_data.append(TABLE_graph % tuple(all_taxons))
            img_data.append(DATA_TABLE_HTML % ''.join(all_counts))

        else:
            # if there is no category cutoff generate plots, without other cat
            title = TITLE % (l, total, total, len(fracs_labels_other))
            all_taxons = [l]

            # make pie chart image
            pie = make_pie_chart(fracs_labels_other, dir_path, level,
                                 prefs, pref_colors, background_color, label_color,
                                 generate_image_type,
                                 plot_width, plot_height, bar_width, dpi,
                                 include_html_legend,
                                 props={'title': title})

            all_taxons.extend(pie)

            # put the charts into the html image data
            img_data.append(TABLE_graph % tuple(all_taxons))
            img_data.append(DATA_TABLE_HTML % ''.join(all_counts))

    # generate html for bar and area charts
    elif chart_type == 'area' or chart_type == 'bar':

        taxa_percents = fracs_labels_other
        sample_ids = l
        taxa = other_cat

        all_categories = []
        title = TITLE % (label, total, total, len(fracs_labels_other))
        all_taxons = [label]

        # make area chart image
        area = make_area_bar_chart(sample_ids, taxa_percents, taxa, dir_path,
                                   level, prefs, pref_colors,
                                   background_color, label_color, chart_type,
                                   generate_image_type,
                                   plot_width, plot_height, bar_width, dpi,
                                   resize_nth_label, label_type,
                                   include_html_legend,
                                   props={'title': title})

        all_taxons.extend(area)

        # put the charts into the html image data
        img_data.append(TABLE_graph % tuple(all_taxons))

    return img_data


def get_counts(label, colorby, num_categories, dir_path, level, color_data,
               prefs, pref_colors, background_color, label_color, chart_type,
               generate_image_type, plot_width, plot_height,
               bar_width, dpi, raw_fpath, resize_nth_label, label_type,
               include_html_legend):
    """gets all the counts for one input file"""

    img_data = []
    labels = []
    level_counts = []

    sample_ids, otu_ids, otu_table = color_data
    labels = sample_ids

    # iterate over the counts table and cleanup taxa labels
    for idx, counts in enumerate(otu_table):
        taxonomy = otu_ids[idx]
        split_label = [i for i in taxonomy.strip().split(";")]
        taxonomy = ';'.join(split_label)
        level_counts.append((sum(map(float, counts)), taxonomy,
                            '<br>'.join(split_label)))
    all_sum = sum([c_over[0] for c_over in level_counts])

    # get the fractions for all samples
    fracs_labels_other, fracs_labels, all_counts, other_cat, red, other_frac = \
        get_fracs(level_counts, num_categories, all_sum,
                  chart_type, True)

    # if making pie charts we perform a couple extra steps, such as making a
    # total pie chart
    if chart_type == 'pie':
        # make the total pie chart
        img_data.extend(
            make_HTML_table(label, other_frac, all_sum, red, other_cat,
                            fracs_labels_other, fracs_labels, dir_path,
                            all_counts, level, prefs, pref_colors,
                            background_color, label_color, chart_type,
                            label, generate_image_type, plot_width,
                            plot_height, bar_width, dpi, resize_nth_label,
                            label_type, include_html_legend))

        if colorby is not None:
            # in the case the user specifies only certain samples we need to
            # handle that case
            for i, l in enumerate(sample_ids):
                if l not in colorby:
                    continue
                total = 0
                sample_counts = []
                for idx, counts in enumerate(otu_table):
                    taxonomy = otu_ids[idx]
                    split_label = [j for j in taxonomy.strip().split(";")]
                    taxonomy = ';'.join(split_label)
                    c = float(counts[i])
                    if c > 0:
                        total += c
                        sample_counts.append((c, taxonomy,
                                              '<br>'.join(split_label)))

                # get fractions for specific samples
                fracs_labels_other, fracs_labels, all_counts,\
                    other_cat, red, other_frac = get_fracs(sample_counts,
                                                           num_categories,
                                                           total, chart_type, True)

                # make the per sample pie charts
                img_data.extend(make_HTML_table('_'.join([label, l.strip()]),
                                                other_frac, total, red, other_cat, fracs_labels_other,
                                                fracs_labels, dir_path, all_counts, level,
                                                prefs, pref_colors, background_color, label_color,
                                                chart_type, l.strip(
                                                ), generate_image_type,
                                                plot_width, plot_height, bar_width, dpi,
                                                resize_nth_label, label_type, include_html_legend))

    # if making an area/bar chart we do not make per sample images, instead
    # we make a total chart only
    elif chart_type == 'area' or chart_type == 'bar':
        area_plot_arr = []
        area_plot_sample_ids = []
        area_plot_taxa_arr = []
        taxa_html = []
        total_area_table_out = []
        total_sums = []

        if colorby is not None:
            # in the case the user specifies only certain samples we need to
            # handle that case
            for i, l in enumerate(sample_ids):
                if l not in colorby:
                    continue
                total = 0
                area_plot_sample_ids.append(l)
                sample_counts = []

                # iterate over the counts and cleanup taxa for this particular
                # fxn

                #add in zip
                for idx, counts in enumerate(otu_table):
                    taxonomy = otu_ids[idx]
                    split_label = [j for j in taxonomy.strip().split(";")]
                    taxonomy = ';'.join(split_label)

                    c = float(counts[i])
                    total += c
                    sample_counts.append(
                        (c, taxonomy, '<br>'.join(split_label)))

                # get fractions for specific samples
                fracs_labels_other, fracs_labels, all_counts,\
                    other_cat, red, other_frac = get_fracs(sample_counts,
                                                           len(sample_counts), total,
                                                           chart_type, False)

                total_area_table_out.append(all_counts)
                total_sums.append(sum([float(i) for i in all_counts]))

                # get the percents for each taxa and sample
                area_plot_per = []
                area_plot_taxa = []
                for i in fracs_labels_other:
                    area_plot_per.append(i[1])
                    area_plot_taxa.append(i[0])

                area_plot_arr.append(area_plot_per)
                area_plot_taxa_arr.append(area_plot_taxa)

        # write out the data table html, since it is different than pie chart
        # data table
        taxa_html.append('<tr><th>' + l.strip() +
                         '</th></tr>' + ''.join(all_counts) + '')

        data_table = zip(*total_area_table_out)

        # create link for raw data file
        data_html_str = '<table><tr class=ntitle><td><a href="%s" target="_blank">View Table (%s)</a></td></tr></table>' % \
            (os.path.join('raw_data', os.path.split(raw_fpath)[-1]),
             os.path.splitext(raw_fpath)[-1])

        # create the output table
        data_html_str += '<table cellpadding=1 cellspacing=1 border=1 ' + \
            'style=\"text-align:center;border-color:white;' +\
            'border-style:groove;\">' + \
            '<tr class=\"ntitle\"><td class=\"header\" colspan="2"></td><td' +\
            ' valign=\"bottom\" class=\"header\">Total</td>'

        ct_head_row = '<tr class=ntitle>' + \
            '<td valign=\"bottom\" ' + \
            'class=\"header\">Legend</td><td ' + \
            'valign=\"bottom\" class=\"header\">Taxonomy</td>' + \
            '<td class=\"header\">%</td>'

        for i in area_plot_sample_ids:
            data_html_str += '<td valign=bottom class=header>%s</td>' % (i)
            ct_head_row += '<td class=\"header\">%</td>'

        data_html_str += '</tr>'
        ct_head_row += '</tr>'
        data_html_str += ct_head_row
        table_sum = sum(total_sums)

        # list taxa in first row
        for ct, dat in enumerate(otu_ids):
            tax = dat
            split_label = [i for i in tax.strip().split(";")]
            split_label[-1] = "<a href=javascript:gg(\'%s\');>%s</a>" % \
                (split_label[-1].replace(' ', '+'),
                 split_label[-1].replace(' ', '&nbsp;'))
            joined_label = ';'.join(split_label).replace('"', '')
            row_sum = sum([float(i) for i in data_table[ct]])
            data_html_str += "<tr><td class=\"normal\" bgcolor=\"%s\">&nbsp;&nbsp;</td><td style=\"text-align:left;\" class=\"normal\">%s</td><td class=\"normal\">%5.1f&#37;</td>"\
                % (data_colors[pref_colors[tax]].toHex(), joined_label,
                   row_sum / table_sum * 100)

            # add the percent taxa for each sample
            for i, per_tax in enumerate(data_table[ct]):
                if float(per_tax) > 0:
                    data_html_str += '<td class=\"normal\" style=\"border-color:%s;\">%5.1f&#37;</td>' %\
                        (data_colors[pref_colors[tax]].toHex(),
                         (float(per_tax) / total_sums[i] * 100))
                else:
                    data_html_str += '<td class=\"normal\">%5.1f&#37;</td>' % \
                        (float(per_tax) / total_sums[i] * 100)

            data_html_str += '</tr>'

        data_html_str += '</table>'

        # make sure that the taxa array is in the proper order
        for i in range(len(area_plot_taxa_arr) - 1):
            if area_plot_taxa_arr[i] != area_plot_taxa_arr[i + 1]:
                raise ValueError('The taxonomies are out of order!')

        # add data to the html output
        img_data.extend(make_HTML_table(area_plot_sample_ids,
                                        other_frac, all_sum, red, otu_ids, area_plot_arr,
                                        fracs_labels, dir_path, [' '.join(
                                            taxa_html)], level,
                                        prefs, pref_colors, background_color, label_color, chart_type,
                                        label, generate_image_type,
                                        plot_width, plot_height, bar_width, dpi, resize_nth_label,
                                        label_type, include_html_legend))
        img_data.append(data_html_str)

    return img_data


def make_all_charts(data, dir_path, filename, num_categories, colorby, args,
                    color_data, prefs, background_color, label_color,
                    chart_type, generate_image_type, plot_width, plot_height,
                    bar_width, dpi, resize_nth_label, label_type,
                    include_html_legend):
    """Generate interactive charts in one HTML file"""

    # iterate over the preferences and assign colors according to taxonomy
    img_data = []
    for label, f_name in data:
        raw_fpath = os.path.join(
            dir_path,
            'raw_data',
            os.path.split(f_name)[-1])
        # move raw file to output directory
        shutil.copyfile(f_name, raw_fpath)

        f = color_data['counts'][f_name]
        level = max([len(t.split(';')) - 1 for t in f[1]])

        for key in prefs.keys():
            if prefs[key]['column'] != str(level):
                continue
            col_name = 'Taxon'
            mapping = [['Taxon']]
            mapping.extend([[m] for m in f[1]])
            if 'colors' in prefs[key]:
                if isinstance(prefs[key]['colors'], dict):
                    pref_colors = prefs[key]['colors'].copy()
                    # copy so we can mutate
                else:
                    pref_colors = prefs[key]['colors'][:]
            else:
                pref_colors = {}
            labelname = prefs[key]['column']

            # Define groups and associate appropriate colors to each group
            groups = group_by_field(mapping, col_name)
            pref_colors, data_colors, data_color_order = \
                get_group_colors(groups, pref_colors)

        updated_pref_colors = {}

        if chart_type == 'area' and len(f[0]) == 1:
            raise ValueError(
                'When generating area charts, the number of samples (or category values) must be greater than 1.  However, you can still produce a pie chart or bar chart with only 1 sample (or category value), but you must remove the area chart value from the input arguments.')

        for key in pref_colors:
            updated_pref_colors[key.replace('"', '')] = pref_colors[key]

        for i, val in enumerate(f[1]):
            f[1][i] = val.replace('"', '')

        # parse the counts and continue processing
        img_data.extend(get_counts(label.strip(), colorby, num_categories,
                        dir_path, level, f, prefs, updated_pref_colors,
                        background_color,
                        label_color, chart_type, generate_image_type,
                        plot_width, plot_height, bar_width, dpi, raw_fpath,
                        resize_nth_label, label_type, include_html_legend))

    # generate html filepath
    outpath = os.path.join(dir_path, '%s_charts.html' % chart_type)
    out_table = ''.join(img_data)
    # write out html file
    write_html_file(out_table, outpath)
