#!/usr/bin/env python
# File created on 09 Feb 2010
# file make_2d_plots.py

__author__ = "Jesse Stombaugh and Micah Hamady"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself
__credits__ = ["Jesse Stombaugh", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

from matplotlib import use
use('Agg', warn=False)
from matplotlib.pyplot import rc
import matplotlib.pyplot as plt
from matplotlib.cbook import iterable
from matplotlib.patches import Ellipse
from matplotlib.font_manager import FontProperties
from numpy import asarray
from qiime.util import summarize_pcoas, isarray, qiime_system_call
from qiime.parse import parse_coords
from qiime.colors import iter_color_groups
import numpy as np
from tempfile import mkdtemp
import os

from qiime.util import MissingFileError, load_pcoa_files

SCREE_TABLE_HTML = """<table cellpadding=0 cellspacing=0 border=1>
<tr><th align=center colspan=3 border=0>Scree plot</th></tr>
<tr>
<td class=normal align=center border=0>%s</td>
</tr>
</table>
<br><br>"""

TABLE_HTML = """<table cellpadding=0 cellspacing=0 border=1>
<tr><th align=center colspan=3 border=0>%s</th></tr>
<tr>
<td class=normal align=center border=0>%s</td>
<td class=normal align=center border=0>%s</td>
<td class=normal align=center border=0>%s</td>
</tr>
</table>
<br><br>"""

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
    's',  # : square
    'o',  # : circle
    '^',  # : triangle up
    '>',  # : triangle right
    'v',  # : triangle down
    '<',  # : triangle left
    'd',  # : diamond
    'p',  # : pentagon
    'h',  # : hexagon
]
'''
data_colors={'blue':'#0000FF','lime':'#00FF00','red':'#FF0000', \
             'aqua':'#00FFFF','fuchsia':'#FF00FF','yellow':'#FFFF00', \
             'green':'#008000','maroon':'#800000','teal':'#008080', \
             'purple':'#800080','olive':'#808000', \
             'silver':'#C0C0C0','gray':'#808080'}
'''
default_colors = ['blue', 'lime', 'red', 'aqua', 'fuchsia', 'yellow', 'green',
                  'maroon', 'teal', 'purple', 'olive', 'silver', 'gray']


# This function used to live in make_3d_plots.py but in the Boulder sk-bio
# code sprint it got moved here to remove the 3D files.
def get_coord(coord_fname, method="IQR"):
    """Opens and returns coords location matrix and metadata.
       Also two spread matrices (+/-) if passed a dir of coord files.
       If only a single coord file, spread matrices are returned as None.
    """
    if not os.path.isdir(coord_fname):
        try:
            coord_f = open(coord_fname, 'U')
        except (TypeError, IOError):
            raise MissingFileError('Coord file required for this analysis')
        coord_header, coords, eigvals, pct_var = parse_coords(coord_f)
        return [coord_header, coords, eigvals, pct_var, None, None]
    else:
        master_pcoa, support_pcoas = load_pcoa_files(coord_fname)

        # get Summary statistics
        coords, coords_low, coords_high, eigval_average, coord_header = \
            summarize_pcoas(master_pcoa, support_pcoas, method=method)
        pct_var = master_pcoa[3]  # should be getting this from an average

        # make_3d_plots expects coord_header to be a python list
        coord_header = list(master_pcoa[0])
        return (
            [coord_header,
             coords,
             eigval_average,
             pct_var,
             coords_low,
             coords_high]
        )


def make_line_plot(
        dir_path, data_file_link, background_color, label_color, xy_coords,
        props, x_len=8, y_len=4, draw_axes=False, generate_eps=True):
    """ Write a line plot

    xy_coords: a dict of form
       {series_label:([x data], [y data], point_marker, color)}

    (code adapted from Micah Hamady's code)
    """
    rc('font', size='8')
    rc('axes', linewidth=.5, edgecolor=label_color)
    rc('axes', labelsize=8)
    rc('xtick', labelsize=8)
    rc('ytick', labelsize=8)
    fig, ax = plt.subplots(figsize=(x_len, y_len))
    mtitle = props.get("title", "Groups")
    x_label = props.get("xlabel", "X")
    y_label = props.get("ylabel", "Y")

    ax.set_title('%s' % mtitle, fontsize='10', color=label_color)
    ax.set_xlabel(x_label, fontsize='8', color=label_color)
    ax.set_ylabel(y_label, fontsize='8', color=label_color)

    sorted_keys = sorted(xy_coords.keys())

    for s_label in sorted_keys:
        s_data = xy_coords[s_label]
        c = s_data[3]
        m = s_data[2]
        ax.plot(s_data[0], s_data[1], c=c, marker=m, label=s_label,
                linewidth=.1, ms=5, alpha=1.0)

    fp = FontProperties()
    fp.set_size('8')
    ax.legend(prop=fp, loc=0)

    img_name = 'scree_plot.png'
    fig.savefig(
        os.path.join(dir_path,
                     img_name),
        dpi=80,
        facecolor=background_color)

    # Create zipped eps files
    eps_link = ""
    if generate_eps:
        eps_img_name = str('scree_plot.eps')
        fig.savefig(os.path.join(dir_path, eps_img_name), format='eps')
        out, err, retcode = qiime_system_call(
            "gzip -f " + os.path.join(dir_path, eps_img_name))
        eps_link = DOWNLOAD_LINK % ((os.path.join(data_file_link,
                                                  eps_img_name) +
                                     ".gz"), "Download Figure")

    return os.path.join(data_file_link, img_name), eps_link


def draw_scree_graph(dir_path, data_file_link, background_color, label_color,
                     generate_eps, data):
    """Draw scree plot

    (code adapted from Micah Hamady's code)
    """

    dimensions = len(data['coord'][3])

    props = {}
    props["title"] = "PCoA Scree Plot (First %s dimensions)" % dimensions
    props["ylabel"] = "Fraction of Variance"
    props["xlabel"] = "Principal component"

    xy_coords = {}
    x_points = [x for x in range(dimensions)]
    c_data = [float(x) / 100.0 for x in data['coord'][3]]
    xy_coords['Variance'] = (x_points, c_data, 'o', 'r')

    cum_var = [c_data[0]]
    for ix in range(dimensions - 1):
        cum_var.append(cum_var[ix] + c_data[ix + 1])
    xy_coords['Cumulative variance'] = (x_points, cum_var, 's', 'b')

    img_src, eps_link = make_line_plot(
        dir_path, data_file_link, background_color, label_color,
        xy_coords=xy_coords, props=props, x_len=4.5,
        y_len=4.5, generate_eps=generate_eps)

    return IMG_SRC % img_src, eps_link


def make_interactive_scatter(plot_label, dir_path, data_file_link,
                             background_color, label_color, sample_location,
                             alpha, xy_coords,
                             props, x_len=8, y_len=4, size=10,
                             draw_axes=False, generate_eps=True):
    """Write interactive plot

    xy_coords: a dict of form {series_label:([x data], [y data], \
    [xy point label],[color])}
    """
    my_axis = None
    rc('font', size='8')
    rc('patch', linewidth=0)
    rc('axes', linewidth=.5, edgecolor=label_color)
    rc('axes', labelsize=8)
    rc('xtick', labelsize=8, color=label_color)
    rc('ytick', labelsize=8, color=label_color)

    sc_plot = draw_scatterplot(props, xy_coords, x_len, y_len, size,
                               background_color, label_color, sample_location,
                               alpha)

    mtitle = props.get("title", "Groups")
    x_label = props.get("xlabel", "X")
    y_label = props.get("ylabel", "Y")
    ax = plt.gca()
    fig = ax.figure
    ax.set_title('%s' % mtitle, fontsize='10', color=label_color)
    ax.set_xlabel(x_label, fontsize='8', color=label_color)
    ax.set_ylabel(y_label, fontsize='8', color=label_color)

    if draw_axes:
        ax.axvline(linewidth=.5, x=0, color=label_color)
        ax.axhline(linewidth=.5, y=0, color=label_color)

    if my_axis is not None:
        ax.axis(my_axis)
    img_name = x_label[0:3] + '_vs_' + y_label[0:3] + '_plot.png'
    fig.savefig(os.path.join(dir_path, img_name),
                dpi=80, facecolor=background_color)

    # Create zipped eps files
    eps_link = ""
    if generate_eps:
        eps_img_name = str(x_label[0:3] + 'vs' + y_label[0:3] + 'plot.eps')
        fig.savefig(os.path.join(dir_path, eps_img_name), format='eps')
        out, err, retcode = qiime_system_call(
            "gzip -f " + os.path.join(dir_path, eps_img_name))
        eps_link = DOWNLOAD_LINK % ((os.path.join(data_file_link, eps_img_name)
                                     + ".gz"), "Download Figure")

    all_cids, all_xcoords, all_ycoords = transform_xy_coords(
        xy_coords, sc_plot)

    xmap, img_height, img_width = generate_xmap(
        x_len, y_len, all_cids, all_xcoords,
        all_ycoords)

    points_id = plot_label + x_label[2:3] + y_label[2:3]

    return IMG_MAP_SRC % (os.path.join(data_file_link, img_name), points_id,
                          img_width, img_height), MAP_SRC % \
        (points_id, ''.join(xmap)), eps_link


def generate_xmap(x_len, y_len, all_cids, all_xcoords, all_ycoords):
    """Generates the html interactive image map"""
    # Determine figure height and width"""
    img_height = x_len * 80
    img_width = y_len * 80

    # Write html script which allows for mouseover of labels
    xmap = []
    for cid, x, y in zip(all_cids, all_xcoords, all_ycoords):
        xmap.append(AREA_SRC % (x, img_height - y, cid, cid))

    return xmap, img_height, img_width


def draw_scatterplot(props, xy_coords, x_len, y_len, size, background_color,
                     label_color, sample_location, alpha):
    """Create scatterplot figure"""

    fig = plt.figure(figsize=(x_len, y_len))
    xPC = int(props['xlabel'][2:3])
    yPC = int(props['ylabel'][2:3])
    sorted_keys = xy_coords.keys()
    scatters = {}
    size_ct = shape_ct = 0

    xPC = xPC - 1
    yPC = yPC - 1
    # Iterate through coords and add points to the scatterplot
    for s_label in sorted_keys:
        s_data = xy_coords[s_label]
        if s_data[0] == []:
            pass
        else:
            c = s_data[3]
            m = s_data[4][0]

            ax = fig.add_subplot(111, axisbg=background_color)
            # set tick colors and width
            for line in ax.yaxis.get_ticklines():
                # line is a matplotlib.lines.Line2D instance
                line.set_color(label_color)
                line.set_markeredgewidth(1)

            for line in ax.xaxis.get_ticklines():
                # line is a matplotlib.lines.Line2D instance
                line.set_color(label_color)
                line.set_markeredgewidth(1)

            if isarray(s_data[5][0]) and isarray(s_data[6][0]) and \
                    isarray(s_data[7][0]):
                matrix_low = s_data[5][0]
                matrix_high = s_data[6][0]
                ellipse_ave = s_data[7][0]

                ellipse_x = [ellipse_ave[sample_location[s_label], xPC]]
                ellipse_y = [ellipse_ave[sample_location[s_label], yPC]]

                width = [np.fabs(matrix_high[sample_location[s_label], xPC] -
                                 matrix_low[sample_location[s_label], xPC])]
                height = [np.fabs(matrix_high[sample_location[s_label], yPC] -
                                  matrix_low[sample_location[s_label], yPC])]

                sc_plot = scatter_ellipse(ax, ellipse_x,
                                          ellipse_y, width, height, c=c, a=0.0,
                                          alpha=alpha)
                sc_plot.scatter(ellipse_x, ellipse_y, c=c, marker=m,
                                alpha=1.0)
            else:
                sc_plot = ax.scatter(s_data[0], s_data[1], c=c, marker=m,
                                     alpha=1.0, s=size, linewidth=1,
                                     edgecolor=c)
            size_ct += 1
            shape_ct += 1
            scatters[s_label] = sc_plot

    return sc_plot


def transform_xy_coords(xy_coords, sc_plot):
    """Transform the coords from the scatterplot into coords that can be \
       referenced in the html page"""
    sorted_keys = xy_coords.keys()
    all_cids = []
    all_xcoords = []
    all_ycoords = []
    sc_plot.set_transform(sc_plot.axes.transData)
    trans = sc_plot.get_transform()

    for s_label in sorted_keys:
        s_data = xy_coords[s_label]
        if s_data[0] == []:
            pass
        else:
            icoords = trans.transform(zip(s_data[0], s_data[1]))
            xcoords, ycoords = zip(*icoords)
            all_cids.extend(s_data[2])
            all_xcoords.extend(xcoords)
            all_ycoords.extend(ycoords)

    return all_cids, all_xcoords, all_ycoords


def draw_pcoa_graph(plot_label, dir_path, data_file_link, coord_1, coord_2,
                    coord_1r, coord_2r, mat_ave, sample_location,
                    data, prefs, groups, colors, background_color, label_color,
                    data_colors, data_color_order,
                    generate_eps=True,
                    pct_variation_below_one=False):
    """Draw PCoA graphs"""

    coords, pct_var = convert_coord_data_to_dict(data)

    if coord_1 not in coords:
        raise ValueError("Principal coordinate: %s not available." % coord_1)

    if coord_2 not in coords:
        raise ValueError("Principal coordinate: %s not available." % coord_2)

    # Handle matplotlib scale bug when all coords are 0.0
    if not len([x for x in map(float, coords[coord_2]) if x != 0.0]):
        for ix in range(len(coords[coord_2])):
            coords[coord_2][ix] = '1e-255'
    if not len([x for x in map(float, coords[coord_1]) if x != 0.0]):
        for ix in range(len(coords[coord_1])):
            coords[coord_1][ix] = '1e-255'

    # Write figure labels
    pct_exp1 = float(pct_var[coord_1])
    pct_exp2 = float(pct_var[coord_2])
    if float(pct_var['1']) < 1 and not pct_variation_below_one:
        pct_exp1 *= 100
        pct_exp2 *= 100
    props = {}
    props["title"] = "PCoA - PC%s vs PC%s" % (coord_1, coord_2)
    props["ylabel"] = "PC%s - Percent variation explained %.2f%%" \
        % (coord_2, pct_exp2)
    props["xlabel"] = "PC%s - Percent variation explained %.2f%%" \
        % (coord_1, pct_exp1)

    labels = coords['pc vector number']
    p1 = map(float, coords[coord_2])
    p2 = map(float, coords[coord_1])
    if isarray(coord_1r) and isarray(coord_2r) and isarray(mat_ave):
        p1r = coord_2r
        p2r = coord_1r
    else:
        p1r = None
        p2r = None
        mat_ave = None

    if len(p1) != len(p2):
        raise ValueError("Principal coordinate vectors unequal length.")
    p1d = dict(zip(labels, p1))
    p2d = dict(zip(labels, p2))

    alpha = data['alpha']

    xy_coords = extract_and_color_xy_coords(
        p1d, p2d, p1r, p2r, mat_ave, colors,
        data_colors, groups, coords)

    img_src, img_map, eps_link = make_interactive_scatter(
        plot_label, dir_path,
        data_file_link, background_color, label_color,
        sample_location, alpha,
        xy_coords=xy_coords, props=props, x_len=4.5,
        y_len=4.5, size=20, draw_axes=True,
        generate_eps=generate_eps)

    return img_src + img_map, eps_link


def extract_and_color_xy_coords(
        p1d, p2d, p1dr, p2dr, mat_ave, colors, data_colors,
        groups, coords):
    """Extract coords from appropriate columns and attach their \
       corresponding colors based on the group"""

    xy_coords = {}
    shape_ct = 0
    for group_name, ids in (groups.items()):

        color = data_colors[colors[group_name]].toHex()
        m = shape[shape_ct % len(shape)]
        shape_ct += 1
        for id_ in (ids):
            cur_labs = []
            cur_x = []
            cur_y = []
            cur_color = []
            cur_shape = []
            cur_1r = []
            cur_2r = []
            new_mat_ave = []
            if id_ in coords['pc vector number']:
                cur_labs.append(id_ + ': ' + group_name)
                cur_x.append(p2d[id_])
                cur_y.append(p1d[id_])
                cur_color.append(color)
                cur_shape.append(m)

                if isarray(p2dr) and isarray(p1dr) and isarray(mat_ave):
                    cur_1r.append(p1dr)
                    cur_2r.append(p2dr)
                    new_mat_ave.append(mat_ave)
                else:
                    cur_1r = [None]
                    cur_2r = [None]
                    new_mat_ave = [None]

            xy_coords["%s" % id_] = (cur_x, cur_y, cur_labs, cur_color,
                                     cur_shape, cur_1r, cur_2r, new_mat_ave)

    return xy_coords


def create_html_filename(coord_filename, name_ending):
    """Generate html filename using the given coord filename"""
    outpath = coord_filename.split('/')[-1] + name_ending
    return outpath


def convert_coord_data_to_dict(data):
    """Convert the coord data into a dictionary"""
    coord_header = data['coord'][0]
    coords = data['coord'][1]
    pct_var = data['coord'][3]
    coords_dict = {}
    pct_var_dict = {}
    coords_dict['pc vector number'] = coord_header
    for x in range(len(coords)):
        coords_dict[str(x + 1)] = coords[0:, x]
        pct_var_dict[str(x + 1)] = pct_var[x]

    return coords_dict, pct_var_dict


def write_html_file(out_table, outpath):
    """Write 2D plots into an html file"""
    page_out = PAGE_HTML % (outpath, out_table)
    out = open(outpath, "w+")
    out.write(page_out)
    out.close()


def generate_2d_plots(prefs, data, html_dir_path, data_dir_path, filename,
                      background_color, label_color, generate_scree,
                      pct_variation_below_one):
    """Generate interactive 2D scatterplots"""
    coord_tups = [("1", "2"), ("3", "2"), ("1", "3")]
    mapping = data['map']
    out_table = ''
    # Iterate through prefs and generate html files for each colorby option
    # Sort by the column name first
    sample_location = {}

    groups_and_colors = iter_color_groups(mapping, prefs)
    groups_and_colors = list(groups_and_colors)

    for i in range(len(groups_and_colors)):
        labelname = groups_and_colors[i][0]
        groups = groups_and_colors[i][1]
        colors = groups_and_colors[i][2]
        data_colors = groups_and_colors[i][3]
        data_color_order = groups_and_colors[i][4]

        data_file_dir_path = mkdtemp(dir=data_dir_path)

        new_link = os.path.split(data_file_dir_path)
        data_file_link = os.path.join('.', os.path.split(new_link[-2])[-1],
                                      new_link[-1])

        img_data = {}
        plot_label = labelname

        if 'support_pcoas' in data:
            matrix_average, matrix_low, matrix_high, eigval_average, m_names = \
                summarize_pcoas(data['coord'], data['support_pcoas'],
                                method=data['ellipsoid_method'])
            data['coord'] = \
                (m_names, matrix_average, data['coord'][2], data['coord'][3])
            for i in range(len(m_names)):
                sample_location[m_names[i]] = i
        else:
            matrix_average = None
            matrix_low = None
            matrix_high = None

            m_names = None

        for coord_tup in coord_tups:
            if isarray(matrix_low) and isarray(matrix_high) and \
                    isarray(matrix_average):
                coord_1r = asarray(matrix_low)
                coord_2r = asarray(matrix_high)
                mat_ave = asarray(matrix_average)
            else:
                coord_1r = None
                coord_2r = None
                mat_ave = None
                sample_location = None

            coord_1, coord_2 = coord_tup
            img_data[coord_tup] = draw_pcoa_graph(
                plot_label, data_file_dir_path,
                data_file_link, coord_1, coord_2,
                coord_1r, coord_2r, mat_ave,
                sample_location,
                data, prefs, groups, colors,
                background_color, label_color,
                data_colors, data_color_order,
                generate_eps=True,
                pct_variation_below_one=pct_variation_below_one)

        out_table += TABLE_HTML % (labelname,
                                   "<br>".join(img_data[("1", "2")]),
                                   "<br>".join(img_data[("3", "2")]),
                                   "<br>".join(img_data[("1", "3")]))

    if generate_scree:
        data_file_dir_path = mkdtemp(dir=data_dir_path)
        new_link = os.path.split(data_file_dir_path)
        data_file_link = os.path.join(
            '.',
            os.path.split(new_link[-2])[-1],
            new_link[-1])

        img_src, download_link = draw_scree_graph(
            data_file_dir_path, data_file_link, background_color,
            label_color, generate_eps=True, data=data)

        out_table += SCREE_TABLE_HTML % ("<br>".join((img_src, download_link)))

    outfile = create_html_filename(filename, '.html')
    outfile = os.path.join(html_dir_path, outfile)

    write_html_file(out_table, outfile)


def scatter_ellipse(axis_ob, x, y, w, h, c='b', a=0.0, alpha=0.5):
    """
    SCATTER_ELLIPSE(x, y, w=None, h=None, c='b', a=0.0)

    Make a scatter plot of x versus y with ellipses surrounding the
    center point.  w and h represent the width
    and height of the ellipse that surround each x,y coordinate.
    They are arrays of the same length as x or y.  c is
    a color and can be a single color format string or an length(x) array
    of intensities which will be mapped by the colormap jet. a is the
    angle or rotation in degrees of each ellipse (anti-clockwise). It is
    also an array of the same length as x or y or a single value to be
    iterated over all points.


    """
    if not axis_ob._hold:
        axis_ob.cla()

    if not iterable(a):
        a = [a] * len(x)

    if not iterable(alpha):
        alpha = [alpha] * len(x)
    if len(c) != len(x):
        raise ValueError('c and x are not equal lengths')
    if len(w) != len(x):
        raise ValueError('w and x are not equal lengths')

    if len(h) != len(x):
        raise ValueError('h and x are not equal lengths')
    if len(a) != len(x):
        raise ValueError('a and x are not equal lengths')
    # if len(alpha)!=len(x):
    #    raise ValueError, 'alpha and x are not equal lengths'
    patches = []
    for thisX, thisY, thisW, thisH, thisC, thisA, thisAl in \
            zip(x, y, w, h, c, a, alpha):
        ellip = Ellipse((thisX, thisY), width=thisW, height=thisH,
                        angle=thisA)

        ellip.set_facecolor(thisC)
        ellip.set_alpha(thisAl)
        axis_ob.add_patch(ellip)
        patches.append(ellip)
    axis_ob.autoscale_view()
    return axis_ob
