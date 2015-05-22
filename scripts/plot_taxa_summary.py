#!/usr/bin/env python
# File created on 19 Jan 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jesse Stombaugh", "Julia Goodrich", "Justin Kuczynski",
               "John Chase", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
"""
This script generates taxonomy charts
"""

from qiime.util import parse_command_line_parameters, get_qiime_project_dir
from qiime.util import make_option
from qiime.util import create_dir
from qiime.plot_taxa_summary import make_all_charts
from tempfile import mkdtemp
from qiime.colors import taxonomy_color_prefs_and_map_data_from_options
import re
import matplotlib
import os
import shutil

plot_filetype_choices = ['pdf', 'svg', 'png']

script_info = {}
script_info['brief_description'] = """Make taxaonomy summary charts based on\
 taxonomy assignment"""
script_info['script_description'] = """This script automates the construction\
 of pie, bar and area charts showing the breakdown of taxonomy by given levels.\
 The script creates an html file for each chart type for easy visualization. It\
 uses the taxonomy or category counts from summarize_taxa.py for combined\
 samples by level (-i) and user specified labels for each file passed in (-l).\
 Output will be written to the user specified folder (-o) the, where the\
 default is the current working directory. The user can also specify the number\
 of categories displayed for within a single pie chart, where the rest are\
 grouped together as the 'other category' using the (-n) option, default is 20.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Examples:""",
                                    """If you wish to run the code using default parameters, you must supply a\
 counts file (phylum.txt) along with the taxon level label (Phylum), the\
 type(s) of charts to produce, and an output directory, by using the following\
 command:""",
                                    """%prog -i phylum.txt -l phylum -c pie,bar,area -o phylum_charts/"""))

script_info['script_usage'].append(("""""",
                                    """If you want to make charts for multiple levels at a time\
 (phylum.txt,class.txt,genus.txt) use the following command:""",
                                    """%prog -i phylum.txt,class.txt,genus.txt -l Phylum,Class,Genus\
 -c pie,bar,area -o phylum_class_genus_charts/"""))

script_info['script_usage'].append(("""""",
                                    """Additionally, if you would like to display on a set number of taxa ("-n 10")\
 in the pie charts, you can use the following command:""",
                                    """%prog -i class.txt -l Class -c pie -n 10 -o class_pie_n10_charts/"""))

script_info['script_usage'].append(("""""",
                                    """If you would like to display generate pie charts for specific samples, i.e.\
 sample 'PC.636' and sample 'PC.635' that are in the counts file header, you\
 can use the following command:""",
                                    """%prog -i class.txt -l Class -b PC.636,PC.635 -o sample_charts/"""))

script_info['output_description'] = """The script generates an output folder,\
 which contains several files. For each pie chart there is a png and a pdf\
 file. The best way to view all of the pie charts is by opening up the file\
 taxonomy_summary_pie_chart.html."""

script_info['required_options'] = [
    # dest should equal long-form parameter names! Can you clean this up?
    # Also note that you don't need to pass type='string' - that's the default
    make_option('-i', '--counts_fname',
                help='Input comma-separated list of summarized taxa filepaths' +
                ' (i.e results from summarize_taxa.py) [REQUIRED]',
                type='existing_filepaths'),
]
script_info['optional_options'] = [
    # changed this from type='string' (default) to type='int'
    make_option('-l', '--labels',
                help='Comma-separated list of taxonomic levels (e.g.' +
                ' Phylum,Class,Order)  [default=%default]', default=None),
    make_option('-n', '--num_categories', dest='num_categories',
                help='The maximum number of taxonomies to show in each pie chart.' +
                ' All additional taxonomies are grouped into an "other" category.' +
                ' NOTE: this functionality only applies to the pie charts.' +
                ' [default: %default]', default=20, type='int'),
    make_option('-o', '--dir_path',
                help='Output directory',
                type='new_dirpath'),
    make_option('-b', '--colorby', dest='colorby', type='string',
                help='This is the categories to color by in the plots from the' +
                ' metadata mapping file. The categories must match the name of a ' +
                ' column header in the mapping file exactly and multiple categories' +
                ' can be list by comma separating them without spaces.' +
                ' [default=%default]'),
    make_option('-p', '--prefs_path',
                help='Input user-generated preferences filepath. NOTE: This is a' +
                ' file with a dictionary containing preferences for the analysis.' +
                ' The key taxonomy_coloring is used for the coloring.' +
                ' [default: %default]',
                type='existing_filepath'),
    make_option('-k', '--background_color',
                help='This is the background color to use in the plots' +
                ' (black or white) [default: %default]', default='white',
                type='choice', choices=['black', 'white'],),
    make_option('-d', '--dpi',
                help='This is the resolution of the plot. [default: %default]',
                type='int', default=80),
    make_option('-x', '--x_width',
                help='This is the width of the x-axis to use in the plots.' +
                ' [default: %default]', default=12, type='int'),
    make_option('-y', '--y_height',
                help='This is the height of the y-axis to use in the plots.' +
                ' [default: %default]', default=6, type='int'),
    make_option('-w', '--bar_width',
                help='This the width of the bars in the bar graph and should be a' +
                ' number between 0 and 1. NOTE: this only applies to the bar charts.' +
                ' [default: %default]', default=0.75, type='float'),
    make_option('-t', '--type_of_file', type='choice',
                help='This is the type of image to produce (i.e. ' +
                ','.join(plot_filetype_choices) + '). [default: %default]',
                choices=plot_filetype_choices, default='pdf'),
    make_option('-c', '--chart_type', type='multiple_choice',
                mchoices=['pie', 'bar', 'area'],
                help='This is the type of chart to plot (i.e. pie, bar or area).' +
                ' The user has the ability to plot multiple types, by using a' +
                ' comma-separated list (e.g. area,pie) [default: %default]',
                default='area,bar'),
    make_option('-r', '--resize_nth_label', type='int',
                help='Make every nth label larger than the other lables.' +
                ' This is for large area and bar charts where the font on the x-axis' +
                ' is small. This requires an integer value greater than 0.' +
                ' [default: %default]', default=0),
    make_option('-s', '--include_html_legend', action='store_true',
                dest='include_html_legend', default=False,
                help='Include HTML legend. If present, the writing of the legend' +
                ' in the html page is included. [default: %default]'),
    make_option('-a', '--label_type', type='choice',
                help='Label type ("numeric" or "categorical"). ' +
                ' If the label type is defined as numeric, the x-axis will be' +
                ' scaled accordingly. Otherwise the x-values will treated' +
                ' categorically and be evenly spaced [default: %default].',
                choices=['categorical', 'numeric'], default='categorical'),
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # get QIIME directory
    qiime_dir = get_qiime_project_dir()

    if not opts.counts_fname:
        option_parser.error("A list of input files must be specified")

    # get color preferences
    color_prefs, color_data, background_color, label_color = \
        taxonomy_color_prefs_and_map_data_from_options(opts)

    colorby = opts.colorby
    if colorby is None:
        colorby = []
        for c in color_data['counts'].values():
            colorby.extend(c[0])
    else:
        colorby = colorby.strip().strip("'").split(',')

    counts_fname = opts.counts_fname

    # Define labels to use
    labels = opts.labels

    if not opts.labels:
        new_labels = []
        # create an empty list since the user didn't specify labels
        for i in counts_fname:
            new_labels.append("")
        labels = ','.join(new_labels)

    data = [(label, f.strip())
            for f, label in zip(counts_fname, labels.split(","))]
    filepath = data[0][1]

    filename = filepath.strip().rpartition('/')[0]
    num_categories = int(opts.num_categories)
    if num_categories <= 0:
        raise ValueError('The number of categories has to be greater than 0!')

    # create directory path
    dir_path = os.getcwd()
    if opts.dir_path:
        dir_path = opts.dir_path
        try:
            create_dir(opts.dir_path)
        except OSError:
            pass

    # make javascript output directory
    javascript_path = os.path.join(dir_path, 'js')
    try:
        create_dir(javascript_path)
    except OSError:  # raised if dir exists
        pass

    # make raw_data output directory
    raw_data_path = os.path.join(dir_path, 'raw_data')
    try:
        create_dir(raw_data_path)
    except OSError:  # raised if dir exists
        pass

    # move javascript file to javascript output directory
    shutil.copyfile(os.path.join(qiime_dir, 'qiime', 'support_files',
                    'js/overlib.js'),
                    os.path.join(javascript_path, 'overlib.js'))

    # make css output directory
    css_path = os.path.join(dir_path, 'css')
    try:
        create_dir(css_path)
    except OSError:  # raised if dir exists
        pass

    # move css file to css output directory
    shutil.copyfile(os.path.join(qiime_dir, 'qiime', 'support_files',
                    'css/qiime_style.css'),
                    os.path.join(css_path, 'qiime_style.css'))

    # verify all parameters are valid
    plot_width = float(opts.x_width)
    if plot_width <= 0:
        raise ValueError('The width of the plot has to be greater than 0!')

    plot_height = float(opts.y_height)
    if plot_height <= 0:
        raise ValueError('The height of the plot has to be greater than 0!')

    bar_width = float(opts.bar_width)
    if bar_width <= 0 or bar_width > 1:
        raise ValueError(
            'The bar width of the plot has to be between 0 and 1!')

    dpi = float(opts.dpi)
    if dpi <= 0:
        raise ValueError('The dpi of the plot has to be greater than 0!')

    resize_nth_label = int(opts.resize_nth_label)
    if resize_nth_label < 0:
        raise ValueError('The resize_nth_label of the plot has to be greater\
 than 0!')

    generate_image_type = opts.type_of_file
    label_type = opts.label_type
    include_html_legend = opts.include_html_legend
    plots_to_make = opts.chart_type
    for chart_type in plots_to_make:

        # make pie chart output path
        charts_path = os.path.join(dir_path, 'charts')
        try:
            create_dir(charts_path)
        except OSError:  # raised if dir exists
            pass

        make_all_charts(data, dir_path, filename, num_categories,
                        colorby, args, color_data, color_prefs, background_color, label_color,
                        chart_type, generate_image_type, plot_width, plot_height, bar_width, dpi,
                        resize_nth_label, label_type, include_html_legend)


if __name__ == "__main__":
    main()
