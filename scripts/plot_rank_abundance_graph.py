#!/usr/bin/env python
# File created on 16 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"
 

from qiime.util import make_option

from qiime.util import get_tmp_filename

from qiime.plot_rank_abundance_graph import plot_rank_abundance_graphs
from qiime.util import parse_command_line_parameters, get_options_lookup, \
    create_dir
from biom.parse import parse_biom_table

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "plot rank-abundance curve"
script_info['script_description'] = "Plot a set of rank-abundance graphs from an OTU table and a set of sample names. Multiple graphs will be plotted into the same figure, in order to allow for an easy comparison across samples."
script_info['script_usage'] = [("Single graph example",
                                "Plot the rank-abundance curve of one sample using a linear scale for the x_axis:",
                                "%prog -i otu_table.biom  -s 'PC.354' -x -v -o single_plot/"),
                               ("multiple graph example",
                                "Plot the rank-abundance curve of several sample:",
                                "%prog -i otu_table.biom  -s 'PC.354,PC.481,PC.636' -x -v -o multi_plot/"),
                               ("multiple graph example",
                                "Plot the rank-abundance curve of all samples in an OTU table:",
                                "%prog -i otu_table.biom  -s '*' -x -f eps -v -o all_plot/"),
                               ]

script_info['output_description']= ""

script_info['required_options'] = [
 options_lookup['otu_table_as_primary_input'],
 make_option('-s','--sample_name',help='name of the sample to plot. Use "*" to plot all.'),
 make_option('-o','--output_dir',help='name of output directory',type='new_dirpath'),
 ]

#could basically allow all of matplotlib format here
file_types = ['pdf','svg','png','eps']

script_info['optional_options'] = [\

    make_option('-a','--absolute_counts', help='plot absolute abundance values instead of ' 
                +'relative [default: %default]', action='store_true', default=False),\

    make_option('-n','--no-legend', action='store_true', default=False,
                help='do not draw a legend [default: %default]'),\

    make_option('-x','--x_linear_scale',help='draw x axis in linear scale '
                +'[default: %default]',
                action='store_true', default=False),\

    make_option('-y','--y_linear_scale',help='draw y axis in linear scale '
                +'[default: %default]',
                action='store_true', default=False),

    make_option('-f','--file_type',help='save plot using this image type. Choice of '+
                ', '.join(file_types) +' [default: %default]',
                default='pdf')

]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)
    
    output_dir = opts.output_dir
    otu_table_fp = opts.otu_table_fp
    sample_name = opts.sample_name
    
    create_dir(output_dir)
    
    if opts.verbose:
        log_fh = open(output_dir+"/plot_rank_abundance_log.txt",'w')
        log_fh.write("OTU table file: %s\n"% otu_table_fp)
        log_fh.write("sample names: %s\n" % sample_name)
    else:
        log_fh=None
    
    otu_table = parse_biom_table(open(opts.otu_table_fp,"U"))
    
    plot_rank_abundance_graphs(opts.sample_name, otu_table,
                               output_dir, opts.file_type,
                               opts.absolute_counts, opts.x_linear_scale,
                               opts.y_linear_scale, opts.no_legend, log_fh)
    
if __name__ == "__main__":
    main()
