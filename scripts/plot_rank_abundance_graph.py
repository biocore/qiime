#!/usr/bin/env python
# File created on 16 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"
 

from optparse import make_option
from itertools import cycle

from matplotlib.pyplot import ylim, xlim, show, legend, \
    savefig
from cogent.app.util import get_tmp_filename

from qiime.util import parse_command_line_parameters, get_options_lookup, \
    parse_otu_table, create_dir
from qiime.colors import data_color_order
from qiime.plot_rank_abundance_curve import plot_rank_abundance_graph


options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "plot rank-abundance curve"
script_info['script_description'] = "Plot a set of rank-abundance graphs from an OTU table and a set of sample names. Multiple graphs will be plotted into the same figure, in order to allow for an easy comparison across samples."
script_info['script_usage'] = [("Single graph example",
                                "Plot the rank-abundance curve of one sample using a linear scale for the x_axis:",
                                "plot_rank_abundance_graph.py -i otu_table.txt  -s 'Sample1' -x -v"),
                               ("multiple graph example",
                                "Plot the rank-abundance curve of several sample:",
                                "plot_rank_abundance_graph.py -i otu_table.txt  -s 'Sample1,Sample3,Sample5' -x  -v"),
                               ("multiple graph example",
                                "Plot the rank-abundance curve of all samples in an OTU table:",
                                "plot_rank_abundance_graph.py -i otu_table.txt  -s '*' -x -f eps -v"),
                               ]

script_info['output_description']= ""

script_info['required_options'] = [
 options_lookup['otu_table_as_primary_input'],
 make_option('-s','--sample_name',help='name of the sample to plot. Use "*" to plot all.'),
 ]

#could basically allow all of matplotlib format here
file_types = ['pdf','svg','png','eps']

script_info['optional_options'] = [\

    make_option('-o','--output_dir',help='name of output directory. '
                +'[default: random]'),\

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
    
    sample_ids, otu_ids, otu_table, lineages = parse_otu_table(open(opts.otu_table_fp, "U"))

    if opts.output_dir:  
        output_dir = opts.output_dir
    else:
        output_dir = get_tmp_filename(tmp_dir='./', prefix='rank_abundance_', suffix='')
    create_dir(output_dir, fail_on_exist=True)
    
    if opts.verbose:
        log_fh = open(output_dir+"/plot_rank_abundance_log.txt",'w')
        log_fh.write("OTU table file: %s\n"% opts.otu_table_fp)
        log_fh.write("sample names: %s\n" % opts.sample_name)
        
    #figure out which samples to draw
    if opts.sample_name=='*':
        user_sample_names = sample_ids
    else:
        user_sample_names = opts.sample_name.split(',')
        if len(user_sample_names)<1:
            raise ValueError, "sample IDs must be comma separated list of "\
            +"sample names - found %s" % opts.sample_name 

    # do the actual drawing
    ax=None
    for sample_name,color in zip(user_sample_names, cycle(data_color_order)):
        try:
            index = sample_ids.index(sample_name)
        except ValueError:
            if opts.verbose:
                log_fh.write("Warning: Sample name %s not in OTU table - skipping." % sample_name)
            continue     
        ax = plot_rank_abundance_graph(otu_table[:,index], color= color,
                                       absolute=opts.absolute_counts,
                                       label=sample_name)
        ax.set_label(sample_name)
    
    if ax==None:
        #ax should be defined if at least one series has been drawn
        raise ValueError("No data series drawn. Check your OTU table and sample names")

    #settings for all series
    ax.grid()      
    ax.set_xlabel('Species rank')
    ax.set_ylabel('Relative abundance')

    if not opts.x_linear_scale:
        ax.set_xscale('log')
    if not opts.y_linear_scale:
        ax.set_yscale('log')
  
    if not opts.no_legend:
        legend()

    #build output fp    
    output_fp = output_dir+ "/rank_abundance"
    if len(user_sample_names) < 6:
        output_fp += '_'.join(user_sample_names) 
    output_fp += ".%s" % opts.file_type

    savefig(output_fp, format=opts.file_type)
    
    
if __name__ == "__main__":
    main()
