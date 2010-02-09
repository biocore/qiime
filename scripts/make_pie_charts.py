#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.make_pie_charts import make_all_pie_charts, create_dir

script_description = """Creates 2D pie charts of using taxonomy counts from \
summarize_taxa.py or category counts from summarize_otu_by_cat.py"""

script_usage = """Create pie charts that are viewed in an html file for \
easy visualization of all of the pie charts. It uses the taxonomy or category\
counts for combined samples by level phylum.txt,class.txt,genus.txt (-i) and\
user specified labels for each file passed in (-l). Output will be in a \
randomly generated folder name within the user specified folder (-o) the \
default is the current working directory. There is also an option for making \
seperate pie charts for each sample (-s). User can also specify the number \
of categories displayed in pie charts the rest are grouped together as \
other category (-n) default is 20.

python make_pie_charts.py -i phylum.txt,class.txt,genus.txt -l \
phylum,class,genus -o ./webfiles

Create pie charts using counts for combined samples and individual \
samples:

python ~/code/Qiime/trunk/qiime/make_pie_charts.py -i \
phylum.txt,class.txt,genus.txt -l phylum,class,genus -o ./webfiles -s
"""

required_options = [\
make_option('-i', '--input_files', dest='counts_fname',\
	action='store',type='string',
	help='list of files with sample counts by taxonomy [REQUIRED]'),
make_option('-l', '--labels', dest='labels',action='store',type='string',
            help='list of labels for pie chart(i.e. Phylum,Class)[REQUIRED]')
]

optional_options = [\
make_option('-s', '--sample_flag', dest='do_sample',
	help='if True pie charts will be created for each sample',default=False,
                      action = 'store_true'),
make_option('-n', '--num', dest='num_categories', \
                help='name of file containing metadata [default: %default]', \
                      default='20'),
make_option('-o', '--dir-prefix', dest='dir_path',\
               help='directory prefix for all analyses')
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
    if not opts.counts_fname:
        option_parser.error("A list of input files must be specified")
    if not opts.labels:
        option_parser.error("A list of label names cooresponding to files must \
		be specified")
	
    dir_path = opts.dir_path
    if dir_path == './':
        dir_path = os.getcwd()
    do_sample = opts.do_sample
    counts_fname = opts.counts_fname
    labels = opts.labels
    data = [(label,f.strip()) \
            for f,label in zip(counts_fname.split(","),labels.split(","))]
    filepath=data[0][1]
    filename=filepath.strip().rpartition('/')[0]
    num_categories = int(opts.num_categories)
    make_all_pie_charts(data,dir_path,filename,num_categories, do_sample,args)

if __name__ == "__main__":
    main()
