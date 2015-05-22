#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

from qiime.util import parse_command_line_parameters, make_option
from qiime.differential_abundance import (DA_fitZIG, multiple_file_DA_fitZIG,
    DA_DESeq2, multiple_file_DA_DESeq2, algorithm_list)

import os

__author__ = "Sophie Weiss"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Sophie Weiss"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Sophie Weiss"
__email__ = "sophie.sjw@gmail.com"

script_info = {}
script_info['brief_description'] = "Identify OTUs that are differentially abundance across two sample categories"
script_info['script_description'] = \
"""OTU differential abundance testing is commonly used to identify OTUs that
differ between two mapping file sample categories (i.e. Palm and Tongue body
sites).  We would recommend having at least 5 samples in each category.  These 
methods can be used in comparison to group_significance.py on a rarefied matrix, 
and we would always recommend comparing the results of these approaches to the
rarefied/group_significance.py approaches.  We would also recommend treating the
differentially abundant OTUs identified by these (metagenomeSeq zero-inflated Gaussian,
or ZIG, and DESeq2 negative binomial Wald test) techniques with caution, as they assume a
distribution and are therefore parametric tests.  Parametric tests can do poorly if
the assumptions about the data are not met.  These tests are also newer techniques
that are less well tested compared to rarefying with a group_significance.py test.  

The input is a raw (not normalized, not rarefied) matrix having uneven column sums.
With these techniques, we would still recommend removing low depth samples (e.g. below
1000 sequences per sample), and low abundance/rare OTUs from the data set.  The DESeq2 method
should NOT be used if the fit line on the dispersion plot (one of the diagnostic plots
output by the -d, or --DESeq2_diagnostic_plots option) does not look smooth, there are big
gaps in the point spacings, and the fitted line does not look appropriate for the data.

DESeq2 is stronger at very small/smaller data sets, but the run-time beyond 100 
total samples becomes very long.  MetagenomeSeq's fitZIG is a better algorithm for larger
library sizes and over 50 samples per category (e.g. 50 Palm samples), the more the better. 
In simulation, these techniques have higher sensitivity, but sometimes higher false positive
rate compared to the non-parametric tests (e.g. Wilcoxon rank sum) in group_significance.py, 
especially with very uneven library sizes (starting at 2-3 fold difference).  In practice 
and with real data, we do not observe much of a difference between these results and 
the tests in group_significance.py.  

For more on these techniques please see McMurdie, P. and Holmes, S. 'Waste not want not 
why rarefying microbiome data is inadmissible.' PLoS Comp. Bio. 2014.  For more on 
metagenomeSeq and fitZIG, please read Paulson, JN, et al. 'Differential abundance analysis
for microbial marker-gene surveys.'  Nature Methods 2013.  For DESeq2/DESeq please read
Love, MI et al. 'Moderated estimation of fold change and dispersion for RNA-Seq data 
with DESeq2,' Genome Biology 2014.  Anders S, Huber W. 'Differential expression analysis 
for sequence count data.' Genome Biology 2010.  Additionally, you can also read the
vignettes for each of the techniques on the Bioconductor/R websites.  Also, if you use
these methods, please CITE the proper sources above (metagenomeSeq and DESeq2) as well as QIIME."""

script_info['script_usage'] = []
script_info['script_usage'].append((
    "Single File OTU Differential Abundance Testing with metagenomeSeq_fitZIG",
    """Apply metagenomeSeq_fitZIG differential OTU abundance testing to a """
    """raw (NOT normalized) BIOM table to test for differences in OTU """
    """abundance between samples in the Treatment:Control and """
    """Treatment:Fast groups.""",
    "%prog -i otu_table.biom -o diff_otus.txt -m map.txt -a metagenomeSeq_fitZIG -c Treatment -x Control -y Fast")
    )
script_info['script_usage'].append((
    "Single File OTU Differential Abundance Testing with DESeq2_nbinom",
    """Apply DESeq2_nbinom differential OTU abundance testing to a """
    """raw (NOT normalized) BIOM table to test for differences in OTU """
    """abundance between samples in the Treatment:Control and """
    """Treatment:Fast groups, including output of plots.""",
    "%prog -i otu_table.biom -o diff_otus.txt -m map.txt -a DESeq2_nbinom -c Treatment -x Control -y Fast -d")
    )
script_info['script_usage'].append((
    "Multiple File OTU Differential Abundance Testing with metagenomeSeq_fitZIG",
    """Apply metagenomeSeq_fitZIG differential OTU abundance testing to a """
    """folder of raw (NOT normalized) BIOM tables to test for differences """
    """in OTU abundance between samples in the Treatment:Control and """
    """Treatment:Fast groups.""",
    "%prog -i otu_tables/ -o diff_otus/ -m map.txt -a metagenomeSeq_fitZIG -c Treatment -x Control -y Fast")
    )
script_info['output_description']= "The resulting output OTU txt file contains a list of all the OTUs in the input matrix, along with their associated statistics and FDR p-values."
script_info['required_options']=[
]
script_info['optional_options']=[
    make_option('-i', '--input_path', type='existing_path',
    help='path to the input BIOM file (e.g., the output '
    'from OTU picking) or directory containing input BIOM files '
    'for batch processing [REQUIRED if not passing -l]'),
    make_option('-o', '--out_path', type='new_path',
    help='output filename for single file operation, or output '
    'directory for batch processing [REQUIRED if not passing -l]'),
make_option('-a', '--algorithm', default='metagenomeSeq_fitZIG', type='choice',
    choices=algorithm_list(), help='differential abundance algorithm to '
    'apply to input BIOM table(s) [default: %default]' + ' Available options are: '
    '%s' % ', '.join(algorithm_list())),
make_option('-m', '--mapping_file_path', type='existing_filepath',
    help='path to mapping file [REQUIRED if not passing -l]'),
make_option('-c', '--mapping_file_category',
    help='mapping file category [REQUIRED if not passing -l]'),
make_option('-x', '--mapping_file_subcategory_1',
    help='mapping file subcategory [REQUIRED if not passing -l]'),
make_option('-y', '--mapping_file_subcategory_2',
    help='mapping file subcategory [REQUIRED if not passing -l]'),
make_option('-l', '--list_algorithms', action='store_true', default=False,
    help='show available differential abundance algorithms and exit '
    '[default: %default]'),
make_option('-d', '--DESeq2_diagnostic_plots', default=False,
    action='store_true', help='show a MA plot - y axis: log2 fold change, '
    'x axis: average size factor normalized OTU value. Also show a Dispersion '
    'Estimate plot - visualize the fitted dispersion vs. mean relationship '
    '[default: %default]'),
 ]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    input_path = opts.input_path
    out_path = opts.out_path
    algorithm = opts.algorithm
    mapping_fp = opts.mapping_file_path
    mapping_category = opts.mapping_file_category
    subcategory_1 = opts.mapping_file_subcategory_1
    subcategory_2 = opts.mapping_file_subcategory_2
    list_algorithms = opts.list_algorithms
    DESeq2_diagnostic_plots = opts.DESeq2_diagnostic_plots

    if list_algorithms:
        print 'Available differential abundance algorithms are:\n%s' % ', '.join(algorithm_list())
    else:
        almost_required_options = ['input_path', 'out_path', 'mapping_file_path', 'mapping_file_category', 'mapping_file_subcategory_1', 'mapping_file_subcategory_2']
        for option in almost_required_options:
            if getattr(opts, option) is None:
                option_parser.error('Required option --%s omitted.' % option)        
        if algorithm == 'metagenomeSeq_fitZIG':
            if os.path.isdir(input_path):
                multiple_file_DA_fitZIG(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2)
            elif os.path.isfile(input_path):
                DA_fitZIG(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2)
            else:
                # it shouldn't be possible to get here
                option_parser.error("Unknown input type: %s" % input_path)        
        elif algorithm == 'DESeq2_nbinom':
            if os.path.isdir(input_path):
                multiple_file_DA_DESeq2(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots)
            elif os.path.isfile(input_path):
                DA_DESeq2(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots)
            else:
                # it shouldn't be possible to get here
                option_parser.error("Unknown input type: %s" % input_path)
        else:
            # it shouldn't be possible to get here
            option_parser.error("Unknown normalization algorithm: %s" % algorithm)

if __name__ == "__main__":
    main()
