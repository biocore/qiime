#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

from qiime.util import parse_command_line_parameters, make_option
from qiime.normalize_table import normalize_CSS, normalize_DESeq2, multiple_file_normalize_CSS, multiple_file_normalize_DESeq2, algorithm_list

import os

__author__ = "Sophie Weiss"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Sophie Weiss"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Sophie Weiss"
__email__ = "sophie.sjw@gmail.com"

script_info = {}
script_info['brief_description'] = """Matrix normalization alternatives to rarefaction"""
script_info['script_description'] = \
"""To perform many downstream analyses after OTU picking (besides
metagenomeSeq's fitZIG and DESeq OTU differential abundance testing), the OTU
matrix must be normalized to account for uneven column (sample) sums that are a
result of most modern sequencing techniques.  These methods attempt to correct
for compositionality too.  Rarefying throws away some data by rarefying to a
constant sum and throwing away extremely low depth samples.

Even with these new normalization techniques, we would recommend throwing away
low depth samples (e.g. less that 1000 sequences/sample).  DESeq/DESeq2 outputs
negative values for lower abundant OTUs as a result of its log transformation.
For most ecologically useful metrics (e.g. UniFrac/Bray Curtis) this presents
problems. No good solution exists at the moment for this issue.  Note that one
is added to the matrix to avoid log(0).  It has been shown that clustering
results can be highly dependent upon the choice of the pseudocount (e.g. should
it be 0.01 instead of 1?), for more information read Costea, P. et al. (2014)
"A fair comparison", Nature Methods.

DESeq/DESeq2 can also have a very slow runtime, especially for larger datasets.
In this script, we implement DESeq2's variance stabilization technique. If you do use these
alternatives to rarefying, we would recommend metagenomeSeq's CSS (cumulative sum
scaling) transformation for those metrics that are abundance-based.  It is not
recommended to use these new methods with presence/absence metrics, for example
binary Jaccard or unweighted UniFrac.

For more on metagenomeSeq's CSS, please see Paulson, JN, et al. 'Differential
abundance analysis for microbial marker-gene surveys' Nature Methods 2013.  For DESeq
please see Anders S, Huber W. 'Differential expression analysis for sequence
count data.' Genome Biology 2010.  For DESeq2 please read Love, MI et al.
'Moderated estimation of fold change and dispersion for RNA-Seq data
with DESeq2,' Genome Biology 2014.  If you use these methods, please CITE the
appropriate reference as well as QIIME.  For any of these methods, clustering by
sequence depth MUST BE CHECKED FOR as a confounding variable, e.g. by coloring
by sequences/sample on a PCoA plot and testing for correlations between
taxa abundances and sequencing depth with e.g. adonis in compare_categories.py,
or observation_metadata_correlation.py.

Note: If the input BIOM table contains observation metadata (e.g., taxonomy
metadata for each OTU), this metadata will not be included in the output
normalized BIOM table when using DESeq2. When using CSS the taxonomy metadata
will be included in the output normalized table but it may not be in the same
format as the input table (e.g., "NA" will be added for missing taxonomic
levels). This discrepancy occurs because the underlying R packages used to
perform the normalization store taxonomy metadata in a different format.

As a workaround, the "biom add-metadata" command can be used to add the
original observation metadata to the output normalized table if desired. For
example, to include the original taxonomy metadata on the output normalized
table, "biom add-metadata" can be used with the representative sequence
taxonomic assignment file output by assign_taxonomy.py.

"""

script_info['script_usage']=[]
script_info['script_usage'].append(
      ("Single File CSS Matrix Normalization",
       "Normalize a raw (non-normalized/non-rarefied) otu_table.biom using CSS:",
       "%prog -i otu_table.biom -a CSS -o CSS_normalized_otu_table.biom"))
script_info['script_usage'].append(
       ("Single File DESeq2 Matrix Normalization",
       "Normalize a raw (non-normalized/non-rarefied) otu_table.biom using DESeq2:",
       "%prog -i otu_table.biom -a DESeq2 -o DESeq2_normalized_otu_table.biom"))
script_info['script_usage'].append(
      ("Multiple File Matrix Normalization",
       "Normalize a folder of raw (non-normalized/non-rarefied) otu tables using e.g. DESeq2:",
       "%prog -i otu_tables/ -a DESeq2 -o normalized_tables/"))

script_info['output_description']= \
"""BIOM table with normalized counts."""
script_info['required_options']=[]
script_info['optional_options']=[
make_option('-i', '--input_path', type='existing_path',
            help='path to the input BIOM file (e.g., the output '
            'from OTU picking) or directory containing input BIOM files '
            'for batch processing [REQUIRED if not passing -l]'),
make_option('-o', '--out_path', type='new_path',
            help='output filename for single file operation, or output '
            'directory for batch processing [REQUIRED if not passing -l]'),
make_option('-s', '--output_CSS_statistics', default=False,
            action='store_true', help='output CSS statistics file. This '
            'will be a directory for batch processing, and a filename for '
            'single file operation [default: %default]'),
make_option('-z', '--DESeq_negatives_to_zero', default=False,
            action='store_true', help='replace negative numbers produced by '
            'the DESeq normalization technique with zeros [default: '
            '%default]'),
make_option('-a', '--algorithm', default='CSS', type='choice',
             choices=algorithm_list(), help='normalization algorithm to apply to input '
             'BIOM table(s). [default: %default]' + ' Available options are: '
             '%s' % ', '.join(algorithm_list())),
make_option('-l', '--list_algorithms', action='store_true', default=False,
             help='show available normalization algorithms and exit '
             '[default: %default]'),
         ]
script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    input_path = opts.input_path
    out_path = opts.out_path
    output_CSS_statistics = opts.output_CSS_statistics
    DESeq_negatives_to_zero = opts.DESeq_negatives_to_zero
    algorithm = opts.algorithm
    list_algorithms = opts.list_algorithms


    if list_algorithms:
        print 'Available normalization algorithms are:\n%s' % ', '.join(algorithm_list())
    else:
        almost_required_options = ['input_path', 'out_path']
        for option in almost_required_options:
            if getattr(opts, option) is None:
                option_parser.error('Required option --%s omitted.' % option)
        if algorithm == 'CSS':
            if os.path.isdir(input_path):
                multiple_file_normalize_CSS(input_path, out_path, output_CSS_statistics)
            elif os.path.isfile(input_path):
                normalize_CSS(input_path, out_path, output_CSS_statistics)
            else:
                # it shouldn't be possible to get here
                option_parser.error("Unknown input type: %s" % input_path)
        elif algorithm == 'DESeq2':
            if os.path.isdir(input_path):
                multiple_file_normalize_DESeq2(input_path, out_path, DESeq_negatives_to_zero)
            elif os.path.isfile(input_path):
                normalize_DESeq2(input_path, out_path, DESeq_negatives_to_zero)
            else:
                # it shouldn't be possible to get here
                option_parser.error("Unknown input type: %s" % input_path)
        else:
            # it shouldn't be possible to get here
            option_parser.error("Unknown normalization algorithm: %s" % algorithm)

if __name__ == "__main__":
    main()
