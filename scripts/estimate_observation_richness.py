#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

from os.path import join

from biom.parse import parse_biom_table

from cogent.util.misc import create_dir

from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option)

from qiime.estimate_observation_richness import (
        Chao1MultinomialPointEstimator,
        ObservationRichnessEstimator)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Estimates the observation richness of samples in a BIOM table"
script_info['script_description'] = """
This script provides estimates of the observation richness (i.e. observation count)

References:

Colwell, R. K., A. Chao, N. J. Gotelli, S. Y. Lin, C. X. Mao, R. L. Chazdon, and J. T. Longino. 2012. Models and estimators linking individual-based and sample-based rarefaction, extrapolation and comparison of assemblages. Journal of Plant Ecology 5:3-21.

Chao, A., N. J. Gotelli, T. C. Hsieh, E. L. Sander, K. H. Ma, R. K. Colwell, and A. M. Ellison 2013. Rarefaction and extrapolation with Hill numbers: a unified framework for sampling and estimation in biodiversity studies, Ecological Monographs (under revision).

Hsieh, T. C., K. H. Ma, and A. Chao. 2013. iNEXT online: interpolation and extrapolation (Version 1.0) [Software]. Available from http://chao.stat.nthu.edu.tw/inext/.

Colwell, R. K. 2013. EstimateS: Statistical estimation of species richness and shared species from samples. Version 9. User's Guide and application published at: http://purl.oclc.org/estimates.
"""

script_info['script_usage'] = []
script_info['script_usage'].append(("", "", ""))

script_info['output_description'] = ""

script_info['required_options'] = [
    options_lookup['otu_table_as_primary_input'],
    options_lookup['output_dir']
]
script_info['optional_options'] = [
    make_option('-m', '--min', type='int',
        help='the number of individuals (e.g. sequences) per sample to start performing estimations at [default: %default]', default=1),
    make_option('-x', '--max', type='int',
        help='the number of individuals (e.g. sequences) per sample to stop performing estimations at. By default, the base sample size will be used, which is defined in Chao et al. (2013) as "double the smallest reference sample size or the maximum reference sample size, whichever is larger" [default: base sample size]', default=None),
    make_option('-n', '--num_steps', type='int', help='the number of steps to make between -m/--min and -x/--max.  Increasing this number will result in smoother curves, but will also increase the amount of time needed to run the script. Note that reference sample size for each sample will be included if it does not fall within the min/max/num_steps range [default: %default]', default=10),
    make_option('-c', '--confidence_level', type='float',
        help='the confidence level of the unconditional confidence interval for each estimate. Must be a value between 0 and 1 (exclusive). For example, a 95% unconditional confidence interval would be 0.95 [default: %default]', default=0.95)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create the output dir if it doesn't already exist.
    out_dir = opts.output_dir
    try:
        create_dir(out_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o/--output_dir option.")

    table_fp = opts.otu_table_fp
    with open(table_fp, 'U') as table_f:
        table = parse_biom_table(table_f)

    estimator = ObservationRichnessEstimator(table,
                                             Chao1MultinomialPointEstimator)
    results = estimator(opts.min, opts.max, opts.num_steps,
                        opts.confidence_level)

    out_fp = join(out_dir, 'estimates_table.txt')
    with open(out_fp, 'w') as out_f:
        results.toTable(out_f)


if __name__ == "__main__":
    main()
