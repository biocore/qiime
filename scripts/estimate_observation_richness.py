#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import join
from biom import load_table
from skbio.util import create_dir
from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option)

from qiime.estimate_observation_richness import (
    Chao1MultinomialPointEstimator,
    ObservationRichnessEstimator)

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = ("Estimates the observation (e.g., OTU) richness "
                            "of samples in a BIOM table")
script_info['script_description'] = """
This script provides estimates of the observation (e.g., OTU) richness (i.e.
number of observations) given a sampling depth (i.e. number of
individuals/sequences per sample). Estimators are provided for both
interpolation/rarefaction and extrapolation.

Interpolation/rarefaction applies when the richness is estimated for a
*smaller* number of individuals than the original number of individuals in that
sample. We refer to this original sampling depth as the "reference sampling
depth" or "reference sample size".

Extrapolation applies when the richness is estimated for a *larger* number of
individuals than the reference sample size.

This script currently only provides a single unified estimation model for
interpolation and extrapolation. This model is the individual-based multinomial
model, which uses Chao1 to estimate the full richness of the sample. Please
refer to Colwell et al. (2012) for more details; equations 4, 5, 9, 10, 15a,
and 15b are used in this script.

For each interpolation/extrapolation point, the estimate, its unconditional
standard error, and confidence interval are reported. The script currently only
outputs this information to a table, which can be easily viewed in a program
such as Excel. Other output formats, such as plots, may be added in the future.

If an estimate is reported as "N/A", not enough information was present to
compute an estimate. This can occur when extrapolating if a sample does not
contain any singletons or doubletons, or if there is exactly one singleton and
no doubletons. A singleton is defined as an observation with exactly one
individual/sequence in the sample. A doubleton is defined as an observation
with exactly two individuals/sequences in the sample.

IMPORTANT: If you use the results of this script in any published works, please
be sure to cite the Colwell et al. (2012) paper, as well as QIIME (see
http://qiime.org for details).

In addition to Colwell et al. (2012), the following resources were extremely
useful while implementing and testing these estimators, so it is appropriate to
also acknowledge them here:

- Hsieh et al. (2013)
- Shen et al. (2003)
- Colwell (2013)

References:

Chao, A., N. J. Gotelli, T. C. Hsieh, E. L. Sander, K. H. Ma, R. K. Colwell, and A. M. Ellison 2013. Rarefaction and extrapolation with Hill numbers: a unified framework for sampling and estimation in biodiversity studies, Ecological Monographs (under revision).

Colwell, R. K. 2013. EstimateS: Statistical estimation of species richness and shared species from samples. Version 9. User's Guide and application published at: http://purl.oclc.org/estimates.

Colwell, R. K., A. Chao, N. J. Gotelli, S. Y. Lin, C. X. Mao, R. L. Chazdon, and J. T. Longino. 2012. Models and estimators linking individual-based and sample-based rarefaction, extrapolation and comparison of assemblages. Journal of Plant Ecology 5:3-21.

Hsieh, T. C., K. H. Ma, and A. Chao. 2013. iNEXT online: interpolation and extrapolation (Version 1.0) [Software]. Available from http://chao.stat.nthu.edu.tw/inext/.

Shen T-J, Chao A, Lin C- F. Predicting the number of new species in further taxonomic sampling. Ecology 2003;84:798-804.
"""

script_info['script_usage'] = []
script_info['script_usage'].append((
    "Interpolation and extrapolation of richness",
    "Estimate the richness of each sample in the input BIOM table using the "
    "default sampling depth range, which includes interpolation and "
    "extrapolation.",
    "%prog -i otu_table.biom -o estimates_out"))

script_info['output_description'] = """
A single file containing tabular data in TSV format is created in the output
directory. Other output formats may be added in the future.
"""

script_info['required_options'] = [
    make_option('-i', '--otu_table_fp', type='existing_filepath',
                help='path to the input BIOM table (e.g., the output from '
                'make_otu_table.py). IMPORTANT: This table should contain '
                'observation *counts* (integers), NOT relative abundances '
                '(fractions)'),
    options_lookup['output_dir']
]
script_info['optional_options'] = [
    make_option('-m', '--min', type='int',
                help='the number of individuals (e.g. sequences) per sample '
                'to start performing estimations at [default: %default]',
                default=1),
    make_option('-x', '--max', type='int',
                help='the number of individuals (e.g. sequences) per sample '
                'to stop performing estimations at. By default, the base '
                'sample size will be used, which is defined in Chao et al. '
                '(2013) as "double the smallest reference sample size or the '
                'maximum reference sample size, whichever is larger '
                '[default: base sample size]', default=None),
    make_option('-n', '--num_steps', type='int',
                help='the number of steps to make between -m/--min and '
                '-x/--max. Increasing this number will result in smoother '
                'curves, but will also increase the amount of time needed to '
                'run the script. Note that the reference sample size for each '
                'sample will be included if it does not fall within the '
                'min/max/num_steps range [default: %default]',
                default=10),
    make_option('-c', '--confidence_level', type='float',
                help='the confidence level of the unconditional confidence '
                'interval for each estimate. Must be a value between 0 and 1 '
                '(exclusive). For example, a 95% unconditional confidence '
                'interval would be 0.95 [default: %default]', default=0.95)
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create the output dir if it doesn't already exist.
    output_dir = opts.output_dir
    try:
        create_dir(output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o/--output_dir option.")

    otu_table_fp = opts.otu_table_fp
    table = load_table(otu_table_fp)

    estimator = ObservationRichnessEstimator(table,
                                             Chao1MultinomialPointEstimator)
    results = estimator(opts.min, opts.max, opts.num_steps,
                        opts.confidence_level)

    out_fp = join(output_dir, 'estimates_table.txt')
    with open(out_fp, 'w') as out_f:
        results.toTable(out_f)


if __name__ == "__main__":
    main()
