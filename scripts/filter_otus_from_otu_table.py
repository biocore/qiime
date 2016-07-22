#!/usr/bin/env python
# File created on 21 Dec 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from itertools import izip
from numpy import inf, isinf
from skbio.parse.sequences import parse_fasta
from biom import load_table

from qiime.util import (parse_command_line_parameters, make_option,
                        write_biom_table, EmptyBIOMTableError)
from qiime.filter import filter_otus_from_otu_table

script_info = {}
script_info[
    'brief_description'] = "Filter OTUs from an OTU table based on their observation counts or identifier."
script_info['script_description'] = ""
script_info[
    'script_usage'] = [("Singleton filtering", "Discard all OTUs that are observed fewer than 2 times (i.e., singletons)", "%prog -i otu_table.biom -o otu_table_no_singletons.biom -n 2"),

                       ("Abundance filtering", "Discard all OTUs that are observed greater than 100 times (e.g., if you want to look at low abundance OTUs only)",
                        "%prog -i otu_table.biom -o otu_table_low_abundance.biom -x 100"),

                       ("Chimera filtering", "Discard all OTUs listed in chimeric_otus.txt (e.g., to remove chimeric OTUs from an OTU table)", "%prog -i otu_table.biom -o otu_table_non_chimeric.biom -e chimeric_otus.txt")]

script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i', '--input_fp', type="existing_filepath",
                help='the input otu table filepath in biom format'),
    make_option('-o', '--output_fp', type="new_filepath",
                help='the output filepath in biom format'),
]
script_info['optional_options'] = [
    make_option('--negate_ids_to_exclude',
                action='store_true', default=False,
                help='keep OTUs in otu_ids_to_exclude_fp rather than discard them [default:%default] '),
    make_option('-n',
                '--min_count',
                type='int',
                default=0,
                help="the minimum total observation count of an otu for that otu to be retained [default: %default]"),
    make_option('--min_count_fraction',
                type='float',
                default=0,
                help="fraction of the total observation (sequence) count to apply as the minimum total observation count of an otu for that otu to be retained. this is a fraction, not percent, so if you want to filter to 1%, you specify 0.01. [default: %default]"),
    make_option('-x',
                '--max_count',
                type='int',
                default=inf,
                help="the maximum total observation count of an otu for that otu to be retained [default: infinity]"),
    make_option('-s',
                '--min_samples',
                type='int',
                default=0,
                help="the minimum number of samples an OTU must be observed in for that otu to be retained [default: %default]"),
    make_option('-y',
                '--max_samples',
                type='int',
                default=inf,
                help="the maximum number of samples an OTU must be observed in for that otu to be retained [default: infinity]"),

    make_option('-e', '--otu_ids_to_exclude_fp',
                type='existing_filepath',
                default=None,
                help="file containing list of OTU ids to exclude: can be a text file with one id per line, a text file where id is the first value in a tab-separated line, or can be a fasta file (extension must be .fna or .fasta) [default: %default]")

]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    input_fp = opts.input_fp
    output_fp = opts.output_fp

    min_count = opts.min_count
    max_count = opts.max_count
    min_count_fraction = opts.min_count_fraction
    if min_count_fraction < 0. or min_count_fraction > 1.:
        option_parser.error("min_count_fraction must be between 0 and 1")
    if min_count != 0 and min_count_fraction != 0:
        option_parser.error(
            "cannot specify both min_count and min_count_fraction")

    min_samples = opts.min_samples
    max_samples = opts.max_samples

    otu_ids_to_exclude_fp = opts.otu_ids_to_exclude_fp
    negate_ids_to_exclude = opts.negate_ids_to_exclude

    if not (min_count != 0 or
            min_count_fraction != 0 or
            not isinf(max_count) or
            otu_ids_to_exclude_fp is not None or
            min_samples != 0 or not isinf(max_samples)):
        option_parser.error("No filtering requested. Must provide either "
                            "min counts, max counts, min samples, max samples, min_count_fraction, "
                            "or exclude_fp (or some combination of those).")

    otu_table = load_table(opts.input_fp)

    if min_count_fraction > 0:
        min_count = otu_table.sum() * min_count_fraction

    otu_ids_to_keep = set(otu_table.ids(axis='observation'))

    if otu_ids_to_exclude_fp:
        if otu_ids_to_exclude_fp.endswith('.fasta') or \
           otu_ids_to_exclude_fp.endswith('.fna'):
            otu_ids_to_exclude = set([id_.strip().split()[0]
                                      for id_, seq in parse_fasta(open(otu_ids_to_exclude_fp, 'U'))])
        else:
            otu_ids_to_exclude = set([l.strip().split('\t')[0]
                                      for l in open(otu_ids_to_exclude_fp, 'U')])

        otu_ids_to_keep -= otu_ids_to_exclude

    filtered_otu_table = filter_otus_from_otu_table(otu_table,
                                                    otu_ids_to_keep,
                                                    min_count,
                                                    max_count,
                                                    min_samples,
                                                    max_samples,
                                                    negate_ids_to_exclude)

    try:
        write_biom_table(filtered_otu_table, opts.output_fp)
    except EmptyBIOMTableError:
        option_parser.error(
            "Filtering resulted in an empty BIOM table. "
            "This indicates that no OTUs remained after filtering.")

if __name__ == "__main__":
    main()
