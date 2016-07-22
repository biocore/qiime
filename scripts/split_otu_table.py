#!/usr/bin/env python
# File created on 15 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Antonio Gonzalez Pena",
               "Yoshiki Vazquez Baeza", "Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from os.path import split, splitext, join
from biom import load_table
from biom.table import Table
from numpy import array

from qiime.util import (parse_command_line_parameters, make_option, create_dir,
                        write_biom_table)
from qiime.parse import parse_mapping_file
from qiime.format import format_mapping_file
from qiime.split import (make_non_empty_sample_lists, subset_mapping_data)

script_info = {}
script_info['brief_description'] = "Split a biom table into one table per combination of values found in the specified fields in the mapping file."
script_info['script_description'] = """\
This script splits a biom table based on the cartesian product of the values
found in the mapping fields specified. It accepts any number of mapping fields
to split on. As an example assume the following was your mapping file data:

SampleID       Color       Habitat       Age
S1             Red         Stream        10
S2             Blue        Stream        20
S3             Blue        Lake          30
S4             Red         Stream        30

If we wanted to split a corresponding biom table by the 'Color' and 'Habitat'
fields simultanesouly, we would return 3 biom tables with the following samples
corresponding to the following groups:

(S1, S4): (Red, Stream)
(S2): (Blue, Stream)
(S3): (Blue, Lake)

Combinations which would result in no samples -- in our example (Red, Lake) -- 
are excluded and do not produce (empty) biom tables. The script optionally
produces split mapping files as well. 

The naming convention for split files is (assuming two fields):

input_table.biom -> input_table__field1_value1_field2_value2__.biom
input_mapping.txt -> input_mapping__field1_value1_field2_value2__.txt

So, from our example above:

input_table.biom -> (input_table__Color_Red_Habitat_Stream__.biom,
                     input_table__Color_Blue_Habitat_Stream__.biom,
                     input_table__Color_Blue_Habitat_Lake__.biom)
"""

script_info['script_usage'] = [\
    ("",
     "Split otu_table.biom into per-study OTU tables, and store the results in ./per_study_otu_tables/",
     "%prog -i otu_table.biom -m Fasting_Map.txt -f Treatment -o per_study_otu_tables"),
    ("",
     "Split otu_table.biom into multiple biom tables based on the Treatment and Color of the samples",
     "%prog -i otu_table.biom -m Fasting_Map.txt -f Treatment,Color -o ./per_study_otu_tables/")
    ]

script_info['output_description'] = ""

script_info['required_options'] = [
    make_option(
        '-i',
        '--biom_table_fp',
        type="existing_filepath",
        help='The input biom table file path.'),
    make_option(
        '-m',
        '--mapping_fp',
        type="existing_filepath",
        help='The mapping file path.'),
    make_option(
        '-f',
        '--fields',
        type='string',
        help="Mapping columns to split biom table on, comma separated."),
    make_option(
        '-o',
        '--output_dir',
        type="new_dirpath",
        help='File path to the output directory to be created.')
    ]

script_info['optional_options'] = [\
    make_option('--suppress_mapping_file_output', action='store_true', 
                default=False, help='Do not write out split mapping files.')
    ]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    biom_table_fp = opts.biom_table_fp
    mapping_fp = opts.mapping_fp
    fields = opts.fields.split(',')
    output_dir = opts.output_dir
    suppress_mf = opts.suppress_mapping_file_output
    # column_rename_ids = opts.column_rename_ids
    # include_repeat_cols = opts.include_repeat_cols

    bt = load_table(biom_table_fp)
    mdata, mheaders, mcomments = parse_mapping_file(mapping_fp)
    mdata = array(mdata)

    # check that biom file and mapping file have matching sample names. discard
    # those samples that do not appear in both. 
    shared_samples = list(set(mdata[:, 0]).intersection(bt.ids(axis='sample')))
    if len(shared_samples) == 0:
        raise ValueError('Mapping file and biom table share no samples.')
    elif len(shared_samples) == len(mdata[:, 0]):
        mdata = array(mdata)
    else:
        # we want to preserve the order of the samples in the biom table
        ss_bt_order = [s for s in bt.ids(axis='sample') if s in
                       shared_samples]
        bt = bt.filter(ss_bt_order, axis='sample', inplace=True)
        mdata = subset_mapping_data(mdata, shared_samples)
    # check that headers in mapping data
    if not all([i in mheaders for i in fields]):
        raise ValueError('One or more of the specified fields was not found ' +\
                         'in the mapping file.')

    # create output directory and create base names
    create_dir(output_dir)
    mf_base_name = join(output_dir, splitext(split(mapping_fp)[1])[0])
    bt_base_name = join(output_dir, splitext(split(biom_table_fp)[1])[0])

    # run code and append output
    sample_groups, value_groups = make_non_empty_sample_lists(fields, mheaders,
                                                              mdata)

    for sg, vg in zip(sample_groups, value_groups):
        name_base = '__' + '%s_%s_' * len(vg) + '_'
        name_tmp = []
        for f, v in zip(fields, vg):
            name_tmp.extend([f, v])
        nb = name_base % tuple(name_tmp)

        tmp_mf_data = subset_mapping_data(mdata, sg)
        tmp_mf_str = format_mapping_file(mheaders, tmp_mf_data, mcomments)
        write_biom_table(bt.filter(sg, axis='sample', inplace=False),
                         bt_base_name + nb + '.biom')
        
        if not suppress_mf:
            o = open(mf_base_name + nb + '.txt', 'w')
            o.writelines(tmp_mf_str)
            o.close()

if __name__ == "__main__":
    main()
