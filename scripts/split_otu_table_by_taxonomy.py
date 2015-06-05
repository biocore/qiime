#!/usr/bin/env python
# File created on 14 Mar 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from os.path import split

from biom import load_table

from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        create_dir, make_option, write_biom_table)

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = "Script to split a single OTU table into multiple tables based on the taxonomy at some user-specified depth."
script_info['script_description'] = ""
script_info['script_usage'] = [(
    "",
    "Split seqs_otu_table.biom into taxon-specific OTU tables based on the third level in the taxonomy, and write the taxon-specific OTU tables to ./L3/",
    "%prog -i otu_table.biom -L 3 -o ./L3/")]
script_info['output_description'] = ""

script_info['required_options'] = [
    make_option('-i', '--input_fp',
                type='existing_filepath',
                help='the input otu table in BIOM format'),
    make_option('-o',
                '--output_dir',
                type='new_dirpath',
                help='the output directory'),
    make_option('-L', '--level', type='int',
                help='the taxonomic level to split at'),
]

script_info['optional_options'] = [
    make_option('--md_identifier', default='taxonomy', type='string',
                help='the relevant observation metadat key [default: %default]'),
    make_option('--md_as_string', default=False, action='store_true',
                help='metadata is included as string [default: metadata is '
                'included as list]'),
]
script_info['version'] = __version__


def process_md_as_string(md, md_identifier, level):
    return '_'.join(md[md_identifier].split(';')[:level]).replace(' ', '')


def process_md_as_list(md, md_identifier, level):
    return '_'.join(md[md_identifier][:level]).replace(' ', '')


def split_otu_table_on_taxonomy_to_files(otu_table_fp, level, output_dir,
                                         md_identifier='taxonomy',
                                         md_processor=process_md_as_list):
    """ Split OTU table by taxonomic level, writing otu tables to output dir
    """
    results = []
    otu_table = load_table(otu_table_fp)
    create_dir(output_dir)

    def split_f(id_, obs_md):
        try:
            result = md_processor(obs_md, md_identifier, level)
        except KeyError:
            raise KeyError("Metadata identifier (%s) is not associated with "
                           "all (or any) observerations. You can modify the "
                           "key with the md_identifier parameter." %
                           md_identifier)
        except TypeError:
            raise TypeError("Can't correctly process the metadata string. If "
                            "your input file was generated from QIIME 1.4.0 or"
                            " earlier you may need to pass --md_as_string.")
        except AttributeError:
            raise AttributeError("Metadata category not found. If your input "
                                 "file was generated from QIIME 1.4.0 or "
                                 "earlier you may need to pass --md_identifier"
                                 " \"Consensus Lineage\".")

        return result

    for bin, sub_otu_table in otu_table.partition(split_f, axis='observation'):
        output_fp = '%s/otu_table_%s.biom' % (output_dir, bin)
        write_biom_table(sub_otu_table, output_fp)

        results.append(output_fp)
    return results


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.md_as_string:
        md_processor = process_md_as_string
    else:
        md_processor = process_md_as_list

    split_otu_table_on_taxonomy_to_files(opts.input_fp, opts.level,
                                         opts.output_dir, opts.md_identifier,
                                         md_processor)


if __name__ == "__main__":
    main()
