#!/usr/bin/env python
# File created on 14 Mar 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from qiime.util import get_first_metadata_entry, make_option
from os.path import split
from qiime.util import parse_command_line_parameters, get_options_lookup, create_dir
from qiime.format import format_biom_table
from biom.parse import parse_biom_table

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Script to split a single OTU table into multiple tables based on the taxonomy at some user-specified depth."
script_info['script_description'] = """
Note: If the input BIOM table has multiple metadata entries (e.g., taxonomy, pathways) associated with its observations (e.g., OTUs), this represents a one-to-many relationship between the observation and its metadata, which is not currently supported by this script. Only the first metadata entry will be used to split the table and the remaining metadata entries will be ignored. For example, if an observation in the BIOM table has more than one taxonomy assignment (e.g., a list of taxonomy assignments), only the *first* taxonomy assignment will be used to split the table. Support for one-to-many relationships is a new addition to QIIME and while support is currently very limited, these types of relationships will be better incorporated in future versions of the software.

Note that if --md_as_string is provided, this will only specify a one-to-one relationship, not one-to-many. One-to-many relationships only exist if the metadata is a list of lists of strings, as opposed to a list of strings denoting different taxonomic levels.
"""
script_info['script_usage'] = [("","Split seqs_otu_table.biom into taxon-specific OTU tables based on the third level in the taxonomy, and write the taxon-specific OTU tables to ./L3/","%prog -i otu_table.biom -L 3 -o ./L3/")]
script_info['output_description']= ""

script_info['required_options'] = [
 make_option('-i','--input_fp',
             type='existing_filepath',
             help='the input otu table in BIOM format'),
 make_option('-o',
             '--output_dir',
             type='new_dirpath',
             help='the output directory'),
 make_option('-L','--level',type='int',
             help='the taxonomic level to split at'),
]

script_info['optional_options'] = [
 make_option('--md_identifier',default='taxonomy', type='string',
             help='the relevant observation metadata key. Currently only '
             'one-to-one relationships are supported; please see the full '
             'script description for details [default: %default]'),
 make_option('--md_as_string',default=False,action='store_true',
             help='metadata is included as string [default: metadata is included as list]'),
]
script_info['version'] = __version__


def process_md_as_string(md,md_identifier,level):
    return '_'.join(md[md_identifier].split(';')[:level]).replace(' ','')

def process_md_as_list(md,md_identifier,level):
    return '_'.join(get_first_metadata_entry(md[md_identifier])[:level]).replace(' ','')

def split_otu_table_on_taxonomy_to_files(otu_table_fp,
                                         level,
                                         output_dir,
                                         md_identifier='taxonomy',
                                         md_processor=process_md_as_list):
    """ Split OTU table by taxonomic level, writing otu tables to output dir
    """
    results = []
    otu_table = parse_biom_table(open(otu_table_fp,'U'))
    create_dir(output_dir)
    
    def split_f(obs_md):
        try:
            result = md_processor(obs_md,md_identifier,level)
        except KeyError:
            raise KeyError,\
             "Metadata identifier (%s) is not associated with all (or any) observerations. You can modify the key with the md_identifier parameter." % md_identifier
        except TypeError:
            raise TypeError,\
             "Can't correctly process the metadata string. If your input file was generated from QIIME 1.4.0 or earlier you may need to pass --md_as_string."
        except AttributeError:
            raise AttributeError,\
             "Metadata category not found. If your input file was generated from QIIME 1.4.0 or earlier you may need to pass --md_identifier \"Consensus Lineage\"."
    
        return result
    
    for bin, sub_otu_table in otu_table.binObservationsByMetadata(split_f):
        output_fp = '%s/otu_table_%s.biom' % (output_dir,bin)
        output_f = open(output_fp,'w')
        output_f.write(format_biom_table(sub_otu_table))
        output_f.close()
        results.append(output_fp)
    return results

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if opts.md_as_string:
        md_processor=process_md_as_string
    else:
        md_processor=process_md_as_list
    
    split_otu_table_on_taxonomy_to_files(opts.input_fp,
                                         opts.level,
                                         opts.output_dir,
                                         opts.md_identifier,
                                         md_processor)


if __name__ == "__main__":
    main()
