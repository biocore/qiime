#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Justin Kuczynski",
               "Julia Goodrich", "Daniel McDonald", "Antonio Gonzalez Pena",
               "Jesse Stombaugh", "Jose Carlos Clemente Litran",
               "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Development"
 
from qiime.util import parse_command_line_parameters
from qiime.util import make_option,get_options_lookup,create_dir
from sys import stdout, stderr
from qiime.parse import parse_mapping_file
from qiime.format import write_add_taxa_summary_mapping, format_biom_table
from os.path import split,splitext,join
from biom.parse import parse_biom_table
from qiime.summarize_taxa import (make_summary, add_summary_mapping,
                                  ONE_TO_MANY_TYPES)

options_lookup = get_options_lookup()

# Define the upper/lower percentage defaults here so that we can check whether
# the user supplied -u/-l on the command line with -m (which is not supported).
lower_percentage_default = 1.0
upper_percentage_default = 0.0

script_info={}
script_info['brief_description']="""Summarize taxa and store results in a new table or appended to an existing mapping file."""
script_info['script_description']="""The summarize_taxa.py script provides summary information of the representation of taxonomic groups within each sample. It takes an OTU table that contains taxonomic information as input. This OTU table should contain absolute abundances (not relative abundances). The taxonomic level for which the summary information is provided is designated with the -L option. The meaning of this level will depend on the format of the taxon strings that are returned from the taxonomy assignment step. The taxonomy strings that are most useful are those that standardize the taxonomic level with the depth in the taxonomic strings. For instance, for the RDP classifier taxonomy, Level 2 = Domain (e.g. Bacteria), 3 = Phylum (e.g. Firmicutes), 4 = Class (e.g. Clostridia), 5 = Order (e.g. Clostridiales), 6 = Family (e.g. Clostridiaceae), and 7 = Genus (e.g. Clostridium). By default, the relative abundance of each taxonomic group will be reported, but the raw counts can be returned if -a is passed.

By default, taxa summary tables will be output in both classic (tab-separated) and BIOM formats. The BIOM-formatted taxa summary tables can be used as input to other QIIME scripts that accept BIOM files.
"""
script_info['script_usage']=[]

script_info['script_usage'].append(("""Examples:""","""Summarize taxa based at taxonomic levels 2, 3, 4, 5, and 6, and write resulting taxa tables to the directory "./tax" ""","""%prog -i otu_table.biom -o ./tax"""))

script_info['script_usage'].append(("""Examples:""","""Summarize taxa based at taxonomic levels 2, 3, 4, 5, and 6, and write resulting mapping files to the directory "./tax" ""","""%prog -i otu_table.biom -o tax_mapping/ -m Fasting_Map.txt"""))

script_info['output_description']="""There are two possible output formats depending on whether or not a mapping file is provided with the -m option. If a mapping file is not provided, a table is returned where the taxonomic groups are each in a row and there is a column for each sample. If a mapping file is provided, the summary information will be appended to this file. Specifically, a new column will be made for each taxonomic group to which the relative abundances or raw counts will be added to the existing rows for each sample. The addition of the taxonomic information to the mapping file allows for taxonomic coloration of Principal coordinates plots in the 3d viewer. As described in the make_3d_plots.py section, principal coordinates plots can be dynamically colored based on any of the metadata columns in the mapping file. Dynamic coloration of the plots by the relative abundances of each taxonomic group can help to distinguish which taxonomic groups are driving the clustering patterns.
"""

script_info['required_options']= [\
    make_option('-i','--otu_table_fp', dest='otu_table_fp',
        help='Input OTU table filepath. Table should contain absolute '
        'abundances (i.e. raw counts), not relative abundances [REQUIRED]',
        type='existing_filepath'),
]
script_info['optional_options'] = [\
    make_option('-L','--level',default='2,3,4,5,6' , type='string',
        help='Taxonomic level to summarize by. [default: %default]'),
    make_option('-m','--mapping', 
        help='Input metadata mapping filepath. If supplied, then the taxon' +\
        ' information will be added to this file. This option is ' +\
        ' useful for coloring PCoA plots by taxon abundance or to ' +\
        ' perform statistical tests of taxon/mapping associations.',
        type='existing_filepath'),
    make_option('--md_identifier',default='taxonomy', type='string',
             help='the relevant observation metadata key [default: %default]'),
    make_option('--md_as_string',default=False,action='store_true',
             help='metadata is included as a string. The string will be split '
             'on the delimiter and used as one-to-one metadata (see '
             '-M/--one_to_many for more details) [default: metadata is '
             'included as a list of strings or a list of list of strings]'),
    make_option('-d','--delimiter',action='store',type='string',
        dest='delimiter',default=';', 
        help='Delimiter separating taxonomy levels. [default: %default]'),
    make_option('-M', '--one_to_many', type='choice',
        choices=ONE_TO_MANY_TYPES, default='first',
        help='the method of handling one-to-many relationships with metadata. '
        'If an observation (e.g., OTU) has more than one piece of metadata '
        'associated with it (identified via --md_identifier; e.g., multiple '
        'taxonomy assignments or gene pathways), this is a one-to-many '
        'relationship between the observation and its metadata. If "first" is '
        'supplied via this option, only the first piece of metadata (e.g., '
        'the first gene pathway) is used to summarize/collapse the '
        'observations, and all other pieces of metadata are ignored (this is '
        'the default behavior). If "add" is supplied via this option, the '
        'resulting summarized table will have the observation abundances '
        'counted for each of the pieces of metadata. This has the side effect '
        'of increasing the total count of observations in the table. '
        'If "divide" is supplied via this option, the observation abundances '
        'are divided evenly between each of the pieces of metadata. The total '
        'count of observations in the table will not be increased when this '
        'mode is used. One-to-many relationships will exist if the metadata '
        'is stored as a list of list of strings, where each inner list of '
        'strings is a single piece of metadata (e.g. taxonomy assignment or '
        'gene pathway). If the metadata is instead a list of strings, or the '
        'nested list only contains a single piece of metadata, this indicates '
        'a one-to-one relationship. Note that one-to-many and one-to-one '
        'relationships may both exist within a single BIOM table, as some '
        'observations may only be associated with a single piece of metadata '
        'while others may be associated with many. Therefore, this option '
        'only specifies how to handle one-to-many relationships and does not '
        'affect how one-to-one relationships are handled [default: %default]'),
    make_option('-a', '--absolute_abundance', action='store_true',\
        dest='absolute_abundance', default=False, \
        help='If present, the absolute abundance of the lineage in ' +\
        ' each sample is reported. By default, this script uses relative' +\
        ' abundance [default: %default]'),
    make_option('-l', '--lower_percentage', type='float', default=None,
        help='If present, taxa having higher *absolute* abundance are '
        'trimmed. IMPORTANT: these taxa will be trimmed based on *absolute '
        'abundance*. Thus, if -a/--absolute_abundance is not provided (the '
        'default), the taxa will be trimmed *before* the table is converted '
        'into relative abundances (so the relative abundance calculation will '
        'only consider the remaining taxa that were *not* trimmed). For '
        'example, to remove taxa that make up more than 5% of the total '
        'dataset you would pass 0.05. [default: no taxa are trimmed]'),
    make_option('-u', '--upper_percentage', type='float', default=None,
        help='If present, taxa having lower *absolute* abundance are '
        'trimmed. IMPORTANT: these taxa will be trimmed based on *absolute '
        'abundance*. Thus, if -a/--absolute_abundance is not provided (the '
        'default), the taxa will be trimmed *before* the table is converted '
        'into relative abundances (so the relative abundance calculation will '
        'only consider the remaining taxa that were *not* trimmed). For '
        'example, to remove taxa that make up less than 45% of the total '
        'dataset you would pass 0.45. [default: no taxa are trimmed]'),
    make_option('-t', '--transposed_output', action='store_true',\
        dest='transposed_output', default=False, \
        help='If present, the output will be written transposed from' +\
        ' the regular output. This is helpful in cases when you want to' +\
        ' use Site Painter to visualize your data [default: %default]'),
    options_lookup['output_dir'],
    make_option('--suppress_classic_table_output', action='store_true',
        default=False, help='If present, the classic (TSV) format taxon table '
        'will not be created in the output directory. This option is ignored '
        'if -m/--mapping is present [default: %default]'),
    make_option('--suppress_biom_table_output', action='store_true',
        default=False, help='If present, the BIOM-formatted taxon table will '
        'not be created in the output directory. This option is ignored if '
        '-m/--mapping is present [default: %default]')
]
script_info['option_label']={'otu_table_fp':'OTU table filepath',
                             'output_fp': 'Output filepath',
                             'mapping':'QIIME-formatted mapping filepath',
                             'level':'Summarize level',
                             'delimiter': 'Taxonomic delimiter',
                             'absolute_abundance':'Use absolute abundance',
                             'lower_percentage':'Top % of taxa to remove',
                             'upper_percentage':'Bottom % of taxa to remove'}

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    lower_percentage = opts.lower_percentage
    upper_percentage = opts.upper_percentage
    otu_table_fp = opts.otu_table_fp
    otu_table = parse_biom_table(open(otu_table_fp, 'U'))
    delimiter = opts.delimiter
    mapping_fp = opts.mapping
    md_as_string = opts.md_as_string
    md_identifier = opts.md_identifier
    levels = opts.level.split(',')
    suppress_classic_table_output = opts.suppress_classic_table_output
    suppress_biom_table_output = opts.suppress_biom_table_output

    if (upper_percentage is not None or lower_percentage is not None) and \
            mapping_fp:
        raise ValueError("upper_percentage and lower_percentage cannot be used with mapping file")

    if upper_percentage is None:
        upper_percentage = upper_percentage_default
    if lower_percentage is None:
        lower_percentage = lower_percentage_default

    if upper_percentage < 0 or upper_percentage > 1.0:
        raise ValueError('upper_percentage should be between 0.0 and 1.0')
    
    if lower_percentage < 0 or lower_percentage > 1.0:
        raise ValueError('lower_percentage should be between 0.0 and 1.0')

    if mapping_fp:
        mapping_file = open(mapping_fp, 'U')
        mapping, header, comments = parse_mapping_file(mapping_file)
        
        # use the input Mapping file for producing the output filenames
        map_dir_path,map_fname=split(mapping_fp)
        map_basename,map_fname_ext=splitext(map_fname)
    else:
        if suppress_classic_table_output and suppress_biom_table_output:
            option_parser.error("Both classic and BIOM output formats were "
                                "suppressed.")

    # introduced output directory to will allow for multiple outputs
    if opts.output_dir:
        create_dir(opts.output_dir,False)
        output_dir_path=opts.output_dir
    else:
        output_dir_path='./'

    # use the input OTU table to produce the output filenames
    dir_path,fname=split(otu_table_fp)
    basename,fname_ext=splitext(fname)
    
    # Iterate over the levels and generate a summarized taxonomy for each
    for level in levels:
        if mapping_fp:
            #define output filename
            output_fname = join(output_dir_path,
                                map_basename+'_L%s.txt' % (level))

            summary, tax_order = add_summary_mapping(otu_table, 
                    mapping, int(level), md_as_string, md_identifier,
                    delimiter=delimiter, one_to_many=opts.one_to_many,
                    absolute_abundance=opts.absolute_abundance)

            write_add_taxa_summary_mapping(summary,tax_order,mapping,
                                           header,output_fname)
        else:
            # define the output filename. The extension will be added to the
            # end depending on the output format
            output_fname = join(output_dir_path, basename + '_L%s' % level)

            ts_table = make_summary(otu_table,
                                    int(level),
                                    upper_percentage,
                                    lower_percentage,
                                    md_as_string,
                                    md_identifier,
                                    delimiter=delimiter,
                                    one_to_many=opts.one_to_many,
                                    absolute_abundance=opts.absolute_abundance)

            if opts.transposed_output:
                ts_table = ts_table.transpose()
                obs_col_name = 'SampleID'
            else:
                obs_col_name = 'Taxon'

            if not suppress_classic_table_output:
                classic_table_str = ts_table.delimitedSelf(
                        observation_column_name=obs_col_name)
                with open(output_fname + '.txt', 'w') as output_f:
                    output_f.write(classic_table_str)
                    output_f.write('\n')

            if not suppress_biom_table_output:
                with open(output_fname + '.biom', 'w') as output_f:
                    output_f.write(format_biom_table(ts_table))


if __name__ == "__main__":
    main()
