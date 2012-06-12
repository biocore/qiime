#!/usr/bin/env python
# File created on 12 Jun 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


from biom.parse import parse_biom_table
from biom.exception import TableException
from qiime.util import parse_command_line_parameters, make_option
from qiime.format import format_biom_table
from qiime.core_microbiome import filter_table_to_core

script_info = {}
script_info['brief_description'] = "Identify the core microbiome."
script_info['script_description'] = ""
script_info['script_usage'] = [("","Identify the core OTUs in in1.biom, defined as the OTUs that are present in at least 50% of the samples. Write the list of core OTUs to a text file, and a new BIOM file containing only the core OTUs.","%prog -i otu_table.biom -o core_50.txt --fraction_for_core 0.5 --output_table_fp otu_table_core_50.biom")]
script_info['output_description']= ""
script_info['required_options'] = [\
 make_option('-i','--input_fp',type="existing_filepath",help='the input otu table in BIOM format'),
 make_option('-o','--output_fp',type="new_filepath",help='the output core microbiome summary'),
 make_option('--fraction_for_core',type="float",
             help='the fraction of samples that an OTU must be observed in to be considered part of the core as a number in the range [0,1] [default: %default]'),
]
script_info['optional_options'] = [
 make_option('--output_table_fp',type="new_filepath",
             help='the otu table filtered to core otus only [default: %default]'),\
 make_option('--otu_md',default='taxonomy',
             help='the otu metadata category to write to the output file [defualt: %default]'),\
]
script_info['version'] = __version__



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    input_fp = opts.input_fp
    output_fp = opts.output_fp
    output_table_fp = opts.output_table_fp
    fraction_for_core = opts.fraction_for_core
    otu_md = opts.otu_md
    
    # get core across all samples
    sample_ids = None
    
    input_table = parse_biom_table(open(input_fp,'U'))
    output_f = open(output_fp,'w')
    
    try:
        core_table = filter_table_to_core(input_table,
                                          sample_ids,
                                          fraction_for_core)
    except TableException:
        output_f.write("# No OTUs present in %1.2f%% of samples." % (fraction_for_core * 100.))
        output_f.close()
        exit()
    
    # write the otu id and corresponding metadata for all core otus
    for value, id_, md in core_table.iterObservations():
        output_f.write('%s\t%s\n' % (id_,md[otu_md]))
    output_f.close()
    
    # write the core biom table, if requested by the user
    if output_table_fp:
        output_table_f = open(output_table_fp,'w')
        output_table_f.write(format_biom_table(core_table))
        output_table_f.close()
        
        
        
        
        


if __name__ == "__main__":
    main()