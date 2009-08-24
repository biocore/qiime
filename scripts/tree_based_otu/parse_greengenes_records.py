#!/usr/bin/env python

"""Parse the Greengenes formatted sequence data records for taxonomy info

The script is intended to be used with the following input:
http://greengenes.lbl.gov/Download/Sequence_Data/Greengenes_format/greengenes16SrRNAgenes.txt.gz
"""

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"

# from cogent.parse.fasta
def is_blank(x):
    """Checks if x is blank."""
    return (not x) or x.isspace()

def split_delim_f(delim):
    """Factory function for delim spliiter"""
    def inner_f(line):
        k,v = line.split(delim, 1)
        if v is '' and k == 'core_set_member':
            v = 0
        elif v is '':
            v = None
        return (k,v)
    return inner_f

def DemarkedParser(lines, delim, recordtype=dict, start='BEGIN', end='END'):
    """Yields successive recordtype objects from lines

    Records start with start and end with end.
    """
    is_start = lambda x: x == start
    is_end = lambda x: x == end
    split_line = split_delim_f(delim)

    for line in lines:
        line = line.strip()
        if is_blank(line):
            continue
        if is_end(line):
            yield curr_record
            curr_record = {}
            continue
        if is_start(line):
            curr_record = {}
            continue
        key, value = split_line(line)
        curr_record[key] = value

"""
prokMSA_id This has the id of the sequence in the tree that Daniel is using
ncbi_tax_id Taxon ID from NCBI, same as in the taxonomy file John has parsed
ncbi_acc_w_ver This has the NCBI accession number
core_set_member This has true/false according to whether the seq is in the "core set" of good sequences
ncbi_tax_string_format_2 This has the NCBI taxonomy string, semicolon-delimited (i.e. the successive semicolon-delimited fields are taxon names)
Hugenholtz_tax_string_format_2 ditto for Hugenholtz taxonomy
Ludwig_tax_string_format_2 ditto for Ludwig taxonomy
Pace_tax_string_format_2 ditto for Pace taxonomy
RDP_tax_string_format_2 ditto for RDP taxonomy
"""

taxonomy_fields = ['prokMSA_id','ncbi_tax_id','ncbi_acc_w_ver',
                   'core_set_member','ncbi_tax_string_format_2',
                   'Hugenholtz_tax_string_format_2',
                   'Ludwig_tax_string_format_2',
                   'Pace_tax_string_format_2',
                   'RDP_tax_string_format_2']

def taxonomy_delim_string(record, delim='\t'):
    """Returns a delimited string in taxonomy_fields order

    Will raise if field is not present. This is intended. Would like to know
    if the input data is bad
    """
    return delim.join(map(str, [record[k] for k in taxonomy_fields]))

def main(input_file, output_file):
    input = open(input_file)
    output = open(output_file, 'w')

    output.write('#' + '\t'.join(taxonomy_fields))

    write_buffer = []
    buffer_max = 1000
    for count, record in enumerate(DemarkedParser(input, '=')):
        write_buffer.append(taxonomy_delim_string(record))

        if count % buffer_max == 0:
            output.write('\n')
            output.write('\n'.join(write_buffer))
            write_buffer = []
    if write_buffer:
        output.write('\n')
        output.write('\n'.join(write_buffer))
    output.close()
    input.close()

if __name__ == '__main__':
    from sys import argv
    main(argv[1], argv[2])

