#! /usr/bin/env python

"""script to convert a sample mapping file to an OTU table. 

Allows for users that have already created sample mapping files for use with theUnifrac web interface to use QIIME, which records this information in an OTU table.

Usage: python convert_SampleMapping_to_OTUtable.py SampleMappingFile.txt OTUtable_outfile.txt 
"""

from copy import deepcopy
from sys import argv

def SampleMapping_To_OTUtable(lines):
    """Convert sample mapping file to OTU table
    """
    out = ["#Full OTU Counts"]
    header = ["#OTU ID"]
   
    OTU_sample_info, all_sample_names = parse_sample_mapping(lines)
    all_sample_names = list(all_sample_names)
    all_sample_names.sort()
    header.extend(all_sample_names)
    out.append('\t'.join(header))
    for OTU in OTU_sample_info:
        new_line = []
        new_line.append(OTU)
        for sample in all_sample_names:
            new_line.append(OTU_sample_info[OTU][sample])
        out.append('\t'.join(new_line))
    return out

def parse_sample_mapping(lines):
    """Stores the info from the sample mapping file

    creates a dict of OTU names mapped to sample:count dict
    """
    lines = [line.strip().split('\t') for line in lines]

    all_sample_names = [line[1] for line in lines]
    all_sample_names = set(all_sample_names)
    #create a dict of dicts with the OTU name mapped to a dictionary of
    #sample names with counts
    OTU_sample_info = {}
    for line in lines:
        OTU_name = line[0]
        if OTU_name not in OTU_sample_info:
            sample_info = dict([(i,'0') for i in all_sample_names])
            OTU_sample_info[OTU_name] = deepcopy(sample_info)
        sample_name = line[1]
        count = line[2]
        OTU_sample_info[OTU_name][sample_name] = count
    return OTU_sample_info, all_sample_names

if __name__ == '__main__':
    sample_mapping_file = open(argv[1], 'U') 
    of = open(argv[2], 'w')
    
    result = SampleMapping_To_OTUtable(sample_mapping_file)
    of.write('\n'.join(result))
    of.close()


