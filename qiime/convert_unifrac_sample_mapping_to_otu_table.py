#! /usr/bin/env python

__author__ = "Cathy Lozupone"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Catherine Lozupone", "Greg Caporaso",
                "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Cathy Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Development"

from biom.table import table_factory
from string import letters, digits, maketrans
from copy import deepcopy

def build_sample_ids_transtable():
    """Build translation table for sample ids being MIENS compliant"""
    all_chars = ''.join([chr(i) for i in range(128)])
    valid_sample_id_chars = letters + digits + "."
    non_valid_sample_id_chars = all_chars.translate(maketrans("",""),
        valid_sample_id_chars)
    trans_table = maketrans(non_valid_sample_id_chars,
        "."*len(non_valid_sample_id_chars))
    return trans_table

def parse_sample_mapping(lines):
    """Parses the UniFrac sample mapping file (environment file)

    The sample mapping file is a required input for the UniFrac web interface.
    Returns a dict of OTU names mapped to sample:count dictionaries.
    This code is used to convert this file to an OTU table for QIIME

    Corrects the sample ids to be MIENS compliant
    """
    #add the count of 1 if count info is not supplied
    new_lines = []
    for line in lines:
        line = line.strip().split('\t')
        if len(line) == 2:
            line.append('1')
        new_lines.append(line)

    trans_table = build_sample_ids_transtable()

    all_sample_names = [line[1].translate(trans_table) for line in new_lines]
    all_sample_names = set(all_sample_names)
    #create a dict of dicts with the OTU name mapped to a dictionary of
    #sample names with counts
    OTU_sample_info = {}
    for line in new_lines:
        OTU_name = line[0]
        if OTU_name not in OTU_sample_info:
            sample_info = dict([(i,'0') for i in all_sample_names])
            OTU_sample_info[OTU_name] = deepcopy(sample_info)
        sample_name = line[1].translate(trans_table)
        count = line[2]
        OTU_sample_info[OTU_name][sample_name] = count
    return OTU_sample_info, all_sample_names

def sample_mapping_to_otu_table(lines):
    """Converts the UniFrac sample mapping file to an OTU table
    
    The sample mapping file is a required input for the UniFrac web interface.
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

def sample_mapping_to_biom_table(lines):
    """Converts the UniFrac sample mapping file to biom table object
    
    The sample mapping file is a required input for the UniFrac web interface.

    Corrects the sample ids to be MIENS compliant
    """
    trans_table = build_sample_ids_transtable()

    data = []
    sample_ids = []
    observation_ids = []
    for line in lines:
        fields = line.strip().split()
        observation_id = fields[0]
        sample_id = fields[1].translate(trans_table)
        count = float(fields[2])
        
        try:
            sample_idx = sample_ids.index(sample_id)
        except ValueError:
            sample_idx = len(sample_ids)
            sample_ids.append(sample_id)
        try:
            observation_idx = observation_ids.index(observation_id)
        except ValueError:
            observation_idx = len(observation_ids)
            observation_ids.append(observation_id)
            
        data.append([observation_idx, sample_idx, count])
    
    return table_factory(data,sample_ids,observation_ids)