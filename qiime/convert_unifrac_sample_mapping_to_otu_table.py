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

from qiime.parse import parse_sample_mapping
from biom.table import table_factory

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
    """
    data = []
    sample_ids = []
    observation_ids = []
    for line in lines:
        fields = line.strip().split()
        observation_id = fields[0]
        sample_id = fields[1]
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