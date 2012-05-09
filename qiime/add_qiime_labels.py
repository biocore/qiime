#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Release"

from os.path import join, basename
from glob import glob
from string import letters, digits

from cogent.parse.fasta import MinimalFastaParser

from qiime.util import duplicates_indices

def add_qiime_labels(mapping_f,
                     fasta_dir,
                     output_dir="./",
                     count_start=0):
    """ Main function for combining fasta files, writing valid QIIME labels
    
    mapping_f:  open file object of the SampleID to fasta file name.
    fasta_dir:  Directory of fasta files to combine into a single file.
    output_dir:  Directory to write output combined file to.
    count_start:  Number to start enumeration of fasta labels with.
    
    """
    
    fasta_name_to_sample_id = check_mapping_data(mapping_f)
    
    fasta_files = get_fasta_fps(fasta_dir)
    
    check_fasta_fps(fasta_name_to_sample_id, fasta_files)
    
    write_combined_fasta(fasta_name_to_sample_id, fasta_files, output_dir,
     counter=count_start)
    
def check_mapping_data(mapping_data):
    """ Checks mapping data for MIMARKS SampleIDs, unique IDs, fasta file names
    
    Also returns a dict of fasta file name: SampleID
    
    mapping_data:  list of lines of data from mapping file
    """
    
    valid_mimarks = letters + digits + "."
    
    fasta_name_to_sample_id = {}

    fasta_names = []
    for line in mapping_data:
        curr_line = line.strip().split('\t')
        if not curr_line or line.startswith("#"):
            continue
        try:
            fasta_name_to_sample_id[basename(curr_line[1].strip())] =\
             curr_line[0]
        except IndexError:
            raise IndexError,("Found non-tab separated line in mapping "+\
             "data.  Offending line is: %s" % line)
        for curr_char in curr_line[0]:
            if curr_char not in valid_mimarks:
                raise ValueError,("Found invalid character in line: %s\n" %\
                 line + "SampleIDs must be alphanumeric and . characters "+\
                 "only")
        fasta_names.append(curr_line[1].strip())
       
    fasta_name_dups = duplicates_indices(fasta_names)
    if fasta_name_dups:
        raise ValueError,("Found duplicate fasta names: %s" %\
         "\t".join([fasta_name for fasta_name in fasta_name_dups.keys()]))
         
    return fasta_name_to_sample_id
        
def get_fasta_fps(fasta_dir):
    """ Returns list of fasta filepaths (only .fna, .fasta, and .fa files)
    
    fasta_dir:  Directory of fasta files to check
    """
    fasta_files = glob(fasta_dir + "/*.fna")
    fasta_files.extend(glob(fasta_dir + "/*.fa"))
    fasta_files.extend(glob(fasta_dir + "/*.fasta"))
    
    return fasta_files
    
def write_combined_fasta(fasta_name_to_sample_id,
                         fasta_files,
                         output_dir="./",
                         counter=0):
    """ Writes combined, enumerated fasta file
    
    fasta_name_to_sample_id:  dict of fasta file name to SampleID
    fasta_files: list of filepaths to iterate through
    output_dir:  output directory to write combined file to
    counter:  Starting number to enumerate sequences with
    """
    
    combined_file_out = open(join(output_dir + "/", "combined_seqs.fna"), "w")
    
    for curr_fasta in fasta_files:
        for label, seq in MinimalFastaParser(open(curr_fasta, "U")):
            combined_file_out.write(">%s_%d %s\n" %\
             (fasta_name_to_sample_id[basename(curr_fasta)], counter, label))
            combined_file_out.write("%s\n" % seq)
            counter += 1
            
def check_fasta_fps(fasta_name_to_sample_id,
                    fasta_files):
    """ Checks that all fasta files found match fasta names in mapping file
    
    fasta_name_to_sample_id:  dict of fasta file names to sample IDs
    fasta_files:  list of fasta filepaths
    """
    
    for fp in fasta_files:
        curr_filename = basename(fp)
        if curr_filename not in fasta_name_to_sample_id.keys():
            raise ValueError,("found %s, but %s not found in mapping file." %\
             (fp, curr_filename))
    
    return True