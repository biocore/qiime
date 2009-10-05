#!/usr/bin/env python
#file filter_otus_by_sample.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2009, the 454 Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Prototype"

"""
Author: Jesse Stombaugh (jesse.stombaugh@colorado.edu)
Status: Prototype

Requirements:
Python 2.5

Example 1: Filter otus and representative sequence set
Usage: filter_otus_by_sample.py [options] {-i input_otu_path -f fasta_file -s samples_to_extract}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:python filter_otus_by_sample.py -i otus.txt -f seqs.fna -s 'ControlA,ControlB' -o './Directory'
"""

from commands import getoutput
from string import strip
from parse import fields_to_dict
from numpy import array
from optparse import OptionParser
from time import strftime
import os
import re
from cogent import LoadSeqs
from make_3d_plots import create_dir

def filter_otus(aln,otus,prefs):
    """filters the otus file based on which samples should be removed and 
       determines which sequences to remove"""
    seqs_to_remove=[]
    new_otus_list=[]

    #create an list containing the seqs to remove
    for i in range(aln.getNumSeqs()):
        aln_name=aln.Names[i]
        is_sample=False
        otu_key=str(i)
        for sample_id in prefs:
           for j in range(len(otus[otu_key])):
               if re.match(prefs[sample_id],otus[otu_key][j]):
                   is_sample=True
                   break
           if is_sample:
               break
               
        #determine if sampleid is in otu and add to appropriate list
        if is_sample:
            seqs_to_remove.append(otu_key)
        else:
            new_otus_list.append((otu_key,otus[otu_key]))

    return seqs_to_remove,new_otus_list
    
def filter_aln_by_otus(new_otus_list,aln,seqs_to_remove):
    """filters the representative set of seqs based on which samples should
        be removed"""
    filtered_seqs=[]
    removed_seqs=[]
    for j in range(aln.getNumSeqs()):
        remove=False
        aln_name=aln.Names[j]
        for i in range(len(seqs_to_remove)):
            seq_identifier=seqs_to_remove[i]+' '
            if re.match(seq_identifier,aln_name):
                remove=True
                break

        #Create list of filtered and removed seqs
        if remove:
            removed_seqs.append((aln_name,aln.getSeq(aln_name)))
        else:
            filtered_seqs.append((aln_name,aln.getSeq(aln_name)))
    
    return filtered_seqs,removed_seqs

def process_extract_samples(samples_to_extract):
    """Parses the samples to extract option from the command line"""
    prefs = {}

    if samples_to_extract:
        samples = samples_to_extract.strip().strip("'").split(',')

    for j, col in enumerate(samples):
        key = str(j)   
        prefs[key]={}
        prefs[key]=col

    return prefs

def _do_sample_filter(prefs, data, dir_path='', filename=None):
    """processes the filtering of the otus file and representative seq set, then
        writes filtered otus and filtered representative seq set files""" 
        
    aln=data['aln']
    otus=data['otus']

    #filter the otus file based on which samples to remove
    seqs_to_remove, new_otus_list=filter_otus(aln,otus,prefs)

    filtered_otus_output_filepath = '%s/%s_sfiltered_otus.txt' \
                                    % (dir_path,filename)
    filtered_otus_output_filepath=open(filtered_otus_output_filepath,'w')
    
    # Write out a new otus file
    for key in range(len(new_otus_list)):
        filtered_otus_output_filepath.write(str(key)+'\t')
        for j in range(len(new_otus_list[key][1])):
            if j<len(new_otus_list[key][1])-1:
                filtered_otus_output_filepath.write(new_otus_list[key][1][j]+\
                                                    '\t')
            else:
                filtered_otus_output_filepath.write(new_otus_list[key][1][j])
        if key<len(new_otus_list):
            filtered_otus_output_filepath.write('\n')
    filtered_otus_output_filepath.close()

    #filter representative seq set
    filtered_seqs,removed_seqs=filter_aln_by_otus(new_otus_list,aln,\
                                                    seqs_to_remove)

    #write a fasta containing list of sequences removed from 
    #representative set
    removed_seqs=LoadSeqs(data=removed_seqs,aligned=False)
    output_filepath2 = '%s/%s_sremoved.fasta' % (dir_path,filename)
    output_file2=open(output_filepath2,'w')
    output_file2.write(removed_seqs.toFasta())
    output_file2.close()

    #write a fasta containing the filtered representative seqs
    filtered_seqs=LoadSeqs(data=filtered_seqs,aligned=False)
    output_filepath = '%s/%s_sfiltered.fasta' % (dir_path,filename)
    output_file=open(output_filepath,'w')
    output_file.write(filtered_seqs.toFasta())
    output_file.close()

def _make_cmd_parser():
    """Returns the command-line options"""
    parser = OptionParser(usage="Usage: this_file.py -i <otus file> \
-f <fasta file> -s <SampleIDs> -o <write to directory>")

    parser.add_option('-i', '--input_otu_path',\
        help='otu file path [REQUIRED]')
    parser.add_option('-f', '--fasta_file', \
        help='name of fasta file [REQUIRED]')
    parser.add_option('-s', '--samples_to_extract',\
        help='samples to extract [REQUIRED]')
    parser.add_option('-o', '--dir_path',\
        help='directory prefix for all analyses [default=%default]',default='')
    options, args = parser.parse_args()
    return options

def _process_prefs(options):
    """opens files as necessary based on prefs"""
    data = {}

    fasta_file = options.fasta_file

    # load the input alignment
    data['aln'] = LoadSeqs(fasta_file,aligned=False)

    #Load the otu file
    otu_path=options.input_otu_path
    otu_f = open(otu_path, 'U')
    otus = fields_to_dict(otu_f)
    otu_f.close()
    data['otus']=otus

    #Determine which which samples to extract from representative seqs
    #and from otus file
    if options.samples_to_extract:
        prefs=process_extract_samples(options.samples_to_extract)

    filepath=options.fasta_file
    filename=filepath.strip().split('/')[-1]
    filename=filename.split('.')[0]

    file_path=__file__.split('/')
    if len(file_path)==1:
        qiime_dir='./'
    else:
        qiime_dir='';
        for i in range(len(file_path)-1):
            qiime_dir+=file_path[i]+'/'

    dir_path = create_dir(options.dir_path,'filtered_by_otus')

    action_str = '_do_sample_filter'
    try:
        action = eval(action_str)
    except NameError:
        action = None
    #Place this outside try/except so we don't mask NameError in action
    if action:
        action(prefs, data, dir_path,filename)
    
if __name__ == "__main__":
    from sys import argv, exit
    options = _make_cmd_parser()
    
    #Kept, just in case we allow for reading a prefs file
    #prefs = eval(open(options.pref_fname, 'U').read()) 

    _process_prefs(options)