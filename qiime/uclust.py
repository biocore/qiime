#!/usr/bin/env python
"""Application controller for uclust version 1.0.50

Includes application controllers for different functions of uclust
(sorting fasta files, finding clusters, converting to cd-hit format) and
a parser for the resulting .clstr file, along with 
convienence functions for calling these controllers.

Modified from cogent.app.cd_hit.py on 1-21-10, written by Daniel McDonald.
Currently only implementing DNA clustering, and percent ID for OTU as
parameters. """

__author__ = "William Walters"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.5.0.dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Development"

import shutil
from os import remove, makedirs
from os.path import split, splitext, basename, isdir, abspath
from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath,\
 get_tmp_filename

class UclustFastaSort(CommandLineApplication):
    """ ApplicationController for sorting a fasta file according to seq lens
    
    """
    
    _command = 'uclust'
    _parameters = {\
     # Fasta input file
     '--mergesort':ValuedParameter('--',Name='mergesort',Delimiter=' ',\
      IsPath=True),
     # Sorted fasta output file by length; fasta input file for uclust
     # needs to be arranged in order of largest to shortest seq lens.
     '--output':ValuedParameter('--',Name='output',Delimiter=' ',IsPath=True)
     }
     
    _input_handler = '_input_as_parameters'
    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self,data):
        """ Set the input path (a fasta filepath)
        """
        for param_id, param_value in data.items():
            self.Parameters[param_id].on(param_value)
        return ''
        
    def _get_result_paths(self,data):
        """ Set the result paths """
        
        result = {}
        result['SortedFasta'] = ResultPath(\
         Path=self.Parameters['--output'].Value,\
         IsWritten=True)
        return result
        
    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
        """
        UCLUST is hosted as an open source project at:
        http://www.drive5.com/uclust/

        The following papers should be cited if this resource is used:

        Paper pending-check with Robert Edgar who is writing the paper
        for uclust as of 1-21-2010.  Cite the above URL for the time being.
        """
        return help_str
         
def uclust_fasta_sort_from_filepath(fasta_filepath,output_filepath=None):
    """Generates sorted fasta file via uclust --mergesort."""

    # prefix for tmp files should indicate what created it
    app = UclustFastaSort()
    
    output_filepath = output_filepath or \
     get_tmp_filename(tmp_dir=app.WorkingDir, prefix='uclust_fasta_sort', \
     suffix='.fasta') 
    app_result = app(data={'--mergesort':fasta_filepath,\
                           '--output':output_filepath})
                           
    return app_result


class UclustCreateClusterFile(CommandLineApplication):
    """uclust Application Controller for creating a uclust cluster file

    Currently only using DNA clustering with uclust
    """

    _command = 'uclust'
    _input_handler = '_input_as_parameters'
    _parameters = {
        # input filename, fasta format, required parameter
        '--input':ValuedParameter('--',Name='input',Delimiter=' ',IsPath=True),
        
        # Output filename will be in uclust (.uc) format
        # Output cluster file, required parameter
        '--uc':ValuedParameter('--',Name='uc',Delimiter=' ',IsPath=True),
        
        # ID percent for OTU, by default is 97%
     	'--id':ValuedParameter('--',Name='id',Delimiter=' ',IsPath=False),
     	
     	# Disable reverse comparison option, if norev is disabled
     	# memory usage is expected to double for uclust
     	'--norev':FlagParameter('--',Name='norev'),
     	

    }
    
    
    _suppress_stdout = False
    _suppress_stderr = False
    
    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
        """
        UCLUST is hosted as an open source project at:
        http://www.drive5.com/uclust/

        The following papers should be cited if this resource is used:

        Paper pending-check with Robert Edgar who is writing the paper
        for uclust as of 1-21-2010.  Cite the above URL for the time being.
        """
        return help_str

    def _input_as_parameters(self,data):
        """ Set the input path (fasta) and output path (.uc), other parameters
        """
        for param_id, param_value in data.items():
          	if param_id == "--norev":
          		self.Parameters[param_id].on()
          	else:
          		self.Parameters[param_id].on(param_value)
        return ''
        
    def _get_result_paths(self,data):
        """ Set the result paths """
        
        result = {}
        result['ClusterFilepath'] = ResultPath(\
         Path = self.Parameters['--uc'].Value,\
         IsWritten=True)
        return result
        
   
def uclust_cluster_from_sorted_fasta_filepath(fasta_filepath, \
 output_filepath=None, percent_ID=0.97, enable_rev_strand_matching=False):
    """ Returns clustered uclust file from sorted fasta"""
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='uclust_clusters',suffix='.uc') 
    # prefix for tmp files should indicate what created it
    app = UclustCreateClusterFile()
    if enable_rev_strand_matching:
    	data={'--input':fasta_filepath,\
    	      '--uc':output_filepath,\
    	      '--id':percent_ID}
    else:
    	data={'--input':fasta_filepath,\
    	      '--uc':output_filepath,\
    	      '--id':percent_ID,\
    	      '--norev':True}
    app_result = app(data)
    return app_result


class UclustConvertToCdhit(CommandLineApplication):
    """uclust Application Controller for conversion to cdhit format

    Converts .uc files to .clstr (cd-hit format)
    """

    _command = 'uclust'
    _input_handler = '_input_as_parameters'
    _parameters = {
        # input filename, .uc format, required parameter
        '--uc2clstr':ValuedParameter('--', Name='uc2clstr', Delimiter=' ', \
         IsPath=True),
        
        # Output will be in cd-hit (.clstr) format
        
        # Output cluster file, required parameter
        '--output':ValuedParameter('--',Name='output',Delimiter=' ',IsPath=True)
        
    }
    
    
    _suppress_stdout = False
    _suppress_stderr = False
    
    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
        """
        UCLUST is hosted as an open source project at:
        http://www.drive5.com/uclust/

        The following papers should be cited if this resource is used:

        Paper pending-check with Robert Edgar who is writing the paper
        for uclust as of 1-21-2010.  Cite the above URL for the time being.
        """
        return help_str

    def _input_as_parameters(self,data):
        """ Set the input/output paths (.uc and .clstr files)
        """
        for param_id, param_value in data.items():
            self.Parameters[param_id].on(param_value)
        return ''
        
    def _get_result_paths(self,data):
        """ Set the result paths """
        
        result = {}
        result['CdhitFilepath'] = ResultPath(\
         Path = self.Parameters['--output'].Value,\
         IsWritten=True)
        return result
        
def uclust_convert_uc_to_cdhit_from_filepath(uc_filepath, \
 output_filepath=None):
    """ Returns cdhit (.clstr) file from input uclust (.uc) file"""
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='uclust_to_cdhit',suffix='.clstr') 
    # prefix for tmp files should indicate what created it
    app = UclustConvertToCdhit()
    app_result = app(data={'--uc2clstr':uc_filepath,\
                           '--output':output_filepath })
    return app_result
    
    
def parse_uclust_clstr_file(lines):
    """Returns a list of list of sequence ids representing clusters"""
    clusters = []
    curr_cluster = []

    for l in lines:
        if l.startswith('>Cluster'):
            if not curr_cluster:
                continue
            clusters.append(curr_cluster)
            curr_cluster = []
        else:
            # Slice off the '>' character
            seq_label = l.split()[2][1:]
            # Check for artifact "..." characters added by uclust
            # Remove if found
            if seq_label.endswith("..."):
            	seq_label = seq_label[:-3]
            curr_cluster.append(seq_label)

    if curr_cluster:
        clusters.append(curr_cluster)

    return clusters
    

    
def get_output_filepaths(output_dir, fasta_filepath):
    """ Returns filepaths for intermediate files to be kept """
    
    output_dir, output_filename = split(output_dir)
    output_dir = output_dir or './'
    output_file_basename, output_file_ext = splitext(fasta_filepath)
    fasta_output_filepath = '%s/%s_sorted.fasta' % \
     (output_dir,output_file_basename)
    uc_output_filepath = '%s/%s_sorted.uc' % \
     (output_dir,output_file_basename)
    cd_hit_filepath = '%s/%s_cdhit.clstr' % \
     (output_dir,output_file_basename)
     
    return fasta_output_filepath, uc_output_filepath, cd_hit_filepath, \
     output_dir
    
    

def get_clusters_from_fasta_filepath(fasta_filepath, percent_ID=0.97, \
 output_dir=None, enable_rev_strand_matching=False):
    """ Main convenience wrapper for uclust, returns list of lists of clusters.
    
    A source fasta file is required for the fasta_filepath.  This will be 
    sorted to be in order of longest to shortest length sequences.  Following
    this, the sorted fasta file is used to generate a cluster file in the
    uclust (.uc) format.  Next the .uc file is converted to cd-hit format
    (.clstr).  Finally this file is parsed and returned as a list of lists, 
    where each sublist a cluster of sequences.  If an output_dir is
    specified, the intermediate files will be preserved, otherwise all
    files created are temporary and will be deleted at the end of this 
    function
    
    The percent_ID parameter specifies the percent identity for a clusters,
    i.e., if 99% were the parameter, all sequences that were 99% identical
    would be grouped as a cluster.
    """
    
    # Create readable intermediate filenames if they are to be kept
    if output_dir:
        if not (output_dir.endswith("/")):
            output_dir += "/"
        fasta_output_filepath, uc_output_filepath, cd_hit_filepath, \
         output_dir = get_output_filepaths(output_dir, fasta_filepath)
        if not isdir(output_dir):
            makedirs(output_dir)
        
    else:
        fasta_output_filepath = None
        uc_output_filepath = None
        cd_hit_filepath = None
        

    # Sort fasta input file from largest to smallest sequence 
    sort_fasta = uclust_fasta_sort_from_filepath(fasta_filepath, \
     fasta_output_filepath)
    # Get sorted fasta name from application wrapper
    sorted_fasta_filepath = sort_fasta['SortedFasta'].name

    # Generate uclust cluster file (.uc format)
    uclust_cluster = \
     uclust_cluster_from_sorted_fasta_filepath(sorted_fasta_filepath, \
     uc_output_filepath, percent_ID = percent_ID, \
     enable_rev_strand_matching = False)
    # Get cluster file name from application wrapper
    uc_filepath = uclust_cluster['ClusterFilepath'].name

    # Convert the .uc file to a cdhit (.clstr) format
    cdhit_conversion = \
     uclust_convert_uc_to_cdhit_from_filepath(uc_filepath, cd_hit_filepath)

    # Get the open file object from the cdhit conversion wrapper
    clstr_filepath = cdhit_conversion['CdhitFilepath']
    
    # Get list of lists for each cluster
    clusters = parse_uclust_clstr_file(clstr_filepath)
    
    # Remove temp files unless user specifies output filepath
    if not output_dir:
        sort_fasta.cleanUp()
        uclust_cluster.cleanUp()
        cdhit_conversion.cleanUp()
    
    return clusters

