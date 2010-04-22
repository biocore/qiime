#!/usr/bin/env python
"""Application controller for uclust version 1.1.579

Includes application controllers for uclust and
convenience wrappers for different functions of uclust including
sorting fasta files, finding clusters, converting to cd-hit format and
searching and aligning against a database. Also contains
a parser for the resulting .clstr file.

Modified from cogent.app.cd_hit.py on 1-21-10, written by Daniel McDonald.
"""

__author__ = "William Walters"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["William Walters","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0.dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Development"

from os import remove, makedirs
from os.path import split, splitext, basename, isdir, abspath, isfile
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath,\
 get_tmp_filename, ApplicationError, ApplicationNotFoundError
from cogent.util.misc import reverse_complement, remove_files

class UclustParseError(Exception):
    pass

class Uclust(CommandLineApplication):
    """ Uclust ApplicationController
    
    """
    
    _command = 'uclust'
    _input_handler = '_input_as_parameters'
    _parameters = {\
        
        # Fasta input file for merge-sort function
        '--mergesort':ValuedParameter('--',Name='mergesort',Delimiter=' ',
            IsPath=True),
    
        # Output file, used by several difference functions
        '--output':ValuedParameter('--',Name='output',Delimiter=' ',
            IsPath=True),
    
        # Sets temp directory for uclust to create temp fasta file
        '--tmpdir':ValuedParameter('--',Name='tmpdir',Delimiter=' ',
            IsPath=True),
    
        # input filename, fasta format
        '--input':ValuedParameter('--',Name='input',Delimiter=' ',
            IsPath=True),
        
        # Output filename will be in uclust (.uc) format
        # Output cluster file, required parameter
        '--uc':ValuedParameter('--',Name='uc',Delimiter=' ',
            IsPath=True),
        
        # ID percent for OTU, by default is 97%
        '--id':ValuedParameter('--',Name='id',Delimiter=' ',IsPath=False),

        # Disable reverse comparison option, if norev is disabled
        # memory usage is expected to double for uclust
        '--rev':FlagParameter('--',Name='rev'),

        # 'library' file -- a reference of sequences representing pre-existing
        # clusters
        '--lib':ValuedParameter('--',Name='lib',Delimiter=' ',IsPath=True),
        
        # only compare sequences to the library file, don't add new clusters
        # for sequences which don't hit the library
        '--libonly':FlagParameter('--',Name='libonly'),
        
        # the max number of matches to review when looking for the best match
        '--maxaccepts':ValuedParameter('--',Name='maxaccepts',Delimiter=' '),
        
        # 
        '--maxrejects':ValuedParameter('--',Name='maxrejects',Delimiter=' '),
        
        # output fp for pairwise aligned sequences
        '--fastapairs':ValuedParameter('--',Name='fastapairs',Delimiter=' ',
            IsPath=True),
        
        # input filename, .uc format
        '--uc2clstr':ValuedParameter('--', Name='uc2clstr', Delimiter=' ',
            IsPath=True),

        # Don't assume input is sorted by length (default assume sorted).
        '--usersort':FlagParameter('--',Name='usersort'),
        
        # Same as --maxrejects 0 --maxaccepts 0 --nowordcountreject -- 
        # comes with a performance hit.
        '--optimal':FlagParameter('--',Name='optimal'),
        
    }
     
    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self,data):
        """ Set the input path (a fasta filepath)
        """
        # The list of values which can be passed on a per-run basis
        allowed_values = ['--input','--uc','--fastapairs','--lib',\
                           '--uc2clstr','--output','--mergesort']
        
        unsupported_parameters = set(data.keys()) - set(allowed_values)
        if unsupported_parameters:
            raise ApplicationError,\
             "Unsupported parameter(s) passed when calling uclust: %s" %\
              ' '.join(unsupported_parameters)
        
        for v in allowed_values:
            # turn the parameter off so subsequent runs are not
            # affected by parameter settings from previous runs
            self.Parameters[v].off()
            if v in data:
                # turn the parameter on if specified by the user
                self.Parameters[v].on(data[v])
        
        return ''
        
    def _get_result_paths(self,data):
        """ Set the result paths """
        
        result = {}
        
        result['Output'] = ResultPath(\
         Path=self.Parameters['--output'].Value,\
         IsWritten=self.Parameters['--output'].isOn())
         
        result['ClusterFile'] = ResultPath(
         Path = self.Parameters['--uc'].Value,
         IsWritten=self.Parameters['--uc'].isOn())
         
        result['PairwiseAlignments'] = ResultPath(
         Path = self.Parameters['--fastapairs'].Value,
         IsWritten=self.Parameters['--fastapairs'].isOn())
         
        return result
        
    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
        """
        UCLUST is hosted as an open source project at:
        http://www.drive5.com/uclust/

        The following papers should be cited if this resource is used:

        Paper pending. Check with Robert Edgar who is writing the paper
        for uclust as of March 2010.  Cite the above URL for the time being.
        """
        return help_str
         
def uclust_fasta_sort_from_filepath(
    fasta_filepath,
    output_filepath=None,
    HALT_EXEC=False):
    """Generates sorted fasta file via uclust --mergesort."""

    output_filepath = output_filepath or \
     get_tmp_filename(prefix='uclust_fasta_sort', suffix='.fasta')
    tmp_working_dir = split(output_filepath)[0]
    
    app = Uclust(params={'--tmpdir':tmp_working_dir},HALT_EXEC=HALT_EXEC)
    
    app_result = app(data={'--mergesort':fasta_filepath,\
                           '--output':output_filepath})
    
    return app_result


def get_next_hit_record(lines):
    for line in lines:
        line = line.strip()
        if line.startswith('H'):
            yield line.split('\t')
    return
            
def get_next_two_fasta_records(lines):
    result = []
    for record in MinimalFastaParser(lines):
        result.append(record)
        if len(result) == 2:
            yield result
            result = []
    return

def process_uclust_pw_alignment_results(fasta_pairs_lines,uc_lines):
    """ Process results of uclust search and align """
    alignments = get_next_two_fasta_records(fasta_pairs_lines)
    for hit in get_next_hit_record(uc_lines):
        matching_strand = hit[4]
        if matching_strand == '-':
            strand_id = '-'
            target_rev_match = True
        elif matching_strand == '+':
            strand_id = '+'
            target_rev_match = False
        elif matching_strand == '.':
            # protein sequence, so no strand information
            strand_id = ''
            target_rev_match = False
        else:
            raise UclustParseError, "Unknown strand type: %s" % matching_strand
        uc_query_id = hit[8]
        uc_target_id = hit[9]
        percent_id = float(hit[3])
        
        fasta_pair = alignments.next()
        
        fasta_query_id = fasta_pair[0][0]
        aligned_query = fasta_pair[0][1]
        
        if fasta_query_id != uc_query_id:
            raise UclustParseError,\
             "Order of fasta and uc files do not match."+\
             " Got query %s but expected %s." %\
              (fasta_query_id, uc_query_id)
            
        fasta_target_id = fasta_pair[1][0]
        aligned_target = fasta_pair[1][1]
            
        if fasta_target_id != uc_target_id + strand_id:
            raise UclustParseError, \
             "Order of fasta and uc files do not match."+\
             " Got target %s but expected %s." %\
              (fasta_target_id, uc_target_id + strand_id)
            
        if target_rev_match:
            query_id = uc_query_id + ' RC'
            aligned_query = reverse_complement(aligned_query)
            target_id = uc_target_id
            aligned_target = reverse_complement(aligned_target)
        else:
            query_id = uc_query_id
            aligned_query = aligned_query
            target_id = uc_target_id
            aligned_target = aligned_target
            
        yield (query_id, target_id, aligned_query, aligned_target,percent_id)

def uclust_search_and_align_from_fasta_filepath(
    query_fasta_filepath,
    subject_fasta_filepath,
    percent_ID=0.75,
    enable_rev_strand_matching=True,
    max_accepts=8,
    max_rejects=32,
    HALT_EXEC=False):
    """ query seqs against subject fasta using uclust
    
       return global pw alignment of best match
    """
     
    # Explanation of parameter settings
    #  id - min percent id to count a match
    #  maxaccepts = 8 , searches for best match rather than first match 
    #                   (0 => infinite accepts, or good matches before 
    #                    quitting search)
    #  maxaccepts = 32, 
    #  libonly = True , does not add sequences to the library if they don't
    #                   match something there already. this effectively makes
    #                   uclust a search tool rather than a clustering tool
    
    params = {'--id':percent_ID,\
              '--maxaccepts':max_accepts,\
              '--maxrejects':max_rejects,\
              '--libonly':True}
              
    if enable_rev_strand_matching:
        params['--rev'] = True
    
    # instantiate the application controller
    app = Uclust(params,HALT_EXEC=HALT_EXEC)
    
    # apply uclust
    alignment_filepath = \
     get_tmp_filename(prefix='uclust_alignments',suffix='.fasta')
    uc_filepath = \
     get_tmp_filename(prefix='uclust_results',suffix='.uc')
    input_data = {'--input':query_fasta_filepath,\
                  '--fastapairs':alignment_filepath,
                  '--lib':subject_fasta_filepath,\
                  '--uc':uc_filepath}
    app_result = app(input_data)
    
    # yield the pairwise alignments
    for result in process_uclust_pw_alignment_results(
     app_result['PairwiseAlignments'],app_result['ClusterFile']):
        try:
            yield result
        except GeneratorExit:
            break
    
    # clean up the temp files that were generated
    app_result.cleanUp()
    
    return

def uclust_cluster_from_sorted_fasta_filepath(
    fasta_filepath,
    output_filepath=None, 
    percent_ID=0.97, 
    optimal = False,
    suppress_sort = False,
    enable_rev_strand_matching=False,
    HALT_EXEC=False):
    """ Returns clustered uclust file from sorted fasta"""
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='uclust_clusters',suffix='.uc')
    
    params = {'--id':percent_ID}
    app = Uclust(params,HALT_EXEC=HALT_EXEC)
    
    # Set any additional parameters specified by the user
    if enable_rev_strand_matching: app.Parameters['--rev'].on()
    if optimal: app.Parameters['--optimal'].on()
    if suppress_sort: app.Parameters['--usersort'].on()
    
    app_result = app({'--input':fasta_filepath,'--uc':output_filepath})
    return app_result

        
def uclust_convert_uc_to_cdhit_from_filepath(
    uc_filepath,
    output_filepath=None,
    HALT_EXEC=False):
    """ Returns cdhit (.clstr) file from input uclust (.uc) file"""
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='uclust_to_cdhit',suffix='.clstr') 
    # prefix for tmp files should indicate what created it
    app = Uclust(HALT_EXEC=HALT_EXEC)
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
    

def get_clusters_from_fasta_filepath(
    fasta_filepath,
    percent_ID=0.97,
    optimal=False,
    suppress_sort=False,
    output_dir=None,
    enable_rev_strand_matching=False):
    """ Main convenience wrapper for using uclust to generate cluster files
    
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
        
    sorted_fasta_filepath = ""
    uc_filepath = ""
    clstr_filepath = ""


    # Error check in case any app controller fails
    files_to_remove = []
    try:
        if not suppress_sort:
            # Sort fasta input file from largest to smallest sequence 
            sort_fasta = uclust_fasta_sort_from_filepath(fasta_filepath, \
            fasta_output_filepath)
            # Get sorted fasta name from application wrapper
            sorted_fasta_filepath = sort_fasta['Output'].name
            files_to_remove.append(sorted_fasta_filepath)
        else:
            sort_fasta = None
            sorted_fasta_filepath = fasta_filepath
    
        # Generate uclust cluster file (.uc format)
        uclust_cluster = \
         uclust_cluster_from_sorted_fasta_filepath(sorted_fasta_filepath,
         uc_output_filepath, percent_ID = percent_ID,
         optimal = optimal, suppress_sort = suppress_sort,
         enable_rev_strand_matching = enable_rev_strand_matching)
        # Get cluster file name from application wrapper
        uc_filepath = uclust_cluster['ClusterFile'].name
        files_to_remove.append(uc_filepath)

        # Convert the .uc file to a cdhit (.clstr) format
        cdhit_conversion = \
         uclust_convert_uc_to_cdhit_from_filepath(uc_filepath, cd_hit_filepath)

        # Get the open file object from the cdhit conversion wrapper
        clstr_filepath = cdhit_conversion['Output']
        files_to_remove.append(clstr_filepath)
    except ApplicationError:
        remove_files(files_to_remove)
        raise ApplicationError, ('Error running uclust, make sure the proper '+\
         'version (1.1.579 or greater) is installed and the input fasta file '+\
         'is properly formatted.')
    except ApplicationNotFoundError:
        remove_files(files_to_remove)
        raise ApplicationNotFoundError('uclust not found, is it properly '+\
         'installed?')
    
    # Get list of lists for each cluster
    clusters = parse_uclust_clstr_file(clstr_filepath)
    
    # Remove temp files unless user specifies output filepath
    if not output_dir:
        try:
            sort_fasta.cleanUp()
        except AttributeError:
            pass
        uclust_cluster.cleanUp()
        cdhit_conversion.cleanUp()
    
    if len(clusters) == 0:
        raise ApplicationError, ('Clusters result empty, please check source '+\
         'fasta file for proper formatting.')
    
    return clusters

