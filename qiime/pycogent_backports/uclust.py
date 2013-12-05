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
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["William Walters","Greg Caporaso","Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Production"

from os import remove, makedirs
from os.path import split, splitext, basename, isdir, abspath, isfile, join
from tempfile import gettempdir
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath,\
 get_tmp_filename, ApplicationError, ApplicationNotFoundError
from cogent.util.misc import remove_files
from cogent import DNA

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
        
        # Maximum hits before quitting search (default 1, 0=infinity).
        '--maxaccepts':ValuedParameter('--',Name='maxaccepts',Delimiter=' '),
        
        # Maximum rejects before quitting search (default 8, 0=infinity). 
        '--maxrejects':ValuedParameter('--',Name='maxrejects',Delimiter=' '),
        
        # Target nr. of common words (default 8, 0=don't step)
        '--stepwords':ValuedParameter('--',Name='stepwords',Delimiter=' '),
        
        # Word length for windex (default 5 aa.s, 8 nuc.s).
        '--w':ValuedParameter('--',Name='w',Delimiter=' '),
        
        # output fp for pairwise aligned sequences
        '--fastapairs':ValuedParameter('--',Name='fastapairs',Delimiter=' ',
            IsPath=True),
        
        # input filename, .uc format
        '--uc2clstr':ValuedParameter('--', Name='uc2clstr', Delimiter=' ',
            IsPath=True),

        # Don't assume input is sorted by length (default assume sorted).
        '--usersort':FlagParameter('--',Name='usersort'),
        
        # Same as --maxrejects 0 --nowordcountreject.
        # comes with a performance hit.
        '--exact':FlagParameter('--',Name='exact'),

        # Same as --maxrejects 0 --maxaccepts 0 --nowordcountreject -- 
        # comes with a performance hit.
        '--optimal':FlagParameter('--',Name='optimal'),

        '--stable_sort':FlagParameter('--',Name='stable_sort'),

        # From uclust help:
        # Write all accepts to .uc file (default top hit/no match only).
        '--allhits':FlagParameter('--',Name='allhits'),
    }
     
    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self,data):
        """ Set the input path (a fasta filepath)
        """
        # The list of values which can be passed on a per-run basis
        allowed_values = ['--input','--uc','--fastapairs',\
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
        
    def _accept_exit_status(self,exit_status):
        """ Test for acceptable exit status
        
            uclust can seg fault and still generate a parsable .uc file
            so we explicitly check the exit status
        
        """
        return exit_status == 0
        
    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
        """
        UCLUST is hosted at:
        http://www.drive5.com/uclust/

        The following papers should be cited if this resource is used:

        Paper pending. Check with Robert Edgar who is writing the paper
        for uclust as of March 2010.  Cite the above URL for the time being.
        """
        return help_str

## Start functions for processing uclust output files
def get_next_record_type(lines,types):
    for line in lines:
        line = line.strip()
        if line and line[0] in types:
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
    for hit in get_next_record_type(uc_lines,'H'):
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
            aligned_query = DNA.rc(aligned_query)
            target_id = uc_target_id
            aligned_target = DNA.rc(aligned_target)
        else:
            query_id = uc_query_id
            aligned_query = aligned_query
            target_id = uc_target_id
            aligned_target = aligned_target
            
        yield (query_id, target_id, aligned_query, aligned_target, percent_id)

def clusters_from_uc_file(uc_lines,
                          error_on_multiple_hits=True):
    """ Given an open .uc file, return lists (clusters, failures, new_seeds)
    
        uc_lines: open .uc file, or similar object -- this is the output
         generated by uclust's -uc parameter
        error_on_multiple_hits: if True (default), when a single query hits
         to multiple seeds, as can happen when --allhits is passed to uclust,
         throw a UclustParseError. if False, when a single query hits to
         multiple seeds, it will appear in each cluster.
         
        This function processes all hit (H), seed (S), and no hit (N) lines
         to return all clusters, failures, and new_seeds generated in
         a uclust run. failures should only arise when users have passed
         --lib and --libonly, and a sequence doesn't cluster to any existing
         reference database sequences.
    
    """
    clusters = {}
    failures = []
    seeds = []
    all_hits = set()
    # the types of hit lines we're interested in here
    # are hit (H), seed (S), library seed (L) and no hit (N) 
    hit_types={}.fromkeys(list('HSNL'))
    for record in get_next_record_type(uc_lines,hit_types):
        hit_type = record[0]
        # sequence identifiers from the fasta header lines only 
        # (no comment data) are stored to identify a sequence in 
        # a cluster -- strip off any comments here as this value
        # is used in several places
        query_id = record[8].split()[0]
        target_cluster = record[9].split()[0]
        if hit_type == 'H':
            if error_on_multiple_hits and query_id in all_hits:
                raise UclustParseError, \
                 ("Query id " + query_id + " hit multiple seeds. "
                  "This can happen if --allhits is "
                  "enabled in the call to uclust, which isn't supported by default. "
                  "Call clusters_from_uc_file(lines, error_on_multiple_hits=False) to "
                  "allow a query to cluster to multiple seeds.")
            else:
                # add the hit to its existing cluster (either library
                # or new cluster)
                clusters[target_cluster].append(query_id)
                all_hits.add(query_id)
        elif hit_type == 'S':
            # a new seed was identified -- create a cluster with this 
            # sequence as the first instance
            if query_id in clusters:
                raise UclustParseError,\
                 ("A seq id was provided as a seed, but that seq id already "
                  "represents a cluster. Are there overlapping seq ids in your "
                  "reference and input files or repeated seq ids in either? "
                  "Offending seq id is %s" % query_id)
            clusters[query_id] = [query_id]
            seeds.append(query_id)
        elif hit_type == 'L':
            # a library seed was identified -- create a cluster with this 
            # id as the index, but don't give it any instances yet bc the hit
            # line will be specified separately. note we need to handle these
            # lines separately from the H lines to detect overlapping seq ids 
            # between the reference and the input fasta files
            if query_id in clusters:
                raise UclustParseError,\
                 ("A seq id was provided as a seed, but that seq id already "
                  "represents a cluster. Are there overlapping seq ids in your "
                  "reference and input files or repeated seq ids in either? "
                  "Offending seq id is %s" % query_id)
            clusters[query_id] = []
        elif hit_type == 'N':
            # a failure was identified -- add it to the failures list
            failures.append(query_id)
        else:
            # shouldn't be possible to get here, but provided for 
            # clarity
            raise UclustParseError,\
             "Unexpected result parsing line:\n%s" % '\t'.join(record)
    
    # will need to return the full clusters dict, I think, to support
    # useful identifiers in reference database clustering
    #return  clusters.values(), failures, seeds
    return  clusters, failures, seeds

## End functions for processing uclust output files


## Start uclust convenience functions
def uclust_fasta_sort_from_filepath(
    fasta_filepath,
    output_filepath=None,
    tmp_dir=gettempdir(),
    HALT_EXEC=False):
    """Generates sorted fasta file via uclust --mergesort."""
    output_filepath = output_filepath or \
     get_tmp_filename(tmp_dir=tmp_dir,prefix='uclust_fasta_sort',
                      suffix='.fasta')
    
    app = Uclust(params={'--tmpdir':tmp_dir},
                 TmpDir=tmp_dir,HALT_EXEC=HALT_EXEC)
    
    app_result = app(data={'--mergesort':fasta_filepath,\
                           '--output':output_filepath})
    
    return app_result

def uclust_search_and_align_from_fasta_filepath(
    query_fasta_filepath,
    subject_fasta_filepath,
    percent_ID=0.75,
    enable_rev_strand_matching=True,
    max_accepts=8,
    max_rejects=32,
    tmp_dir=gettempdir(),
    HALT_EXEC=False):
    """ query seqs against subject fasta using uclust, 
    
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
    
    params = {'--id':percent_ID,
              '--maxaccepts':max_accepts,
              '--maxrejects':max_rejects,
              '--libonly':True,
              '--lib':subject_fasta_filepath,
              '--tmpdir':tmp_dir}
              
    if enable_rev_strand_matching:
        params['--rev'] = True
    
    # instantiate the application controller
    app = Uclust(params,
                 TmpDir=tmp_dir,HALT_EXEC=HALT_EXEC)
    
    # apply uclust
    alignment_filepath = \
     get_tmp_filename(tmp_dir=tmp_dir,prefix='uclust_alignments',
                      suffix='.fasta')
    uc_filepath = \
     get_tmp_filename(tmp_dir=tmp_dir,prefix='uclust_results',
                      suffix='.uc')
    input_data = {'--input':query_fasta_filepath,
                  '--fastapairs':alignment_filepath,
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
    uc_save_filepath=None, 
    percent_ID=0.97, 
    max_accepts=1,
    max_rejects=8, 
    stepwords=8,
    word_length=8,
    optimal = False,
    exact = False,
    suppress_sort = False,
    enable_rev_strand_matching=False,
    subject_fasta_filepath=None,
    suppress_new_clusters=False,
    stable_sort=False,
    tmp_dir=gettempdir(),
    HALT_EXEC=False):
    """ Returns clustered uclust file from sorted fasta"""
    output_filepath = uc_save_filepath or \
     get_tmp_filename(tmp_dir=tmp_dir,prefix='uclust_clusters',
                      suffix='.uc')
     
    
    params = {'--id':percent_ID,
              '--maxaccepts':max_accepts,
              '--maxrejects':max_rejects,
              '--stepwords':stepwords,
              '--w':word_length,
              '--tmpdir':tmp_dir}
    app = Uclust(params,
                 TmpDir=tmp_dir,HALT_EXEC=HALT_EXEC)
    
    # Set any additional parameters specified by the user
    if enable_rev_strand_matching: app.Parameters['--rev'].on()
    if optimal: app.Parameters['--optimal'].on()
    if exact: app.Parameters['--exact'].on()
    if suppress_sort: app.Parameters['--usersort'].on()
    if subject_fasta_filepath: app.Parameters['--lib'].on(subject_fasta_filepath)
    if suppress_new_clusters: app.Parameters['--libonly'].on()
    if stable_sort: app.Parameters['--stable_sort'].on()
    
    app_result = app({'--input':fasta_filepath,'--uc':output_filepath})
    return app_result

def get_output_filepaths(output_dir, fasta_filepath):
    """ Returns filepaths for intermediate file to be kept """
    return join(output_dir,
                splitext(basename(fasta_filepath))[0] + '_clusters.uc')

def get_clusters_from_fasta_filepath(
    fasta_filepath,
    original_fasta_path,
    percent_ID=0.97,
    max_accepts=1,
    max_rejects=8, 
    stepwords=8,
    word_length=8,
    optimal=False,
    exact=False,
    suppress_sort=False,
    output_dir=None,
    enable_rev_strand_matching=False,
    subject_fasta_filepath=None,
    suppress_new_clusters=False,
    return_cluster_maps=False,
    stable_sort=False,
    save_uc_files=True,
    HALT_EXEC=False):
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
    
    fasta_output_filepath = None
    uc_output_filepath = None
    cd_hit_filepath = None
    
    if output_dir and not output_dir.endswith('/'):
        output_dir += '/'
    
    
    if save_uc_files:
        uc_save_filepath = get_output_filepaths(output_dir, original_fasta_path)
    else:
        uc_save_filepath = None
    
    sorted_fasta_filepath = ""
    uc_filepath = ""
    clstr_filepath = ""


    # Error check in case any app controller fails
    files_to_remove = []
    try:
        if not suppress_sort:
            # Sort fasta input file from largest to smallest sequence 
            sort_fasta = uclust_fasta_sort_from_filepath(fasta_filepath, \
            output_filepath=fasta_output_filepath)
            
            # Get sorted fasta name from application wrapper
            sorted_fasta_filepath = sort_fasta['Output'].name
            files_to_remove.append(sorted_fasta_filepath)
            
        else:
            sort_fasta = None
            sorted_fasta_filepath = fasta_filepath
        
        # Generate uclust cluster file (.uc format)
        uclust_cluster = uclust_cluster_from_sorted_fasta_filepath(
         sorted_fasta_filepath,
         uc_save_filepath, 
         percent_ID=percent_ID,
         max_accepts=max_accepts,
         max_rejects=max_rejects, 
         stepwords=stepwords,
         word_length=word_length,
         optimal=optimal, 
         exact=exact, 
         suppress_sort=suppress_sort,
         enable_rev_strand_matching=enable_rev_strand_matching,
         subject_fasta_filepath=subject_fasta_filepath,
         suppress_new_clusters=suppress_new_clusters,
         stable_sort=stable_sort,
         HALT_EXEC=HALT_EXEC)
        # Get cluster file name from application wrapper
        remove_files(files_to_remove)
    except ApplicationError:
        remove_files(files_to_remove)
        raise ApplicationError, ('Error running uclust. Possible causes are '
         'unsupported version (current supported version is v1.2.22) is installed or '
         'improperly formatted input file was provided')
    except ApplicationNotFoundError:
        remove_files(files_to_remove)
        raise ApplicationNotFoundError('uclust not found, is it properly '+\
         'installed?')

    # Get list of lists for each cluster
    clusters, failures, seeds = \
     clusters_from_uc_file(uclust_cluster['ClusterFile'])
    
    # Remove temp files unless user specifies output filepath
    if not save_uc_files:
        uclust_cluster.cleanUp()
    
    if return_cluster_maps:
        return clusters, failures, seeds
    else:
        return clusters.values(), failures, seeds

## End uclust convenience functions
