#!/usr/bin/env python
"""Application controller for usearch v4.2.66

Includes application controllers for usearch and
convenience wrappers for different functions of uclust including
sorting fasta files, finding clusters, converting to cd-hit format and
searching and aligning against a database. Also contains
a parser for the resulting .clstr file.

Modified from pycogent_backports/uclust.py, written by 
Greg Caporaso/William Walters
"""

__author__ = "William Walters"
__copyright__ = "Copyright 2007-2011, The PyCogent Project"
__credits__ = ["William Walters","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.0.dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Development"

from os import remove, makedirs
from os.path import split, splitext, basename, isdir, abspath, isfile
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath,\
 get_tmp_filename, ApplicationError, ApplicationNotFoundError
from cogent.util.misc import remove_files
from cogent import DNA

class UsearchParseError(Exception):
    pass

class Usearch(CommandLineApplication):
    """ Usearch ApplicationController
    
    """
    
    _command = 'usearch'
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
        
        # log filepath
        '--log':ValuedParameter('--', Name='log', Delimiter=' ', IsPath=True),
        
        # cluster command
        '--cluster':ValuedParameter('--', Name='cluster', Delimiter=' ',
         IsPath=True),
        
        
        # Size of compressed index table. Should be prime, e.g. 40000003.
        '--slots':ValuedParameter('--', Name='slots', Delimiter=' ',
         IsPath=False),
         
        # Not specified in usearch helpstring...
        '--sizein':FlagParameter('--', Name='sizein'),
         
        # Not specified in usearch helpstring...
        '--sizeout':FlagParameter('--', Name='sizeout'),
         
        # Not specified in usearch helpstring...
        '--minlen':ValuedParameter('--', Name='minlen', Delimiter=' ',
         IsPath=False),
         
        # output filepath for dereplicated fasta file
        '--seedsout':ValuedParameter('--', Name='seedsout', Delimiter=' ',
         IsPath=True),
         
        # Dereplicate exact subsequences
        '--derep_subseq':FlagParameter('--',Name='derep_subseq'),
        
        # Sort by abundance
        '--sortsize':ValuedParameter('--', Name='sortsize', Delimiter=' ',
         IsPath=True),
         
        # usearch search plus clustering
        '--consout':ValuedParameter('--', Name='consout', Delimiter=' ',
         IsPath=True),
         
        # Abundance skew setting for uchime de novo chimera detection
        '--abskew':ValuedParameter('--', Name='abskew', Delimiter=' ',
         IsPath=False),
         
        # input fasta filepath for uchime chimera
        '--uchime':ValuedParameter('--', Name='uchime', Delimiter=' ',
         IsPath=True),
        
        # output chimera filepath
        '--chimeras':ValuedParameter('--', Name='chimeras', Delimiter=' ',
         IsPath=True),
         
        # output non-chimera filepath
        '--nonchimeras':ValuedParameter('--', Name='nonchimeras',
         Delimiter=' ', IsPath=True),
         
        # reference sequence database for ref based chimera detection
        '--db':ValuedParameter('--', Name='db', Delimiter=' ', IsPath=True),
        
        # output clusters filepath for chimera detection
        '--uchimeout':ValuedParameter('--', Name='uchimeout', Delimiter=' ',
         IsPath=True),
         
        # minimum cluster size for quality filtering
        '--minsize':ValuedParameter('--', Name='minsize', Delimiter=' ',
         IsPath=False),
         
        # input fasta for blast alignments
        '--query':ValuedParameter('--', Name='query', Delimiter=' ',
         IsPath=True),
         
        # global alignment flag
        '--global':FlagParameter('--', Name='global')
         
        
    }
    
     
    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self,data):
        """ Set the input path (a fasta filepath)
        """
        # The list of values which can be passed on a per-run basis
        allowed_values = ['--input','--uc','--fastapairs',\
                           '--uc2clstr','--output','--mergesort', '--log',\
                           '--cluster', '--seedsout', '--sortsize',\
                           '--consout', '--uchime', '--chimeras',\
                           '--nonchimeras', '--db', '--uchimeout',\
                           '--query']
                           
        unsupported_parameters = set(data.keys()) - set(allowed_values)
        if unsupported_parameters:
            raise ApplicationError,\
             "Unsupported parameter(s) passed when calling usearch: %s" %\
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
        
            usearch can seg fault and still generate a parsable .uc file
            so we explicitly check the exit status
        
        """
        return exit_status == 0
        
    def getHelp(self):
        """Method that points to documentation"""
        help_str =\
        """
        USEARCH is hosted at:
        http://www.drive5.com/usearch/

        The following papers should be cited if this resource is used:

        Paper pending. Check with Robert Edgar who is writing the paper
        for usearch as of Aug. 2011
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


def clusters_from_blast_uc_file(uc_lines):
    """ Parses out hit/miss sequences from usearch blast uc file
    
    All lines should be 'H'it or 'N'o hit.  Returns a dict of OTU ids: sequence
    labels of the hits, and a list of all sequence labels that miss.
    
    uc_lines = open file object of uc file
    """
    
    hit_miss_index = 0
    cluster_id_index = 1
    seq_label_index = 8
    
    otus = {}
    unassigned_seqs = []
    
    for line in uc_lines:
        # skip empty, comment lines
        if line.startswith('#') or len(line.strip()) == 0:
            continue
        
        curr_line = line.split('\t')
        
        if curr_line[hit_miss_index] == 'N':
            # only retaining actual sequence label
            unassigned_seqs.append(curr_line[seq_label_index].split()[0])
        
        if curr_line[hit_miss_index] == 'H':
            
            curr_seq_label = curr_line[seq_label_index].split()[0]
            curr_otu_id = curr_line[cluster_id_index]
            # Append sequence label to dictionary, or create key
            try:
                otus[curr_otu_id].append(curr_seq_label)
            except KeyError:
                otus[curr_otu_id] = [curr_seq_label]
                
    return otus, unassigned_seqs
        
                
    

## End functions for processing uclust output files


## Start usearch convenience functions
def usearch_fasta_sort_from_filepath(
    fasta_filepath,
    output_filepath=None,
    log_name = "sortlen.log",
    HALT_EXEC=False,
    save_intermediate_files=False):
    """Generates sorted fasta file via usearch --mergesort.
    
    fasta_filepath: filepath to input fasta file
    output_filepath: filepath for output sorted fasta file.
    log_name: string to specify log filename
    HALT_EXEC: Used for debugging app controller
    save_intermediate_files: Preserve all intermediate files created."""
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='usearch_fasta_sort', suffix='.fasta')
    
    # using abspath to create log filepath, for cluster environments
    tmp_working_dir = split(abspath(output_filepath))[0]
    
    log_filepath = tmp_working_dir + "/" + log_name
    
    app = Usearch(params={},HALT_EXEC=HALT_EXEC)
    
    app_result = app(data={'--mergesort':fasta_filepath,\
                           '--output':output_filepath,\
                           '--log':log_filepath})
    
    return app_result, output_filepath
    
def usearch_dereplicate_exact_subseqs(
    fasta_filepath,
    output_filepath=None,
    minlen=64,
    w=64,
    slots=16769023,
    sizeout=True,
    maxrejects=64,
    log_name = "derep.log",
    usersort=False,
    HALT_EXEC=False,
    save_intermediate_files=False):
    """ Generates clusters and fasta file of dereplicated subsequences
    
    These parameters are those specified by Robert Edgar for optimal use of
    usearch in clustering/filtering sequences.
    
    fasta_filepath = input filepath of fasta file to be dereplicated
    output_filepath = output filepath of dereplicated fasta file
    minlen = (not specified in usearch helpstring)
    w = Word length for U-sorting
    slots = Size of compressed index table. Should be prime, e.g. 40000003.
     Should also specify --w, typical is --w 16 or --w 32.
    sizeout = (not specified in usearch helpstring)
    maxrejects = Max rejected targets, 0=ignore, default 32.
    log_name: string to specify log filename
    usersort = Enable if input fasta not sorted by length purposefully, lest
     usearch will raise an error.
    HALT_EXEC: Used for debugging app controller
    save_intermediate_files: Preserve all intermediate files created."""
     
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='usearch_fasta_dereplicated', suffix='.fasta')
    
    # using abspath to create log filepath, for cluster environments
    tmp_working_dir = split(abspath(output_filepath))[0]
    
    log_filepath = tmp_working_dir + "/" + log_name
    
    uc_filepath = tmp_working_dir + "/derep.uc"
    
    params = {'--derep_subseq':True,
              '--minlen':minlen,
              '--w':w,
              '--slots':slots,
              '--sizeout':sizeout,
              '--maxrejects':maxrejects}
    
    app = Usearch(params,HALT_EXEC=HALT_EXEC)
    
    if usersort:
        app.Parameters['--usersort'].on()
    
    app_result = app({'--cluster':fasta_filepath,
                      '--uc':uc_filepath,
                      '--seedsout':output_filepath,
                      '--log':log_filepath
                      })
                      
    if not save_intermediate_files:
        remove_files([uc_filepath])
    
    # Returning output filepath to delete if specified.
    
    return app_result, output_filepath
    
    
def usearch_sort_by_abundance(
    fasta_filepath,
    output_filepath = None,
    sizein = True,
    sizeout = True,
    minsize = 0,
    log_name = "abundance_sort.log",
    usersort = False,
    HALT_EXEC=False,
    save_intermediate_files=False):
    """ Sorts fasta file by abundance
    
    fasta_filepath = input fasta file, generally a dereplicated fasta
    output_filepath = output abundance sorted fasta filepath
    sizein = not defined in usearch helpstring
    sizeout = not defined in usearch helpstring
    minsize = minimum size of cluster to retain. 
    log_name = string to specify log filename
    usersort = Use if not sorting by abundance or usearch will raise an error
    HALT_EXEC: Used for debugging app controller
    save_intermediate_files: Preserve all intermediate files created.
    """
    
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='usearch_abundance_sorted', suffix='.fasta')
    
    # using abspath to create log filepath, for cluster environments
    tmp_working_dir = split(abspath(output_filepath))[0]
    
    log_filepath = tmp_working_dir + "/" + log_name
    
    params = {}
    
    app = Usearch(params,HALT_EXEC=HALT_EXEC)
    
    if usersort:
        app.Parameters['--usersort'].on()
        
    if minsize:
        app.Parameters['--minsize'].on(minsize)
        
    if sizein:
        app.Parameters['--sizein'].on()
        
    if sizeout:
        app.Parameters['--sizeout'].on()
    
    
    app_result = app({'--sortsize':fasta_filepath,
                      '--output':output_filepath,
                      '--log':log_filepath
                      })
    
    return app_result, output_filepath

def usearch_cluster_error_correction(
    fasta_filepath,
    output_filepath = None,
    percent_id_err = 0.97,
    sizein = True,
    sizeout = True,
    w=64,
    slots=16769023,
    maxrejects=64,
    log_name = "usearch_cluster_err_corrected.log",
    usersort = False,
    HALT_EXEC=False,
    save_intermediate_files=False):
    """ Cluster for err. correction at percent_id_err, output consensus fasta
    
    fasta_filepath = input fasta file, generally a dereplicated fasta
    output_filepath = output error corrected fasta filepath
    percent_id_err = minimum identity percent.
    sizein = not defined in usearch helpstring
    sizeout = not defined in usearch helpstring
    w = Word length for U-sorting
    slots = Size of compressed index table. Should be prime, e.g. 40000003.
     Should also specify --w, typical is --w 16 or --w 32.
    maxrejects = Max rejected targets, 0=ignore, default 32.
    log_name = string specifying output log name
    usersort = Enable if input fasta not sorted by length purposefully, lest
     usearch will raise an error.
    HALT_EXEC: Used for debugging app controller
    save_intermediate_files: Preserve all intermediate files created.
    """
    
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='usearch_cluster_err_corrected', suffix='.fasta')
    
    # using abspath to create log filepath, for cluster environments
    tmp_working_dir = split(abspath(output_filepath))[0]
    
    log_filepath = tmp_working_dir + "/" + log_name
    
    params = {'--sizein':sizein,
              '--sizeout':sizeout,
              '--id':percent_id_err,
              '--w':w,
              '--slots':slots,
              '--maxrejects':maxrejects}
    
    app = Usearch(params,HALT_EXEC=HALT_EXEC)
    
    if usersort:
        app.Parameters['--usersort'].on()
    
    
    app_result = app({'--cluster':fasta_filepath,
                      '--consout':output_filepath,
                      '--log':log_filepath
                      })
    
    return app_result, output_filepath


def usearch_chimera_filter_de_novo(
    fasta_filepath,
    output_chimera_filepath = None,
    output_non_chimera_filepath = None,
    abundance_skew = 2,
    log_name = "uchime_de_novo_chimera_filtering.log",
    usersort = False,
    HALT_EXEC=False,
    save_intermediate_files=False):
    """ Chimera filter de novo, output chimeras and non-chimeras to fastas
    
    fasta_filepath = input fasta file, generally a dereplicated fasta
    output_chimera_filepath = output chimera filepath
    output_non_chimera_filepath = output non chimera filepath
    abundance_skew = abundance skew setting for de novo filtering.
    usersort = Enable if input fasta not sorted by length purposefully, lest
     usearch will raise an error.
    HALT_EXEC: Used for debugging app controller
    save_intermediate_files: Preserve all intermediate files created.
    """
    
    output_chimera_filepath = output_chimera_filepath or \
     get_tmp_filename(prefix='uchime_chimeras_', suffix='.fasta')
     
    output_non_chimera_filepath = output_non_chimera_filepath or \
     get_tmp_filename(prefix='uchime_non_chimeras_', suffix='.fasta')
    
    # using abspath to create log filepath, for cluster environments
    tmp_working_dir = split(abspath(output_chimera_filepath))[0]
    
    log_filepath = tmp_working_dir + "/" + log_name
    
    params = {'--abskew':abundance_skew}
    
    app = Usearch(params,HALT_EXEC=HALT_EXEC)
    
    if usersort:
        app.Parameters['--usersort'].on()
    
    
    app_result = app({'--uchime':fasta_filepath,
                      '--chimeras':output_chimera_filepath,
                      '--nonchimeras':output_non_chimera_filepath,
                      '--log':log_filepath
                      })
    
    return app_result, output_non_chimera_filepath


def usearch_chimera_filter_ref_based(
    fasta_filepath,
    db_filepath,
    output_chimera_filepath = None,
    output_non_chimera_filepath = None,
    rev = True,
    log_name = "uchime_reference_chimera_filtering.log",
    usersort = False,
    HALT_EXEC=False,
    save_intermediate_files=False):
    """ Chimera filter against a reference database.
    
    fasta_filepath = input fasta file, generally a dereplicated fasta
    db_filepath = filepath to reference sequence database
    output_chimera_filepath = output chimera filepath
    output_non_chimera_filepath = output non chimera filepath
    rev = search plus and minus strands of sequences
    abundance_skew = abundance skew setting for de novo filtering.
    log_name = string specifying log filename.
    usersort = Enable if input fasta not sorted by length purposefully, lest
     usearch will raise an error.
    HALT_EXEC: Used for debugging app controller
    save_intermediate_files: Preserve all intermediate files created.
    """
    
    output_chimera_filepath = output_chimera_filepath or \
     get_tmp_filename(prefix='uchime_chimeras_', suffix='.fasta')
     
    output_non_chimera_filepath = output_non_chimera_filepath or \
     get_tmp_filename(prefix='uchime_non_chimeras_', suffix='.fasta')
    
    # using abspath to create log filepath, for cluster environments
    tmp_working_dir = split(abspath(output_chimera_filepath))[0]
    
    log_filepath = tmp_working_dir + "/" + log_name
    
    # clusters filepath created by usearch
    cluster_filepath = tmp_working_dir + "/refdb.uc"
    
    params = {'--rev':rev}
    
    app = Usearch(params,HALT_EXEC=HALT_EXEC)
    
    if usersort:
        app.Parameters['--usersort'].on()
    
    
    app_result = app({'--uchime':fasta_filepath,
                      '--db':db_filepath,
                      '--chimeras':output_chimera_filepath,
                      '--nonchimeras':output_non_chimera_filepath,
                      '--log':log_filepath,
                      '--uchimeout':cluster_filepath
                      })
                      
    if not save_intermediate_files:
        remove_files([cluster_filepath])
    
    return app_result, output_non_chimera_filepath
    

def usearch_cluster_seqs(
    fasta_filepath,
    output_filepath = None,
    percent_id = 0.97,
    sizein = True,
    sizeout = True,
    w=64,
    slots=16769023,
    maxrejects=64,
    log_name = "usearch_cluster_seqs_post_chimera.log",
    usersort = True,
    HALT_EXEC=False,
    save_intermediate_files=False):
    """ Cluster for err. correction at percent_id_err, output consensus fasta
    
    fasta_filepath = input fasta file, generally a dereplicated fasta
    output_filepath = output error corrected fasta filepath
    percent_id = minimum identity percent.
    sizein = not defined in usearch helpstring
    sizeout = not defined in usearch helpstring
    w = Word length for U-sorting
    slots = Size of compressed index table. Should be prime, e.g. 40000003.
     Should also specify --w, typical is --w 16 or --w 32.
    maxrejects = Max rejected targets, 0=ignore, default 32.
    log_name = string specifying output log name
    usersort = Enable if input fasta not sorted by length purposefully, lest
     usearch will raise an error.  In post chimera checked sequences, the seqs
     are sorted by abundance, so this should be set to True.
    HALT_EXEC: Used for debugging app controller
    save_intermediate_files: Preserve all intermediate files created.
    """
    
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='usearch_cluster', suffix='.fasta')
    
    # using abspath to create log filepath, for cluster environments
    tmp_working_dir = split(abspath(output_filepath))[0]
    
    log_filepath = tmp_working_dir + "/" + log_name
    
    uc_filepath = tmp_working_dir + "/clustered_seqs_post_chimera.uc"
    
    params = {'--sizein':sizein,
              '--sizeout':sizeout,
              '--id':percent_id,
              '--w':w,
              '--slots':slots,
              '--maxrejects':maxrejects}
    
    app = Usearch(params,HALT_EXEC=HALT_EXEC)
    
    if usersort:
        app.Parameters['--usersort'].on()
    
    
    app_result = app({'--cluster':fasta_filepath,
                      '--seedsout':output_filepath,
                      '--log':log_filepath,
                      '--uc':uc_filepath
                      })
                      
    if not save_intermediate_files:
        remove_files([uc_filepath])
    
    return app_result, output_filepath
    
def enumerate_otus(fasta_filepath,
                   output_filepath = None,
                   label_prefix = "",
                   label_suffix = "",
                   retain_label_as_comment = False,
                   count_start = 0):
    """ Writes unique, sequential count to OTUs
    
    fasta_filepath = input fasta filepath
    output_filepath = output fasta filepath
    label_prefix = string to place before enumeration
    label_suffix = string to place after enumeration
    retain_label_as_comment = if True, will place existing label in sequence
     comment, after a tab
    count_start = number to start enumerating OTUs with
    
    """
    
    fasta_i = open(fasta_filepath, "U")
    
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='enumerated_seqs_', suffix='.fasta')
     
    fasta_o = open(output_filepath, "w")
    
    for label, seq in MinimalFastaParser(fasta_i):
        curr_label = ">" + label_prefix + str(count_start) + label_suffix
        if retain_label_as_comment:
            curr_label += '\t' + label
        fasta_o.write(curr_label.strip() + '\n')
        fasta_o.write(seq.strip() + '\n')
        count_start += 1
        
    
    return output_filepath
    

    
def assign_reads_to_otus(original_fasta,
                         filtered_fasta,
                         output_filepath,
                         log_name = "assign_reads_to_otus.log",
                         perc_id = 0.97,
                         global_alignment = True,
                         HALT_EXEC=False,
                         save_intermediate_files=False):
    """ Uses original fasta file, blasts to assign reads to filtered fasta
    
    original_fasta = filepath to original query fasta
    filtered_fasta = filepath to enumerated, filtered fasta
    output_filepath = output path to clusters (uc) file
    log_name = string specifying output log name
    usersort = Enable if input fasta not sorted by length purposefully, lest
     usearch will raise an error.  In post chimera checked sequences, the seqs
     are sorted by abundance, so this should be set to True.
    HALT_EXEC: Used for debugging app controller
    save_intermediate_files: Preserve all intermediate files created.
    """
    
    # Not sure if I feel confortable using blast as a way to recapitulate 
    # original read ids....
    
    output_filepath = output_filepath or \
     get_tmp_filename(prefix='assign_reads_to_otus', suffix='.uc')
    
    # using abspath to create log filepath, for cluster environments
    tmp_working_dir = split(abspath(output_filepath))[0]
    
    log_filepath = tmp_working_dir + "/" + log_name
    
    
    params = {'--id':perc_id,
              '--global':global_alignment}
    
    app = Usearch(params,HALT_EXEC=HALT_EXEC)
    
    
    app_result = app({'--query':original_fasta,
                      '--db':filtered_fasta,
                      '--log':log_filepath,
                      '--uc':output_filepath
                      })
                      
    
    return app_result, output_filepath


'''def uclust_cluster_from_sorted_fasta_filepath(
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
    HALT_EXEC=False):
    """ Returns clustered uclust file from sorted fasta"""
    output_filepath = uc_save_filepath or \
     get_tmp_filename(prefix='uclust_clusters',suffix='.uc')
     
    
    params = {'--id':percent_ID,
              '--maxaccepts':max_accepts,
              '--maxrejects':max_rejects,
              '--stepwords':stepwords,
              '--w':word_length}
    app = Uclust(params,HALT_EXEC=HALT_EXEC)
    
    # Set any additional parameters specified by the user
    if enable_rev_strand_matching: app.Parameters['--rev'].on()
    if optimal: app.Parameters['--optimal'].on()
    if exact: app.Parameters['--exact'].on()
    if suppress_sort: app.Parameters['--usersort'].on()
    if subject_fasta_filepath: app.Parameters['--lib'].on(subject_fasta_filepath)
    if suppress_new_clusters: app.Parameters['--libonly'].on()
    if stable_sort: app.Parameters['--stable_sort'].on()
    
    app_result = app({'--input':fasta_filepath,'--uc':output_filepath})
    return app_result '''


def get_output_filepaths(output_dir, fasta_filepath):
    """ Returns filepaths for intermediate file to be kept """
    
    if not output_dir.endswith('/'):
        output_dir += '/'
        
    output_file_basename = "".join(basename(fasta_filepath).split('.')[0:-1])
    us_save_filepath = '%s%s_clusters.uc' % \
     (output_dir, output_file_basename)

    return us_save_filepath




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
        try:
            sort_fasta.cleanUp()
        except AttributeError:
            pass
        uclust_cluster.cleanUp()
    
    if return_cluster_maps:
        return clusters, failures, seeds
    else:
        return clusters.values(), failures, seeds

## End uclust convenience functions
