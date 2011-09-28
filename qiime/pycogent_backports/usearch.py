#!/usr/bin/env python
"""Application controller for usearch v5.0.144

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

from os.path import split, splitext, basename, isdir, abspath, isfile

from cogent.parse.fasta import MinimalFastaParser
from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, ResultPath,\
 get_tmp_filename, ApplicationError, ApplicationNotFoundError
from cogent.util.misc import remove_files

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
        
        # Output filename will be in uclust (.uc) format
        # Output cluster file, required parameter
        '--uc':ValuedParameter('--',Name='uc',Delimiter=' ',
            IsPath=True),
        
        # ID percent for OTU, by default is 97%
        '--id':ValuedParameter('--',Name='id',Delimiter=' ',IsPath=False),

        # Disable reverse comparison option, if norev is disabled
        # memory usage is expected to double for uclust
        '--rev':FlagParameter('--',Name='rev'),
        
        # Maximum hits before quitting search (default 1, 0=infinity).
        '--maxaccepts':ValuedParameter('--',Name='maxaccepts',Delimiter=' '),
        
        # Maximum rejects before quitting search (default 8, 0=infinity). 
        '--maxrejects':ValuedParameter('--',Name='maxrejects',Delimiter=' '),
        
        # Target nr. of common words (default 8, 0=don't step)
        '--stepwords':ValuedParameter('--',Name='stepwords',Delimiter=' '),
        
        # Word length for windex (default 5 aa.s, 8 nuc.s).
        '--w':ValuedParameter('--',Name='w',Delimiter=' '),
        
        # Don't assume input is sorted by length (default assume sorted).
        '--usersort':FlagParameter('--',Name='usersort'),
        
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
        allowed_values =  ['--uc', '--output', '--mergesort', '--log',\
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

## Start functions for processing usearch output files



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
        
                
    

## End functions for processing usearch output files


## Start usearch convenience functions

def usearch_fasta_sort_from_filepath(
    fasta_filepath,
    output_filepath=None,
    log_name = "sortlen.log",
    HALT_EXEC=False,
    save_intermediate_files=False,
    remove_usearch_logs=False):
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
    
    params = {}
            
    app = Usearch(params,HALT_EXEC=HALT_EXEC)
    
    data={'--mergesort':fasta_filepath,\
          '--output':output_filepath,
         }
         
    if not remove_usearch_logs:
        data['--log'] = log_filepath
    
    app_result = app(data)
                           
    
    
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
    save_intermediate_files=False,
    remove_usearch_logs=False):
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
        
    data = {'--cluster':fasta_filepath,
                      '--uc':uc_filepath,
                      '--seedsout':output_filepath
            }
    
    if not remove_usearch_logs:
        data['--log'] = log_filepath
    
    app_result = app(data)

                      
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
    save_intermediate_files=False,
    remove_usearch_logs=False):
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
    
    log_filepath = tmp_working_dir + "/" + "minsize_" + str(minsize) +\
     "_" + log_name
    
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
        
        
    data = {'--sortsize':fasta_filepath,
            '--output':output_filepath
           }
           
    if not remove_usearch_logs:
        data['--log'] = log_filepath
    
    # Can have no data following this filter step, which will raise an 
    # application error, try to catch it here to raise meaningful message.
    try:
        app_result = app(data)
    except ApplicationError:
        raise ValueError, ('No data following filter steps, please check '+\
         'parameter settings for OTUPipe.')
    
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
    save_intermediate_files=False,
    remove_usearch_logs=False):
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
        
    data = {'--cluster':fasta_filepath,
            '--consout':output_filepath
            }
    
    if not remove_usearch_logs:
        data['--log'] = log_filepath
    
    
    app_result = app(data)
    
    return app_result, output_filepath


def usearch_chimera_filter_de_novo(
    fasta_filepath,
    output_chimera_filepath = None,
    output_non_chimera_filepath = None,
    abundance_skew = 2,
    log_name = "uchime_de_novo_chimera_filtering.log",
    usersort = False,
    HALT_EXEC=False,
    save_intermediate_files=False,
    remove_usearch_logs=False):
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
        
    data = {'--uchime':fasta_filepath,
            '--chimeras':output_chimera_filepath,
            '--nonchimeras':output_non_chimera_filepath
            }
        
    if not remove_usearch_logs:
        data['--log'] = log_filepath
    
    
    app_result = app(data)
                      
    if not save_intermediate_files:
        remove_files([output_chimera_filepath])
    
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
    save_intermediate_files=False,
    remove_usearch_logs=False):
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
        
    data = {'--uchime':fasta_filepath,
            '--db':db_filepath,
            '--chimeras':output_chimera_filepath,
            '--nonchimeras':output_non_chimera_filepath,
            '--uchimeout':cluster_filepath
            }
    
    if not remove_usearch_logs:
        data['--log'] = log_filepath

    app_result = app(data)
                      
    if not save_intermediate_files:
        remove_files([cluster_filepath, output_chimera_filepath])
    
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
    log_name = "usearch_cluster_seqs.log",
    usersort = True,
    HALT_EXEC=False,
    save_intermediate_files=False,
    remove_usearch_logs=False):
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
        
    data = {'--cluster':fasta_filepath,
            '--seedsout':output_filepath,
            '--uc':uc_filepath
            }
        
    if not remove_usearch_logs:
        data['--log'] = log_filepath
    
    
    app_result = app(data)
                      
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
                         output_filepath = None,
                         log_name = "assign_reads_to_otus.log",
                         perc_id_blast = 0.97,
                         global_alignment = True,
                         HALT_EXEC=False,
                         save_intermediate_files=False,
                         remove_usearch_logs=False):
    """ Uses original fasta file, blasts to assign reads to filtered fasta
    
    original_fasta = filepath to original query fasta
    filtered_fasta = filepath to enumerated, filtered fasta
    output_filepath = output path to clusters (uc) file
    log_name = string specifying output log name
    perc_id_blast = percent ID for blasting original seqs against filtered set
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
    
    
    params = {'--id':perc_id_blast,
              '--global':global_alignment}
    
    
    
    app = Usearch(params,HALT_EXEC=HALT_EXEC)
    
    data = {'--query':original_fasta,
            '--db':filtered_fasta,
            '--uc':output_filepath
            }
    
    if not remove_usearch_logs:
        data['--log'] = log_filepath
    
    
    app_result = app(data)
                      
    
    return app_result, output_filepath



def otu_pipe(
    fasta_filepath,
    output_dir = None,
    percent_id = 0.97,
    percent_id_err = 0.97,
    minsize = 4,
    abundance_skew = 2,
    db_filepath = None,
    rev = True,
    label_prefix = "",
    label_suffix = "",
    retain_label_as_comment = False,
    count_start = 0,
    perc_id_blast = 0.97,
    save_intermediate_files = False,
    HALT_EXEC = False,
    global_alignment = True,
    sizein = True,
    sizeout = True,
    w=64,
    slots=16769023,
    maxrejects=64,
    minlen=64,
    de_novo_chimera_detection=True,
    reference_chimera_detection=True,
    cluster_size_filtering=True,
    remove_usearch_logs=False,
    usersort=True
    ):
        
    """ Main convenience wrapper for using usearch to filter/cluster seqs
    
    The complete 'OTU Pipe' process is a multistep process with many calls
    to usearch with various parameters.  It is likely to change from the 
    original implementation.  A lot.
    
    fasta_filepath = fasta filepath to filtering/clustering (e.g., output 
     seqs.fna file from split_libraries.py)
    output_dir = directory to store the otu mapping file, as well logs and
     the intermediate files created if save_intermediate_files is True.
    percent_ID = percent ID for clustering sequences.
    percent_ID_err = percent ID for filtering out chimeras
    minsize = Minimum size of cluster for retention after chimera removal.
    abundance_skew = threshold setting for chimera removal with de novo
     chimera detection.
    db_filepath = filepath of reference fasta sequence set for ref based 
     chimera detection.
    rev = search plus and minus strands of sequences, used in ref based chimera
     detection.
    label_prefix = optional prefix added to filtered fasta file.
    label_suffix = optional suffix added to filtered fasta file.
    retain_label_as_comment = option to add usearch generated label to 
     enumerated fasta labels.
    count_start = integer to begin counting at for sequence enumeration.
    perc_id_blast = percent identity setting for using blast algorithm to 
     assign original sequence labels to filtered fasta.
    global_alignment = Setting for assignment of original seq labels to filtered
     seqs.
    sizein = not defined in usearch helpstring
    sizeout = not defined in usearch helpstring
    w = Word length for U-sorting
    slots = Size of compressed index table. Should be prime, e.g. 40000003.
     Should also specify --w, typical is --w 16 or --w 32.
    maxrejects = Max rejected targets, 0=ignore, default 32.
    save_intermediate_files = retain all the intermediate files created during
     this process.
    minlen = (not specified in usearch helpstring), but seems like a good bet
     that this refers to the minimum length of the sequences for dereplication.
    HALT_EXEC = used to debug app controller problems.
    de_novo_chimera_detection = If True, will detect chimeras de novo
    reference_chimera_detection = If True, will detect chimeras ref based
    cluster_size_filtering = If True, will filter OTUs according to seq counts.
    remove_usearch_logs = If True, will not call the --log function for each
     usearch call.
    usersort = Used for specifying custom sorting (i.e., non-length based
     sorting) with usearch/uclust.
    """
    
    
    # Save a list of intermediate filepaths in case they are to be removed.
    intermediate_files = []
    
    
    if output_dir and not output_dir.endswith('/'):
        output_dir += '/'
        


    try:
        
        # Sort seqs by length
        app_result, output_filepath_len_sorted =\
         usearch_fasta_sort_from_filepath(fasta_filepath, output_filepath =\
         output_dir + 'len_sorted.fasta',
         save_intermediate_files=save_intermediate_files,
         remove_usearch_logs=remove_usearch_logs)
         
        intermediate_files.append(output_filepath_len_sorted)
        
        # Dereplicate sequences
        app_result, output_filepath_dereplicated =\
         usearch_dereplicate_exact_subseqs(output_filepath_len_sorted,
         output_filepath = output_dir + 'dereplicated_seqs.fasta',
         minlen=minlen, w=w, slots=slots, sizeout=sizeout, 
         maxrejects=maxrejects, save_intermediate_files=save_intermediate_files,
         remove_usearch_logs=remove_usearch_logs)
        
        intermediate_files.append(output_filepath_dereplicated)
        
        # Sort by abundance, initially no filter based on seqs/otu
        app_result, output_fp =\
         usearch_sort_by_abundance(output_filepath_dereplicated,
         output_filepath = output_dir + 'abundance_sorted.fasta',
         usersort = True, sizein=sizein, sizeout=sizeout, minsize=0,
         remove_usearch_logs=remove_usearch_logs)
         
        intermediate_files.append(output_fp)
        
        
        app_result, output_fp =\
             usearch_cluster_error_correction(output_fp,
             output_filepath = output_dir + 'clustered_error_corrected.fasta',
             usersort = True, percent_id_err=percent_id_err, sizein=sizein,
             sizeout=sizeout, w=w, slots=slots, maxrejects=maxrejects,
             remove_usearch_logs=remove_usearch_logs)

        intermediate_files.append(output_fp)
            
        # Series of conditional tests, using generic 'output_fp' name so the
        # conditional filtering, if any/all are selected, do not matter.
        
        if reference_chimera_detection:
            app_result, output_fp =\
             usearch_chimera_filter_ref_based(output_fp,
             db_filepath=db_filepath, output_chimera_filepath= output_dir +\
             'reference_chimeras.fasta', output_non_chimera_filepath =\
             output_dir + 'reference_non_chimeras.fasta', usersort=True, 
             save_intermediate_files=save_intermediate_files, rev=rev,
             remove_usearch_logs=remove_usearch_logs)
        
            intermediate_files.append(output_fp)
        
        if de_novo_chimera_detection:
            app_result, output_fp =\
             usearch_chimera_filter_de_novo(output_fp, 
             abundance_skew = abundance_skew, output_chimera_filepath =\
             output_dir + 'de_novo_chimeras.fasta',
             output_non_chimera_filepath = output_dir +\
             'de_novo_non_chimeras.fasta', usersort=True,
             save_intermediate_files=save_intermediate_files,
             remove_usearch_logs=remove_usearch_logs)
        
            intermediate_files.append(output_fp)
            

            
        if cluster_size_filtering:
            
            # Test for empty filepath following filters, raise error if all seqs
            # have been removed
            
            app_result, output_fp =\
             usearch_sort_by_abundance(output_fp, output_filepath = output_dir+\
             'abundance_sorted_minsize_' + str(minsize) + '.fasta', 
             minsize=minsize, sizein=sizein, sizeout=sizeout,
             remove_usearch_logs=remove_usearch_logs)
             
             
            intermediate_files.append(output_fp)

        # cluster seqs
        # Should we add in option to use alternative OTU picking here?
        # Seems like it will be a bit of a mess...maybe after we determine
        # if OTU pipe should become standard.
        app_result, output_filepath_clustered_seqs =\
         usearch_cluster_seqs(output_fp, output_filepath = output_dir +\
         'clustered_seqs.fasta', percent_id=percent_id, sizein=sizein,
         sizeout=sizeout, w=w, slots=slots, maxrejects=maxrejects,
         save_intermediate_files=save_intermediate_files,
         remove_usearch_logs=remove_usearch_logs)

        intermediate_files.append(output_filepath_clustered_seqs)
        
        # Enumerate the OTUs in the clusters
        output_filepath_enum_otus =\
         enumerate_otus(output_filepath_clustered_seqs, output_filepath =\
         output_dir + 'enumerated_otus.fasta', label_prefix=label_prefix,
         label_suffix=label_suffix, count_start=count_start,
         retain_label_as_comment=retain_label_as_comment)
        
        intermediate_files.append(output_filepath_enum_otus)
        
        # Get original sequence label identities
        app_result, clusters_file = assign_reads_to_otus(fasta_filepath,
         output_filepath_enum_otus, output_filepath = output_dir +\
         'assign_reads_to_otus.uc', perc_id_blast=perc_id_blast,
         global_alignment=global_alignment,
         remove_usearch_logs=remove_usearch_logs)
         
        intermediate_files.append(clusters_file)


        
    except ApplicationError:
        raise ApplicationError, ('Error running usearch. Possible causes are '
         'unsupported version (current supported version is usearch '+\
         'v5.0.144) is installed or improperly formatted input file was provided')
    except ApplicationNotFoundError:
        remove_files(files_to_remove)
        raise ApplicationNotFoundError('usearch not found, is it properly '+\
         'installed?')

    # Get dict of clusters, list of failures
    clusters, failures = clusters_from_blast_uc_file(open(clusters_file, "U"))
    
    # Remove temp files unless user specifies output filepath
    if not save_intermediate_files:
        remove_files(intermediate_files)
        
    return clusters, failures
    

## End uclust convenience functions
