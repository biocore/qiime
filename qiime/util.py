#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Rob Knight", "Daniel McDonald", "Greg Caporaso", 
"Justin Kuczynski","Jens Reeder","Catherine Lozupone"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"


"""Contains general utility code in support of the Qiime project.

A lot of this might migrate into cogent at some point.
"""

from StringIO import StringIO
from os import getenv, makedirs
from operator import itemgetter
from os.path import abspath, exists, dirname, join, isdir
from collections import defaultdict
import gzip
import sys
import os
from copy import deepcopy
from datetime import datetime
from subprocess import Popen, PIPE, STDOUT
from numpy import min, max, median, mean
import numpy
from numpy.ma import MaskedArray
from numpy.ma.extras import apply_along_axis
from numpy import array, zeros, argsort, shape, vstack,ndarray, asarray, \
        float, where, isnan, mean, std
        
from cogent.util.dict2d import Dict2D
from cogent import LoadSeqs, Sequence,DNA
from cogent.cluster.procrustes import procrustes
from cogent.core.alignment import Alignment
from cogent.core.moltype import MolType, IUPAC_DNA_chars, IUPAC_DNA_ambiguities,\
    IUPAC_DNA_ambiguities_complements, DnaStandardPairs, ModelDnaSequence
from cogent.data.molecular_weight import DnaMW
from cogent.core.sequence import DnaSequence
from cogent.app.blast import Blastall
from cogent.app.util import get_tmp_filename as cogent_get_tmp_filename, FilePath
from cogent.parse.blast import BlastResult
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser
from cogent.util.misc import remove_files
from cogent.util.dict2d import Dict2D
from cogent.app.formatdb import build_blast_db_from_fasta_path,\
    build_blast_db_from_fasta_file
from cogent import LoadSeqs
from cogent.util.misc import (create_dir, 
                              handle_error_codes)
from qiime import __version__ as qiime_library_version
from cogent.util.option_parsing import (parse_command_line_parameters,
                                        make_option)
from qiime.parse import (parse_otu_table,
                         parse_qiime_config_files,
                         parse_coords,
                         parse_newick,
                         PhyloNode,
                         parse_mapping_file,
                         parse_denoiser_mapping)
from qiime.format import format_otu_table

class TreeMissingError(IOError):
    """Exception for missing tree file"""
    pass

class OtuMissingError(IOError):
    """Exception for missing OTU file"""
    pass

class AlignmentMissingError(IOError):
    """Exception for missing alignment file"""
    pass

class MissingFileError(IOError):
    pass

class FileFormatError(IOError):
    """Exception for wrong file format"""
    pass

def make_safe_f(f, allowed_params):
    """Make version of f that ignores extra named params."""
    def inner(*args, **kwargs):
        if kwargs:
            new_kwargs = {}
            for k, v in kwargs.items():
                if k in allowed_params:
                    new_kwargs[k] = v
            return f(*args, **new_kwargs)
        return f(*args, **kwargs)
    return inner

class FunctionWithParams(object):
    """A FunctionWithParams is a replacement for the function factory.
    
    Specifically, the params that will be used in the __call__ method are
    available in a dict so you can keep track of them with the object
    itself.
    """
    Application = None
    Algorithm = None
    Citation = None
    Params = {}
    Name = 'FunctionWithParams' #override in subclasses
    _tracked_properties = []    #properties tracked like params

    def __init__(self, params):
        """Return new FunctionWithParams object with specified params.
        
        Note: expect params to contain both generic and per-method (e.g. for
        cdhit) params, so leaving it as a dict rather than setting
        attributes. 
        
        Some standard entries in params are:

        [fill in on a per-application basis]
        """
        self.Params.update(params)
        self._tracked_properties.extend(['Application','Algorithm','Citation'])

    def __str__(self):
        """Returns formatted key-value pairs from params."""
        res = [self.Name + ' parameters:']
        for t in self._tracked_properties:
            res.append(t + ':' + str(getattr(self, t)))
        for k, v in sorted(self.Params.items()):
            res.append(str(k) + ':' + str(v))
        return '\n'.join(res)

    def writeLog(self, log_path):
        """Writes self.Params and other relevant info to supplied path."""
        f=open(log_path, 'w')
        f.write(str(self))
        f.close()

    def getResult(self, *args, **kwargs):
        """Gets result in __call__. Override in subclasses."""
        return None

    def formatResult(self, result):
        """Formats result as string (for whatever "result" means)."""
        return str(result)

    def writeResult(self, result_path, result):
        """Writes result to result_path. May need to format in subclasses."""
        f = open(result_path, 'w')
        f.write(self.formatResult(result))
        f.close()

    def getOtuTable(self, otu_source):
        """Returns parsed OTU table from putative OTU source."""
        
        #if we have a string starting with #, assume it's an OTU file,
        #otherwise assume it's a path
        # if 4-tuple, just return it
        if type(otu_source) == type((1,3,4,44)):
            return otu_source
        if hasattr(otu_source, 'startswith') and otu_source.startswith('#'):
            try:
                return parse_otu_table(StringIO(otu_source))
            except (TypeError, ValueError), e:
                raise OtuMissingError, \
                    "Tried to read OTUs from string starting with # but got "+e
        else:
            try:
                otu_file = open(otu_source, 'U')
            except (TypeError, IOError):
                raise OtuMissingError, \
                    "Couldn't read OTU file at path: %s" % otu_source
            result = parse_otu_table(otu_file)
            otu_file.close()
            return result

    def getTree(self, tree_source):
        """Returns parsed tree from putative tree source"""
        if isinstance(tree_source, PhyloNode):
            tree = tree_source    #accept tree object directly for tests
        elif tree_source:
            try:
                f = open(tree_source, 'U')
            except (TypeError, IOError):
                raise TreeMissingError, \
                    "Couldn't read tree file at path: %s" % tree_source
            tree = parse_newick(f, PhyloNode)
            f.close()
        else:
            raise TreeMissingError, str(self.Name) + \
                " is a phylogenetic metric, but no tree was supplied."
        return tree

    def getData(self, data_source):
        """Returns data from putative source, which could be a path"""
        if isinstance(data_source, str):
            try:
                return eval(data_source)
            except (NameError, SyntaxError):
                try:
                    data_f = open(data_source, 'U')
                    data = data_f.read()
                    data_f.close()
                    try:
                        return eval(data)
                    except (NameError, SyntaxError, TypeError):
                        pass
                    return data
                except (IOError, NameError, TypeError):
                    pass
        #if we got here, either we didn't get a string or we couldn't read
        #the data source into any other kind of object
        return data_source

    def getAlignment(self, aln_source):
        """Returns parsed alignment from putative alignment source"""
        if isinstance(aln_source, Alignment):
            aln = aln_source
        elif aln_source:
            try:
                aln = LoadSeqs(aln_source, Aligned=True)
            except (TypeError, IOError, AssertionError):
                raise AlignmentMissingError, \
                    "Couldn't read alignment file at path: %s" % aln_source
        else:
            raise AlignmentMissingError, str(self.Name) + \
                " requires an alignment, but no alignment was supplied."
        return aln

    def __call__ (self, result_path=None, log_path=None,\
        *args, **kwargs):
        """Returns the result of calling the function using the params dict.
        
        Parameters:
        [fill in on a per-application basis]
        """
        result = self.getResult(*args, **kwargs)
        if log_path:
            self.writeLog(log_path)
        if result_path:
            self.writeResult(result_path, result)
        else:
            return result

def get_qiime_project_dir():
    """ Returns the top-level QIIME directory
    """
    # Get the full path of util.py
    current_file_path = abspath(__file__)
    # Get the directory containing util.py
    current_dir_path = dirname(current_file_path)
    # Return the directory containing the directory containing util.py
    return dirname(current_dir_path)
    
def get_qiime_scripts_dir():
    """ Returns the QIIME scripts directory 
    
        This value must be stored in qiime_config if the user
        has installed qiime using setup.py. If it is not in
        qiime_config, it is inferred from the qiime_project_dir.
    
    """
    qiime_config = load_qiime_config()
    qiime_config_value = qiime_config['qiime_scripts_dir']
    if qiime_config_value != None:
        result = qiime_config_value
    else:
        result = join(get_qiime_project_dir(),'scripts')
    
    #assert exists(result),\
    # "qiime_scripts_dir does not exist: %s." % result +\
    # " Have you defined it correctly in your qiime_config?"
    
    return result
    
def get_qiime_temp_dir():
    """ Returns the temp directory that should be used by QIIME scripts
    
    """
    qiime_config = load_qiime_config()
    qiime_config_value = qiime_config['temp_dir']
    if qiime_config_value != None:
        result = qiime_config_value
    else:
        result = '/tmp/'
    return result
    
def get_tmp_filename(tmp_dir=None, prefix="tmp", suffix=".txt",
    result_constructor=FilePath):
    """ Wrap cogent.app.util.get_tmp_filename to modify the default tmp_dir """
    if tmp_dir == None:
        tmp_dir = get_qiime_temp_dir()
    return cogent_get_tmp_filename(tmp_dir=tmp_dir,
                                   prefix=prefix,
                                   suffix=suffix,
                                   result_constructor=result_constructor)
    
def load_qiime_config():
    """Return default parameters read in from file"""
    
    qiime_config_filepaths = []
    qiime_project_dir = get_qiime_project_dir()
    qiime_config_filepaths.append(\
     qiime_project_dir + '/qiime/support_files/qiime_config')
    
    qiime_config_env_filepath = getenv('QIIME_CONFIG_FP')
    if qiime_config_env_filepath:
        qiime_config_filepaths.append(qiime_config_env_filepath)
    
    home_dir = getenv('HOME')
    if home_dir:
        qiime_config_home_filepath = home_dir + '/.qiime_config'
        qiime_config_filepaths.append(qiime_config_home_filepath)
    
    qiime_config_files = []
    for qiime_config_filepath in qiime_config_filepaths:
        if exists(qiime_config_filepath):
            qiime_config_files.append(open(qiime_config_filepath))
        
    return parse_qiime_config_files(qiime_config_files)

# The qiime_blast_seqs function should evetually move to PyCogent,
# but I want to test that it works for all of the QIIME functionality that
# I need first. -Greg
def qiime_blast_seqs(seqs,
     blast_constructor=Blastall,
     blast_program='blastn',
     blast_db=None,
     refseqs=None,
     refseqs_fp=None,
     blast_mat_root=None,
     params={},
     WorkingDir=None,
     seqs_per_blast_run=1000,
     HALT_EXEC=False):
    """Blast list of sequences.

    seqs: a list (or object with list-like interace) of (seq_id, seq) 
     tuples (e.g., the output of MinimalFastaParser)
    
    """
    assert blast_db or refseqs_fp or refseqs, \
     'Must provide either a blast_db or a fasta '+\
     'filepath containing sequences to build one.'
     
    if refseqs_fp:
        blast_db, db_files_to_remove =\
         build_blast_db_from_fasta_path(refseqs_fp,output_dir=WorkingDir)
    elif refseqs:
        blast_db, db_files_to_remove =\
         build_blast_db_from_fasta_file(refseqs,output_dir=WorkingDir)
    else:
        db_files_to_remove = []
    
    params["-d"] = blast_db
    params["-p"] = blast_program
           
    blast_app = blast_constructor(
                   params=params,
                   blast_mat_root=blast_mat_root,
                   InputHandler='_input_as_seq_id_seq_pairs',
                   WorkingDir=WorkingDir,
                   SuppressStderr=True,
                   HALT_EXEC=HALT_EXEC)

    current_seqs = []
    blast_results = BlastResult([])
    for seq in seqs:
        current_seqs.append(seq)
        if len(current_seqs) % seqs_per_blast_run == 0:
            if blast_results:
                blast_results.update(\
                 BlastResult(blast_app(current_seqs)['StdOut']))
            else:
                blast_results = BlastResult(blast_app(current_seqs)['StdOut'])
            current_seqs = []
    
    # clean-up run: blast the remaining sequences
    blast_results.update(\
     BlastResult(blast_app(current_seqs)['StdOut']))

    remove_files(db_files_to_remove)
    
    return blast_results

def extract_seqs_by_sample_id(seqs, sample_ids, negate=False):
    """ Returns (seq id, seq) pairs if sample_id is in sample_ids """
    sample_ids = {}.fromkeys(sample_ids)

    if not negate:
        def f(s):
            return s in sample_ids
    else:
        def f(s):
            return s not in sample_ids

    for seq_id, seq in seqs:
        sample_id = seq_id.split('_')[0]
        if f(sample_id):
            yield seq_id, seq
            
def split_fasta_on_sample_ids(seqs):
    """ yields (sample_id, seq_id, seq) for each entry in seqs 
    
        seqs: (seq_id,seq) pairs, as generated by MinimalFastaParser
    
    """
    for seq_id, seq in seqs:
        yield (seq_id.split()[0].rsplit('_',1)[0], seq_id, seq)
    return
        
def split_fasta_on_sample_ids_to_dict(seqs):
    """ return split_fasta_on_sample_ids as {sample_id: [(seq_id, seq), ], }
    
        seqs: (seq_id,seq) pairs, as generated by MinimalFastaParser
    
    """
    result = {}
    for sample_id,seq_id,seq in split_fasta_on_sample_ids(seqs):
        try:
            result[sample_id].append((seq_id,seq))
        except KeyError:
            result[sample_id] = [(seq_id,seq)]
    return result

def write_seqs_to_fasta(fp,seqs,write_mode='w'):
    """Write seqs to fp with specified write mode ('a' or 'w')
    
        seqs: list of (seq_id,seq) tuples, as obtained from
         MinimalFastaParser
    """
    f = open(fp,write_mode)
    for s in seqs:
        f.write('>%s\n%s\n' % (s))
    f.close()

def split_fasta_on_sample_ids_to_files(seqs,
                                       output_dir,
                                       per_sample_buffer_size=500):
    """ output of split_fasta_on_sample_ids to fasta in specified output_dir
    
        seqs: (seq_id,seq) pairs, as generated by MinimalFastaParser
        output_dir: string defining directory where output should be 
         written, will be created if it doesn't exist
         
        This function takes a buffered approach to writing files to 
        avoid hitting errors arising from too many files being open
        when working with large numbers of samples ids (e.g. > 1024 on linux)
    """
    create_dir(output_dir)
    file_lookup = {}
    for sample_id,seq_id,seq in split_fasta_on_sample_ids(seqs):
        # grab or create the list corresponding to the current sample id
        try:
            current_seqs = file_lookup[sample_id][1]
        except KeyError:
            current_fp = '%s/%s.fasta' % (output_dir,sample_id)
            if exists(current_fp):
                raise IOError,\
                 (" %s already exists. Will not perform split -- remove this"
                  " file or specify a different output directory." % current_fp)
            current_seqs = list()
            file_lookup[sample_id] = [current_fp,current_seqs]
        
        # append the current sequence to the current seqs list
        current_seqs.append((seq_id,seq))
        # compare current_seqs length to the buffer size, and write
        # if it has hit the buffer size
        if len(current_seqs) == per_sample_buffer_size:
            current_fp = file_lookup[sample_id][0]
            write_seqs_to_fasta(current_fp,current_seqs,write_mode='a')
            # reset the current seqs buffer
            file_lookup[sample_id][1] = list()
    
    for current_fp,current_seqs in file_lookup.values():
        write_seqs_to_fasta(current_fp,current_seqs,write_mode='a')
    return None

def compute_seqs_per_library_stats(otu_f):
    counts = []
    sample_ids, otu_ids, otu_table, lineages = parse_otu_table(otu_f)
    for i in range(otu_table.shape[1]):
        counts.append(sum(otu_table[:,i]))
        
    return min(counts), max(counts), median(counts), mean(counts),\
     dict(zip(sample_ids,counts))

def median_absolute_deviation(x):
    """ compute the median of the absolute deviations from the median """
    x = asarray(x)
    median_x = median(x)
    median_abs_deviation = median(abs(x - median_x))
    return median_abs_deviation, median_x

def guess_even_sampling_depth(counts_per_sample,num_deviations=2.25):
    """ guess a depth for even sampling 
    
        this is currently computed as the smallest seqs per sample 
         count >= the median seqs per sample count - (2.25 * the median absolute 
         deviation). 2.25 was chosen emprically by seeing how different values
         of num_deviations resulted in a choice that was similar to what
         I've chosen on several real OTU tables.
    """
    counts_per_sample.sort()
    median_abs_dev, median_count = \
     median_absolute_deviation(counts_per_sample)
    min_threshold = median_count - (num_deviations * median_abs_dev)
    for e in counts_per_sample:
        if e >= min_threshold:
            return e
    raise ValueError,\
     "No acceptable even sampling depth identified. "+\
     "It shouldn't be possible to get here, but just in case here's the " +\
     "counts per sample: %s" ' '.join(map(str,counts_per_sample)) 
    
    
def raise_error_on_parallel_unavailable(qiime_config=None):
    """Raise error if no parallel QIIME bc user hasn't set jobs_to_start
    """
    if qiime_config == None:
        qiime_config = load_qiime_config()
    if 'jobs_to_start' not in qiime_config or \
       int(qiime_config['jobs_to_start']) < 2:
       raise RuntimeError,\
        "Parallel QIIME is not available. (Have you set"+\
        " jobs_to_start to greater than 1 in your qiime_config?"
        


def get_options_lookup():
    """ Return dict of commonly used options """
    qiime_config = load_qiime_config()
    result = {}
    result['fasta_as_primary_input'] =\
     make_option('-i','--input_fasta_fp',type="existing_filepath",
      help='path to the input fasta file')
    result['otu_table_as_primary_input'] =\
     make_option('-i','--otu_table_fp',type="existing_filepath",
      help='path to the input OTU table (i.e., the output from make_otu_table.py)')
    result['otu_map_as_primary_input'] =\
     make_option('-i','--otu_map_fp',type="existing_filepath",
      help='path to the input OTU map (i.e., the output from pick_otus.py)')
    result['log_fp'] =\
     make_option('-l','--log_fp',type="new_filepath",
     help='path to write the log file')
    result['input_fasta'] =\
     make_option('-f','--input_fasta_fp',type="existing_filepath",
      help='path to the input fasta file')
    result['output_dir'] =\
     make_option('-o','--output_dir',type="new_dirpath",
     help='path to the output directory')
    result['output_fp'] =\
     make_option('-o','--output_fp',type="new_filepath",
     help='the output filepath')
    result['mapping_fp'] =\
     make_option('-m','--mapping_fp',type="existing_filepath",
     help='the mapping filepath')
     
    ## Define options used by the workflow scripts
    result['jobs_to_start_workflow'] =\
     make_option('-O','--jobs_to_start',type='int',
       help='Number of jobs to start. NOTE: you must also'
       ' pass -a to run in parallel, this defines the number of'
       ' jobs to be started if and only if -a is passed'
       ' [default: %default]',
       default=qiime_config['jobs_to_start'])
    
    ## Define options used by the parallel scripts
    result['jobs_to_start'] =\
     make_option('-O','--jobs_to_start',type='int',\
       help='Number of jobs to start [default: %default]',\
       default=qiime_config['jobs_to_start'])
    result['poller_fp'] =\
     make_option('-P','--poller_fp',action='store',\
       help='full path to '+\
       'qiime/parallel/poller.py [default: %default]',\
       default=join(get_qiime_scripts_dir(),'poller.py'))
    result['retain_temp_files'] =\
     make_option('-R','--retain_temp_files',action='store_true',\
       help='retain temporary files after runs complete '+\
       '(useful for debugging) [default: %default]',\
       default=False)
    result['suppress_submit_jobs'] =\
     make_option('-S','--suppress_submit_jobs',action='store_true',\
       help='Only split input and write commands file - don\'t submit '+\
       'jobs [default: %default]',default=False)
    result['poll_directly'] =\
     make_option('-T','--poll_directly',action='store_true',\
        help='Poll directly for job completion rather than running '+\
        'poller as a separate job. If -T is specified this script will '+\
        'not return until all jobs have completed. [default: %default]',\
        default=False)
    result['cluster_jobs_fp'] =\
     make_option('-U','--cluster_jobs_fp',
        help='path to cluster jobs script (defined in qiime_config) ' +\
        ' [default: %default]',\
        default=qiime_config['cluster_jobs_fp'] or\
         join(get_qiime_scripts_dir(),'start_parallel_jobs.py'))
    result['suppress_polling'] =\
     make_option('-W','--suppress_polling',action='store_true',
        help='suppress polling of jobs and merging of results '+\
        'upon completion [default: %default]',\
        default=False)
    result['job_prefix'] =\
     make_option('-X','--job_prefix',help='job prefix '+\
           '[default: descriptive prefix + random chars]')
    result['python_exe_fp'] =\
     make_option('-Y','--python_exe_fp',
        help='full path to python executable [default: %default]',\
        default=qiime_config['python_exe_fp'])
    result['seconds_to_sleep'] =\
     make_option('-Z','--seconds_to_sleep',type='int',\
        help='Number of seconds to sleep between checks for run '+\
        ' completion when polling runs [default: %default]',\
        default=qiime_config['seconds_to_sleep'] or 60)
     
    return result

def matrix_stats(headers_list, distmats):
    """does, mean, median, stdev on a series of (dis)similarity matrices
    
    takes a series of parsed matrices (list of headers, list of numpy 2d arrays)
    headers must are either row or colunm headers (those must be identical)
    outputs headers (list), means, medians, stdevs (all numpy 2d arrays)
    """
    
    if len(set(map(tuple,headers_list))) > 1:
        raise ValueError("error, not all input matrices have"+\
          " identical column/row headers")
        
    all_mats = numpy.array(distmats) # 3d numpy array: mtx, row, col
    means = numpy.mean(all_mats, axis=0)
    medians = numpy.median(all_mats, axis=0)
    stdevs = numpy.std(all_mats, axis=0)
    
    return deepcopy(headers_list[0]), means, medians, stdevs


def merge_otu_tables(otu_table_f1,otu_table_f2):
    """ Merge two otu tables with non-overlapping sample IDs
    
        WARNING: The OTU ids must refer to the same OTUs, which
         typically only happens when OTUs were picked against a 
         reference database, as with the BLAST OTU picker.
    
    """
    if isinstance(otu_table_f1, tuple):
        sample_ids1, otu_ids1, otu_table1, lineages1 = otu_table_f1
    else:
        sample_ids1, otu_ids1, otu_table1, lineages1 =\
            parse_otu_table(otu_table_f1)
    
    if isinstance(otu_table_f2, tuple):
        sample_ids2, otu_ids2, otu_table2, lineages2 = otu_table_f2
    else:
        sample_ids2, otu_ids2, otu_table2, lineages2 =\
            parse_otu_table(otu_table_f2)
    
    assert set(sample_ids1) & set(sample_ids2) == set(),\
     'Overlapping sample ids detected.'
    sample_ids_result = sample_ids1 + sample_ids2
    sample_ids_result_lookup = dict(
     [(sid,i) for i, sid in enumerate(sample_ids_result)])
    
    # Will need to add support for OTU tables wo tax info at some 
    # point -- in a rush now so don't have time to add it without an
    # immediate use case.
    if lineages1 and lineages2:    
        # map OTU ids to lineages -- in case of conflicts (i.e, OTU assigned)
        # different lineage in different otu tables, the lineage from 
        # OTU table 1 will be taken
        lineages = True
        otu_id_to_lineage = dict(zip(otu_ids1,lineages1))
        otu_id_to_lineage.update(dict([(otu_id,lineage)\
         for otu_id,lineage in zip(otu_ids2,lineages2)\
         if otu_id not in otu_id_to_lineage]))
    elif not (lineages1 or lineages2):
        lineages = False
    else:
      raise ValueError, ('Taxonomic information must be provided either'
       ' for all or none of the OTU tables')
    
    # Get the union of the otu IDs
    otu_ids_result = list(otu_ids1)
    otu_ids_lookup = {}.fromkeys(otu_ids1)
    otu_ids_result.extend([otu_id for otu_id in otu_ids2 \
                                  if otu_id not in otu_ids_lookup])
    otu_ids_result_lookup = dict(
     [(oid,i) for i, oid in enumerate(otu_ids_result)])
    
    otu_table = zeros(shape=(len(otu_ids_result),len(sample_ids_result)),dtype=int)
    for i,sample_id in enumerate(sample_ids1):
        #col_index = sample_ids_result.index(sample_id)
        col_index = sample_ids_result_lookup[sample_id]
        for j,otu_id in enumerate(otu_ids1):
            #row_index = otu_ids_result.index(otu_id)
            row_index = otu_ids_result_lookup[otu_id]
            otu_table[row_index,col_index] = otu_table1[j,i]
        
    for i,sample_id in enumerate(sample_ids2):
        #col_index = sample_ids_result.index(sample_id)
        col_index = sample_ids_result_lookup[sample_id]
        for j,otu_id in enumerate(otu_ids2):
            #row_index = otu_ids_result.index(otu_id)
            row_index = otu_ids_result_lookup[otu_id]
            otu_table[row_index,col_index] = otu_table2[j,i]
    
    if lineages:
        lineages_result = [otu_id_to_lineage[otu_id] 
         for otu_id in otu_ids_result]
    else:
        lineages_result = None
    
    return sample_ids_result, otu_ids_result, otu_table, lineages_result
    
def merge_n_otu_tables(otu_table_fs):
    """ Merge n otu tables """
    if len(otu_table_fs) < 2:
        raise ValueError, "Two or more OTU tables must be provided."
    otu_table_f0 = otu_table_fs[0]
    for otu_table_f in otu_table_fs[1:]:
        sample_names, otu_names, data, taxonomy = \
         merge_otu_tables(otu_table_f0,otu_table_f)
        otu_table_f0 = format_otu_table(sample_names=sample_names, 
                                    otu_names=otu_names,
                                    data=data,
                                    taxonomy=taxonomy).split('\n')
    
    return sample_names, otu_names, data, taxonomy

def convert_otu_table_relative(otu_table):
    """Convert the OTU table to relative abundances

    this method works on a parsed OTU table
    """
    sample_ids, otu_ids, otu_counts, consensus = otu_table
    otu_counts = asarray(otu_counts, float)
    otu_counts = otu_counts / otu_counts.sum(axis=0)
    otu_counts = where(isnan(otu_counts), 0.0, otu_counts)
    return (sample_ids, otu_ids, otu_counts, consensus)

def convert_OTU_table_relative_abundance(otu_table):
    """convert the OTU table to have relative abundances rather than raw counts
    """
    output = []
    data_lines = []
    otu_ids = []
    tax_strings = []
    taxonomy=False
    for line in otu_table:
        line = line.strip().split('\t')
        if line[0].startswith('#OTU ID'):
            output.append('\t'.join(line))
            if line[-1] == 'Consensus Lineage':
                taxonomy=True
        elif line[0].startswith('#'):
            output.append('\t'.join(line))
        else:
            if taxonomy:
                vals = [float(i) for i in line[1:-1]]
                tax_strings.append(line[-1])
            else:
                vals = [float(i) for i in line[1:]]
                tax_string = None
            data = array(vals, dtype=float)
            data_lines.append(data)
            otu_ids.append(line[0])
    data_lines = array(data_lines)
    totals = sum(data_lines)
    new_values = []
    for i in data_lines:
        new_values.append(i/totals)
    for index, i in enumerate(new_values):
        line = [otu_ids[index]]
        line.extend([str(j) for j in i])
        if taxonomy:
            line.append(tax_strings[index])
        output.append('\t'.join(line))
    return output




def load_pcoa_files(pcoa_dir):
    """loads PCoA files from filepaths
    """
    support_pcoas = []
    pcoa_filenames = os.listdir(pcoa_dir)
    #ignore invisible files like .DS_Store
    pcoa_filenames = [fname for fname in pcoa_filenames if not \
        fname.startswith('.')]
    master_pcoa = open(os.path.join(pcoa_dir, pcoa_filenames[0]), 'U')
    master_pcoa = parse_coords(master_pcoa)
    for fname in pcoa_filenames:
        try:
            f = open(os.path.join(pcoa_dir, fname), 'U')
            pcoa_res = parse_coords(f)
            support_pcoas.append(pcoa_res)
            f.close()
        except IOError, err:
            sys.sterr.write('error loading support pcoa ' + fname + '\n')
            exit(1)
    return master_pcoa, support_pcoas

def summarize_pcoas(master_pcoa, support_pcoas, method='IQR', apply_procrustes=True):
    """returns the average PCoA vector values for the support pcoas

    Also returns the ranges as calculated with the specified method. 
    The choices are:
        IQR: the Interquartile Range
        ideal fourths: Ideal fourths method as implemented in scipy
    """
    if apply_procrustes:
        # perform procrustes before averaging
        support_pcoas = [list(sp) for sp in support_pcoas]
        master_pcoa = list(master_pcoa)
        for i, pcoa in enumerate(support_pcoas):
            master_std, pcoa_std, m_squared = procrustes(master_pcoa[1],pcoa[1])
            support_pcoas[i][1] = pcoa_std
        master_pcoa[1] = master_std

    m_matrix = master_pcoa[1]
    m_eigvals = master_pcoa[2]
    m_names = master_pcoa[0]
    jn_flipped_matrices = []
    all_eigvals = []
    for rep in support_pcoas:
        matrix = rep[1]
        eigvals = rep[2]
        all_eigvals.append(eigvals)
        jn_flipped_matrices.append(_flip_vectors(matrix, m_matrix))
    matrix_average, matrix_low, matrix_high = _compute_jn_pcoa_avg_ranges(\
            jn_flipped_matrices, method)
    #compute average eigvals
    all_eigvals_stack = vstack(all_eigvals)
    eigval_sum = numpy.sum(all_eigvals_stack, axis=0)
    eigval_average = eigval_sum / float(len(all_eigvals))
    return matrix_average, matrix_low, matrix_high, eigval_average, m_names

def _compute_jn_pcoa_avg_ranges(jn_flipped_matrices, method):
    """Computes PCoA average and ranges for jackknife plotting

    returns 1) an array of jn_averages
             2) an array of upper values of the ranges
            3) an array of lower values for the ranges

    method: the method by which to calculate the range
        IQR: Interquartile Range
        ideal fourths: Ideal fourths method as implemented in scipy
    """
    x,y = shape(jn_flipped_matrices[0])
    all_flat_matrices = [matrix.ravel() for matrix in jn_flipped_matrices]
    summary_matrix = vstack(all_flat_matrices)
    matrix_sum = numpy.sum(summary_matrix, axis=0)
    matrix_average = matrix_sum / float(len(jn_flipped_matrices))
    matrix_average = matrix_average.reshape(x,y)
    if method == 'IQR':
        result = matrix_IQR(summary_matrix)
        matrix_low = result[0].reshape(x,y)
        matrix_high = result[1].reshape(x,y)
    elif method == 'ideal_fourths':
        result = idealfourths(summary_matrix, axis=0)
        matrix_low = result[0].reshape(x,y)
        matrix_high = result[1].reshape(x,y)
    elif method == "sdev":
        # calculate std error for each sample in each dimension
        sdevs = zeros(shape=[x,y])
        for j in xrange(y):
            for i in xrange(x):
                vals = array([pcoa[i][j] for pcoa in jn_flipped_matrices])
                sdevs[i,j] = vals.std(ddof=1)
        matrix_low = -sdevs/2
        matrix_high = sdevs/2


    return matrix_average, matrix_low, matrix_high

def _flip_vectors(jn_matrix, m_matrix):
    """transforms PCA vectors so that signs are correct"""
    m_matrix_trans = m_matrix.transpose()
    jn_matrix_trans = jn_matrix.transpose()
    new_matrix= zeros(jn_matrix_trans.shape, float)
    for i, m_vector in enumerate(m_matrix_trans):
        jn_vector = jn_matrix_trans[i]
        disT = list(m_vector - jn_vector)
        disT = sum(map(abs, disT))
        jn_flip = jn_vector*[-1]
        disF = list(m_vector - jn_flip)
        disF = sum(map(abs, disF))
        if disT > disF:
            new_matrix[i] = jn_flip
        else:
            new_matrix[i] = jn_vector
    return new_matrix.transpose()

def IQR(x):
    """calculates the interquartile range of x

    x can be a list or an array
    
    returns min_val and  max_val of the IQR"""

    x.sort()
    #split values into lower and upper portions at the median
    odd = len(x) % 2
    midpoint = int(len(x)/2)
    if odd:
        low_vals = x[:midpoint]
        high_vals = x[midpoint+1:]
    else: #if even
        low_vals = x[:midpoint]
        high_vals = x[midpoint:]
    #find the median of the low and high values
    min_val = median(low_vals)
    max_val = median(high_vals)
    return min_val, max_val

def matrix_IQR(x):
    """calculates the IQR for each column in an array
    """
    num_cols = x.shape[1]
    min_vals = zeros(num_cols)
    max_vals = zeros(num_cols)
    for i in range(x.shape[1]):
        col = x[:, i]
        min_vals[i], max_vals[i] = IQR(col)
    return min_vals, max_vals

def idealfourths(data, axis=None):
    """This function returns an estimate of the lower and upper quartiles of the data along
    the given axis, as computed with the ideal fourths. This function was taken
    from scipy.stats.mstat_extra.py (http://projects.scipy.org/scipy/browser/trunk/scipy/stats/mstats_extras.py?rev=6392)
    """
    def _idf(data):
        x = data.compressed()
        n = len(x)
        if n < 3:
            return [numpy.nan,numpy.nan]
        (j,h) = divmod(n/4. + 5/12.,1)
        qlo = (1-h)*x[j-1] + h*x[j]
        k = n - j
        qup = (1-h)*x[k] + h*x[k-1]
        return [qlo, qup]
    data = numpy.sort(data, axis=axis).view(MaskedArray)
    if (axis is None):
        return _idf(data)
    else:
        return apply_along_axis(_idf, axis, data)

def isarray(a):
    """
    This function tests whether an object is an array
    """
    try:
        validity=isinstance(a,ndarray)
    except:
        validity=False

    return validity

#make an alphabet that allows '.' as additional gaps
DNA_with_more_gaps = MolType(
    Sequence = DnaSequence,
    motifset = IUPAC_DNA_chars,
    Ambiguities = IUPAC_DNA_ambiguities,
    label = "dna",
    Gaps = ".",
    MWCalculator = DnaMW,
    Complements = IUPAC_DNA_ambiguities_complements,
    Pairs = DnaStandardPairs,
    make_alphabet_group=True,
    ModelSeq = ModelDnaSequence,
    )

def degap_fasta_aln(seqs):
    """degap a Fasta aligment.

    seqs: list of label,seq pairs
    """
    
    for (label,seq) in seqs:
        degapped_seq = Sequence(moltype=DNA_with_more_gaps,
                                seq=seq, name=label).degap()
        degapped_seq.Name = label
        yield degapped_seq

def write_degapped_fasta_to_file(seqs, tmp_dir="/tmp/"):
    """ write degapped seqs to temp fasta file."""

    tmp_filename = get_tmp_filename(tmp_dir=tmp_dir, prefix="degapped_", suffix=".fasta")
    fh = open(tmp_filename,"w")
    
    for seq in degap_fasta_aln(seqs):
        fh.write(seq.toFasta()+"\n")
    fh.close()
    return tmp_filename


def get_diff_for_otu_maps(otu_map1, otu_map2):
    """return reads in two otu_maps that are not shared

    otu_map1, otu_map2: OTU to seqID mapping as dict of lists
    """
 
    otus1 = set(otu_map1.keys())
    otus2 = set(otu_map2.keys())
    ids1 = set([x for otu in otus1 for x in otu_map1[otu]])
    ids2 = set([x for otu in otus2 for x in otu_map2[otu]])
        
    return ids1-ids2, ids2-ids1

def compare_otu_maps(otu_map1, otu_map2, verbose=False):
    """compare two otu maps and compute fraction of 

    otu_map1, otu_map2: OTU to seqID mapping as dict of lists
    """

    right = 0
    wrong = 0

    otus1 = set(otu_map1.keys())
    otus2 = set(otu_map2.keys())
    shared_otus = otus1.intersection(otus2)
    # check for equal members in shared OTUs
    for otu in shared_otus:
        members1 = set(otu_map1[otu])
        members2 = set(otu_map2[otu])
        
        right += len(members1 & members2)
        missing_in_2  = members1 - members2
        wrong += len(missing_in_2)
        if (verbose and len(missing_in_2)>0):
            print "OTU id: %s" % otu
            print list(missing_in_2)
            print 

    # process OTUs in 1 not in 2
    for otu in otus1 - shared_otus:
        wrong += len(otu_map1[otu])
        if verbose:
            print "OTU id: %s" % otu
            print list(otu_map1[otu])

    return float(wrong)/(right+wrong)
    
def compute_days_since_epoch(day,month,year):
    """ pass day, month, year to compute days since epoch (1/1/1970) 
        
        Note that full years should always be provided: 09 is a 
        different year than 2009!
    """
    d = datetime(int(year),int(month),int(day))
    epoch = datetime.utcfromtimestamp(0)
    return (d - epoch).days
    
def get_interesting_mapping_fields(mapping_data,mapping_headers):
    """ Returns headers for fields that are useful to color by in plots 
    
        These fields are the ones that contain greater than one value 
         and less values than the number of entries (so for example 
         not SampleID) 
    """
    result = []
    num_samples = len(mapping_data)
    num_cols = len(mapping_headers)
    transposed_data = array(mapping_data).T
    for header, datum in zip(mapping_headers, transposed_data):
        d = set(datum)
        len_d = len(d)
        if len_d > 1 and len_d < num_samples:
            result.append(header)
    return result
    
def flowgram_id_to_seq_id_map(seqs):
    """Map flowgram ids to sequence ids from a post-split_libraries fasta file
    """
    result = {}
    for id_, seq in seqs:
        fields = id_.split()
        seq_id = id_
        flowgram_id = fields[1]
        result[flowgram_id] = seq_id
    return result
    
def qiime_system_call(cmd):
    """ Call cmd and return (stdout, stderr, return_value)"""
    proc = Popen(cmd,shell=True,universal_newlines=True,\
                 stdout=PIPE,stderr=PIPE)
    # communicate pulls all stdout/stderr from the PIPEs to 
    # avoid blocking -- don't remove this line!
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return stdout, stderr, return_value
    
def get_qiime_library_version():
    """get QIIME version, including the svn version if applicable"""
    qiime_dir = get_qiime_project_dir()
    qiime_version = qiime_library_version
    cmd_call = 'svn info %s | egrep "Revision: "' % (qiime_dir)
    o,e,r = qiime_system_call(cmd_call)
    svn_revision = o.replace("Revision: ","").strip()
    
    if len(svn_revision) > 0:
        return '%s, svn revision %s' % (__version__, svn_revision)
    else:
        return '%s' % __version__

def inflate_denoiser_output(centroid_seqs,singleton_seqs,denoiser_map,raw_seqs):
    """Expand denoiser fasta files based on denoiser map
    
        The inflation process works as follows: write each centroid 
         sequence n times, where n is the number of reads in that 
         cluster, and write each singleton once. While writing these 
         out map back to original sequence identifiers.
         
        The seqs objects passed in are lists of (seq_id, seq) tuples,
         as returned from MinimalFastaParser.

    
    """
    id_lookup = parse_denoiser_mapping(denoiser_map)
    flowgram_to_seq_id_lookup = flowgram_id_to_seq_id_map(raw_seqs)
    for id_, seq in centroid_seqs:
        #centroid headers look like
        #>FZTHQMS01E140G | cluster size: 4353
        id, cluster_size_str = id_.split(' | ')
        cluster_member_ids = id_lookup[id]
        for c in cluster_member_ids:
            yield flowgram_to_seq_id_lookup[c], seq
    
    for id_, seq in singleton_seqs:
        yield flowgram_to_seq_id_lookup[id_], seq
    
    return
    
## Functions for counting sequences in fasta files

def count_seqs(fasta_filepath):
    """ Count the sequences in fasta_filepath 
    
        fasta_filepath: string indicating the full path to the file
    """
    # Open the file and pass it to py_count_seqs_from_file -- wrapping
    # this makes for easier unit testing
    return count_seqs_from_file(open(fasta_filepath))

def count_seqs_from_file(fasta_file):
    """Return number of sequences in fasta_file (no format checking performed)
    
        fasta_file: an open file object
    
    """
    result = 0
    lens = []
    for seq_id,seq in MinimalFastaParser(fasta_file):
        result += 1
        lens.append(len(seq))
    if result == 0:
        return result, None, None
    else:
        return result, mean(lens), std(lens)
        

def count_seqs_in_filepaths(fasta_filepaths,seq_counter=count_seqs):
    """ Wrapper to apply seq_counter to fasta_filepaths
    
        fasta_filepaths: list of one or more fasta filepaths
        seq_counter: a function which takes a single filepath 
         and returns the count of the number of sequences
         (default: count_seqs) -- this is parameterized to
         facilitate unit testing
    """
    total = 0
    counts = []
    inaccessible_filepaths = []
    # iterate over the input files
    for fasta_filepath in fasta_filepaths:
        try:
            # get the count of sequences in the current file
            current_count = seq_counter(fasta_filepath)
            # store it
            counts.append((current_count,fasta_filepath))
            # and increment the total count
            total += current_count[0]
        except IOError:
            # if the file couldn't be open, keep track of the filepath
            inaccessible_filepaths.append(fasta_filepath)
    
    return counts, total, inaccessible_filepaths
    
## End functions for counting sequences in fasta files

def get_top_fastq_two_lines(open_file):
    """ This function returns the first 4 lines of the open fastq file
    """
    line1 = open_file.readline()
    line2 = open_file.readline()
    line3 = open_file.readline()
    line4 = open_file.readline()
    open_file.seek(0)
    return line1, line2, line3, line4
    

def get_split_libraries_fastq_params_and_file_types(fastq_fps,mapping_fp):
    """ The function takes a list of open fastq files and a mapping file, then 
        returns a recommended parameters string for split_libraries_fastq
    """
    #parse the mapping
    data, headers, run_description= parse_mapping_file(open(mapping_fp,'U'))
    
    #determine the which column of mapping file is the BarcodeSequence
    for i,col_head in enumerate(headers):
        if col_head=='BarcodeSequence':
            barcode_column=i
            
    #create a set of barcodes for easier lookup
    barcode_mapping_column=set(zip(*data)[barcode_column])
    
    #create set of reverse complement barcodes from mapping file
    revcomp_barcode_mapping_column=[]
    for i in barcode_mapping_column:
        revcomp_barcode_mapping_column.append(DNA.rc(i))
    revcomp_barcode_mapping_column=set(revcomp_barcode_mapping_column)
    
    # get the filenames and sort them, so the file1 corresponds to file2
    fastq_fps.sort()
    
    # get the len of the sequence in each of the files, so we can determine
    # which file is the sequence file and which is the barcode sequence
    get_file_type_info={}
    for fastq_file in fastq_fps:
        fastq_fp=open(fastq_file)
        file_lines=get_top_fastq_two_lines(fastq_fp)
        parsed_fastq=MinimalFastqParser(file_lines,strict=False)
        for i,seq_data in enumerate(parsed_fastq):
            if i==0:
                get_file_type_info[fastq_file]=len(seq_data[1])
            else:
                break
        fastq_fp.close()

    # iterate over the sequence lengths and assign each file to either 
    # a sequence list or barcode list
    barcode_files=[]
    sequence_files=[]
    for i in range(0,len(fastq_fps),2):
        if get_file_type_info[fastq_fps[i]]<get_file_type_info[fastq_fps[i+1]]:
            barcode_files.append(fastq_fps[i])
            sequence_files.append(fastq_fps[i+1])
        else:
            barcode_files.append(fastq_fps[i+1])
            sequence_files.append(fastq_fps[i])
    
    # count the number of barcode matches in the forward and reverse direction
    # to determine if the rev_comp_barcode option needs passed
    fwd_count=0
    rev_count=0
    for bfile in barcode_files:
        fastq_fp=open(bfile)
        parsed_fastq=MinimalFastqParser(fastq_fp,strict=False)
        for bdata in parsed_fastq:
            if bdata[1] in barcode_mapping_column:
                fwd_count+=1
            elif bdata[1] in revcomp_barcode_mapping_column:
                rev_count+=1
        fastq_fp.close()
    
    # determine which barcode direction is correct
    if rev_count > fwd_count:
        barcode_orientation='--rev_comp_barcode'
    else:
        barcode_orientation=''
    
    #generate the string to use in command call to split_libraries_fastq
    split_lib_str='-i %s -b %s %s' % (','.join(sequence_files),
                                      ','.join(barcode_files),
                                      barcode_orientation)
    return split_lib_str


def iseq_to_qseq_fields(line,barcode_in_header,barcode_length,barcode_qual_c='b'):
    """ Split an Illumina sequence line into qseq fields"""
    record = line.strip().split(':')
    rec_0_1, rec_0_2 = record[0].split('_')
    rec_4_1, rec_4_23 = record[4].split('#')
    rec_4_2, rec_4_3 = rec_4_23.split('/')
    if barcode_in_header:   
        barcode = rec_4_2[:barcode_length]
        sequence = record[5]
        barcode_qual = barcode_qual_c*barcode_length
        sequence_qual = record[6]
    else:
        barcode = record[5][:barcode_length]
        sequence = record[5][barcode_length:]
        barcode_qual = record[6][:barcode_length]
        sequence_qual = record[6][barcode_length:]
    return (rec_0_1,rec_0_2,record[1],record[2],record[3],\
            rec_4_1,rec_4_2,rec_4_3), sequence, sequence_qual,\
            barcode,barcode_qual
            
def gzip_open(fp):
    return gzip.open(fp,'rb')
    
def make_compatible_distance_matrices(dm1,dm2,lookup=None):
    """ Intersect distance matrices and sort the values """
    dm1_ids = dm1[0]
    dm1_data = dm1[1]
    dm2_ids = dm2[0]
    dm2_data = dm2[1]
    if lookup:
        try:
            dm1_ids = [lookup[e] for e in dm1_ids]
            dm2_ids = [lookup[e] for e in dm2_ids]
        except KeyError,e:
            raise KeyError,\
             ("All entries in both DMs must be in "
              "lookup if a lookup is provided. Missing: %s" % str(e))
    order = [e for e in dm1_ids if e in dm2_ids]
    
    # create Dict2D from dm1
    d1 = {}
    for i,r in enumerate(dm1_ids):
        d1[r] = {}
        for j,c in enumerate(dm1_ids):
            d1[r][c] = dm1_data[i,j]
    result1 = Dict2D(data=d1,RowOrder=order,ColOrder=order)
    # remove entries not in order
    result1.purge()
    # return 2d list in order
    result1 = array(result1.toLists())
    
    # create Dict2D from dm2
    d2 = {}
    for i,r in enumerate(dm2_ids):
        d2[r] = {}
        for j,c in enumerate(dm2_ids):
            d2[r][c] = dm2_data[i,j]
    result2 = Dict2D(data=d2,RowOrder=order,ColOrder=order)
    # remove entries not in order
    result2.purge()
    # return 2d list in order
    result2 = array(result2.toLists())
    
    return (order,result1), (order,result2)
    
def get_rdp_jarpath():
    """ Return jar file name for RDP classifier ($RDP_JAR_PATH)"""
    return getenv('RDP_JAR_PATH')