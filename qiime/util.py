#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight", "Daniel McDonald", "Greg Caporaso", 
"Justin Kuczynski","Jens Reeder","Catherine Lozupone"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"


"""Contains general utility code in support of the Qiime project.

A lot of this might migrate into cogent at some point.
"""

from StringIO import StringIO
from os import getenv, makedirs
#from scipy.stats.mstats import idealfourths
from os.path import abspath, exists, dirname, join, isdir
from numpy import min, max, median, mean
import numpy
from numpy.ma import MaskedArray
from numpy.ma.extras import apply_along_axis
from numpy import array, zeros, argsort, shape, vstack,ndarray
from collections import defaultdict
from optparse import OptionParser, OptionGroup, make_option
import sys
import os
from copy import deepcopy
from cogent import LoadSeqs
from cogent.cluster.procrustes import procrustes
from cogent.parse.tree import DndParser, PhyloNode
from cogent.core.alignment import Alignment
from cogent.app.blast import Blastall
from cogent.app.util import get_tmp_filename
from cogent.parse.blast import BlastResult
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.misc import remove_files
from cogent.util.dict2d import Dict2D
from cogent.app.formatdb import build_blast_db_from_fasta_path,\
    build_blast_db_from_fasta_file
from cogent import LoadSeqs

from cogent.util.misc import curry
from qiime.parse import parse_otu_table, parse_qiime_config_files, parse_coords

class TreeMissingError(IOError):
    """Exception for missing tree file"""
    pass

class OtuMissingError(IOError):
    """Exception for missing OTU file"""
    pass

class AlignmentMissingError(IOError):
    """Exception for missing alignment file"""
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
            tree = DndParser(f, PhyloNode)
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
    
    assert exists(result),\
     "qiime_scripts_dir does not exist: %s." % result +\
     " Have you defined it correctly in your qiime_config?"
    
    return result
    
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

def compute_seqs_per_library_stats(otu_f):
    counts = []
    sample_ids, otu_ids, otu_table, lineages = parse_otu_table(otu_f)
    for i in range(otu_table.shape[1]):
        counts.append(sum(otu_table[:,i]))
        
    return min(counts), max(counts), median(counts), mean(counts),\
     dict(zip(sample_ids,counts))
     
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
        
def sort_fasta_by_abundance(fasta_lines,fasta_out_f):
    """ Sort seqs in fasta_line by abundance, write all seqs to fasta_out_f
    
     Note that all sequences are written out, not just unique ones.
     
     fasta_lines: input file handle (or similar object)
     fasta_out_f: output file handle (or similar object)
     
    ** I am currently doing some work to figure out what the
     best way to do this is. The current implementation is going
     to have problems on very large (e.g., Illumina) files. --Greg
    
    """
    seq_index = {}
    count = 0
    for seq_id,seq in MinimalFastaParser(fasta_lines):
        count += 1
        try:
            seq_index[seq].append(seq_id)
        except KeyError:
            seq_index[seq] = [seq_id]
    
    seqs = []
    for k,v in seq_index.items():
        seqs.append((len(v),k,v))
        del seq_index[k]
    seqs.sort()
    for count, seq, seq_ids in seqs[::-1]:
        for seq_id in seq_ids:
            fasta_out_f.write('>%s\n%s\n' % (seq_id,seq))
    
    
    

## Begin functions for handling command line interfaces
def get_options_lookup():
    """ Return dict of commonly used options """
    qiime_config = load_qiime_config()
    result = {}
    result['fasta_as_primary_input'] =\
     make_option('-i','--input_fasta_fp',help='path to the input fasta file')
    result['otu_table_as_primary_input'] =\
     make_option('-i','--otu_table_fp',\
      help='path to the input OTU table (i.e., the output from make_otu_table.py)')
    result['otu_map_as_primary_input'] =\
     make_option('-i','--otu_map_fp',\
      help='path to the input OTU map (i.e., the output from pick_otus.py)')
    result['log_fp'] =\
     make_option('-l','--log_fp',help='path to write the log file')
    result['input_fasta'] =\
     make_option('-f','--input_fasta_fp',help='path to the input fasta file')
    result['output_dir'] =\
     make_option('-o','--output_dir',help='path to the output directory')
    result['output_fp'] =\
     make_option('-o','--output_fp',help='the output filepath')
    result['mapping_fp'] =\
     make_option('-m','--mapping_fp',help='the mapping filepath')
    
    ## Define options used by the parallel scripts
    result['jobs_to_start'] =\
     make_option('-O','--jobs_to_start',type='int',\
       help='Number of jobs to start [default: %default]',\
       default=qiime_config['jobs_to_start'])
    result['poller_fp'] =\
     make_option('-P','--poller_fp',action='store',\
       type='string',help='full path to '+\
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
        help='path to cluster_jobs.py script ' +\
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

def build_usage_lines(required_options,
    script_description,
    script_usage,
    optional_input_line,
    required_input_line):
    """ Build the usage string from components 
    """
    line1 = 'usage: %prog [options] ' + '{%s}' %\
     ' '.join(['%s %s' % (str(ro),ro.dest.upper())\
               for ro in required_options])
    usage_examples = []
    for title, description, command in script_usage:
        title = title.strip(':').strip()
        description = description.strip(':').strip()
        command = command.strip()
        if title:
            usage_examples.append('%s: %s\n %s' %\
             (title,description,command))
        else:
            usage_examples.append('%s\n %s' % (description,command))
    usage_examples = '\n\n'.join(usage_examples)
    lines = (line1,
             '', # Blank line
             optional_input_line,
             required_input_line,
             '', # Blank line
             script_description,
             '', # Blank line
             'Example usage: ',\
             'Print help message and exit',
             ' %prog -h\n',
             usage_examples)
    return '\n'.join(lines)

def set_parameter(key,kwargs,default=None):
    try:
        return kwargs[key]
    except KeyError:
        return default
        
def set_required_parameter(key,kwargs):
    try:
        return kwargs[key]
    except KeyError:
        raise KeyError,\
         "parse_command_line_parameters requires value for %s" % key
        
def parse_command_line_parameters(**kwargs):
    """ Constructs the OptionParser object and parses command line arguments
    
        parse_command_line_parameters takes a dict of objects via kwargs which
         it uses to build command line interfaces according to QIIME's 
         standards. The currently supported options are listed below with 
         their default values. If no default is provided, the option is 
         required.
        
        script_description
        script_usage
        version
        required_options=None
        optional_options=None
        suppress_verbose=False
        disallow_positional_arguments=True
        help_on_no_arguments=True
        optional_input_line = '[] indicates optional input (order unimportant)'
        required_input_line = '{} indicates required input (order unimportant)'
        
       These values can either be passed directly, as:
        parse_command_line_paramters(script_description="My script",\
                                     script_usage=[('Print help','%prog -h')],\
                                     version=1.0)
                                     
       or they can be passed via a pre-constructed dict, as:
        d = {'script_description'="My script",\
             'script_usage'=[('Print help','%prog -h')],\
             'version'=1.0}
        parse_command_line_paramters(**d)
    
    """
    # Get the options, or defaults if none were provided.
    script_description = set_required_parameter('script_description',kwargs)
    script_usage = set_required_parameter('script_usage',kwargs)
    version = set_required_parameter('version',kwargs)
    required_options = set_parameter('required_options',kwargs,[])
    optional_options = set_parameter('optional_options',kwargs,[])
    suppress_verbose = set_parameter('suppress_verbose',kwargs,False)
    disallow_positional_arguments =\
     set_parameter('disallow_positional_arguments',kwargs,True)
    help_on_no_arguments = set_parameter('help_on_no_arguments',kwargs,True)
    optional_input_line = set_parameter('optional_input_line',kwargs,\
        '[] indicates optional input (order unimportant)')
    required_input_line = set_parameter('required_input_line',kwargs,\
        '{} indicates required input (order unimportant)')
    
    # Build the usage and version strings
    usage = build_usage_lines(required_options,script_description,script_usage,\
                              optional_input_line,required_input_line)
    version = 'Version: %prog ' + version
    
    # Instantiate the command line parser object
    parser = OptionParser(usage=usage, version=version)
    
    # If no arguments were provided, print the help string (unless the
    # caller specified not to)
    if help_on_no_arguments and len(sys.argv) == 1:
        parser.print_usage()
        parser.exit()
    
    # Process the required options
    if required_options:
        # Define an option group so all required options are
        # grouped together, and under a common header
        required = OptionGroup(parser, "REQUIRED options",
         "The following options must be provided under all circumstances.")
        for ro in required_options:
            # if the option doesn't already end with [REQUIRED], 
            # add it.
            if not ro.help.strip().endswith('[REQUIRED]'):
                ro.help += ' [REQUIRED]'
            required.add_option(ro)
        parser.add_option_group(required)

    # Add a verbose parameter (if the caller didn't specify not to)
    if not suppress_verbose:
        parser.add_option('-v','--verbose',action='store_true',\
           dest='verbose',help='Print information during execution '+\
           '-- useful for debugging [default: %default]',default=False)

    # Add the optional options
    map(parser.add_option,optional_options)
    
    # Parse the command line
    opts,args = parser.parse_args()
    
    # If positional arguments are not allowed, and any were provided,
    # raise an error.
    if disallow_positional_arguments and len(args) != 0:
        parser.error("Positional argument detected: %s\n" % str(args[0]) +\
         " Be sure all parameters are identified by their option name.\n" +\
         " (e.g.: include the '-i' in '-i INPUT_DIR')")

    # Test that all required options were provided.
    if required_options:
        required_option_ids = [o.dest for o in required.option_list]
        for required_option_id in required_option_ids:
            if getattr(opts,required_option_id) == None:
                parser.error('Required option --%s omitted.' \
                             % required_option_id)
            
    # Return the parser, the options, and the arguments. The parser is returned
    # so users have access to any additional functionality they may want at 
    # this stage -- most commonly, it will be used for doing custom tests of 
    # parameter values.
    return parser, opts, args

## End functions for handling command line interfaces

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
    """ Merge two otu tables with the same sample IDs
    
        WARNING: The OTU ids must refer to the same OTUs, which
         typically only happens when OTUs were picked against a 
         reference database, as with the BLAST OTU picker.
    
    """
    sample_ids1, otu_ids1, otu_table1, lineages1 =\
        parse_otu_table(otu_table_f1)
    sample_ids2, otu_ids2, otu_table2, lineages2 =\
        parse_otu_table(otu_table_f2)
    
    assert set(sample_ids1) & set(sample_ids2) == set(),\
     'Overlapping sample ids detected.'
    sample_ids_result = sample_ids1 + sample_ids2
    
    # map OTU ids to lineages -- in case of conflicts (i.e, OTU assigned)
    # different lineage in different otu tables, the lineage from 
    # OTU table 1 will be taken
    otu_id_to_lineage = dict(zip(otu_ids1,lineages1))
    otu_id_to_lineage.update(dict([(otu_id,lineage)\
     for otu_id,lineage in zip(otu_ids2,lineages2)\
     if otu_id not in otu_id_to_lineage]))
    
    # Get the union of the otu IDs
    otu_ids_result = list(otu_ids1)
    otu_ids_lookup = {}.fromkeys(otu_ids1)
    otu_ids_result.extend([otu_id for otu_id in otu_ids2 \
                                  if otu_id not in otu_ids_lookup])
    
    otu_table = zeros(shape=(len(otu_ids_result),len(sample_ids_result)),dtype=int)
    for i,sample_id in enumerate(sample_ids1):
        col_index = sample_ids_result.index(sample_id)
        for j,otu_id in enumerate(otu_ids1):
            row_index = otu_ids_result.index(otu_id)
            otu_table[row_index,col_index] = otu_table1[j,i]
        
    for i,sample_id in enumerate(sample_ids2):
        col_index = sample_ids_result.index(sample_id)
        for j,otu_id in enumerate(otu_ids2):
            row_index = otu_ids_result.index(otu_id)
            otu_table[row_index,col_index] = otu_table2[j,i]

    lineages_result = [otu_id_to_lineage[otu_id] for otu_id in otu_ids_result]
    
    return sample_ids_result, otu_ids_result, otu_table, lineages_result

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


#some error codes for creating a dir
NO_ERROR       = 0
DIR_EXISTS     = 1
FILE_EXISTS    = 2
OTHER_OS_ERROR = 3

#Note: We temporarily set the default behaviour of create_dir to overwrite
#      without Error (fail_on_exist=False).
#      After the 1.0 release we will set his back
def create_dir(dir_name, fail_on_exist=False, handle_errors_externally=False):
    """Create a dir safely and fail meaningful.

    dir_name: name of directory to create

    fail_on_exist: if true raise an error if dir already exists
    
    handle_errors_externally: if True do not raise Errors, but return
                   failure codes. This allows to handle errors locally and
                   e.g. hint the user at a --force_overwrite options.
                   
    returns values (if no Error raised):
    
         0:  dir could be safely made
         1:  directory already existed
         2:  a file with the same name exists          
         3:  any other unspecified OSError


    See qiime/denoiser.py for an example of how to use this mechanism.

    Note: Depending  of how thorough we want to be we could add tests,
          e.g. for testing actual write permission in an existing dir.
    """

    #pre-instanciate function with
    ror = curry(handle_error_codes, dir_name, handle_errors_externally)

    if exists(dir_name):
        if isdir(dir_name):
            #dir is there
            if fail_on_exist:
                return ror(DIR_EXISTS)
            else:
                return DIR_EXISTS
        else:
            #must be file with same name
            return ror(FILE_EXISTS)
    else:
        #no dir there, try making it
        try:
            makedirs(dir_name)
        except OSError:
            return ror(OTHER_OS_ERROR)
    
    return NO_ERROR

def handle_error_codes(dir_name, supress_errors=False,
                       error_code=NO_ERROR):
    """Wrapper function for error_handling.

    dir_name: name of directory that raised the error
    suppress_errors: if True raise Errors, otherwise return error_code
    error_code: the code for the error
    """
    error_strings = \
        {DIR_EXISTS :"Directory already exists: %s" % dir_name,
         FILE_EXISTS : "File with same name exists: %s" % dir_name,
         OTHER_OS_ERROR: "Could not create output directory: %s. " % dir_name
         +"Check the permissions."}

    if error_code==NO_ERROR:
        return NO_ERROR
    if supress_errors:
        return error_code
    else:
        raise OSError, error_strings[error_code]

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
    """Returns an estimate of the lower and upper quartiles of the data along
    the given axis, as computed with the ideal fourths.
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
    Test for arrayobjects.
    """
    try:
        result=isinstance(a,ndarray)
    except:
        result=False

    return result
