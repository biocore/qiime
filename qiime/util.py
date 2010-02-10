#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight", "Daniel McDonald", "Greg Caporaso", 
"Justin Kuczynski"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"


"""Contains general utility code in support of the Qiime project.

A lot of this might migrate into cogent at some point.
"""

from StringIO import StringIO
from os import getenv
from os.path import abspath, exists, dirname
from numpy import min, max, median, mean
import numpy
from collections import defaultdict
from optparse import OptionParser, OptionGroup
import sys
from qiime.parse import parse_otus
from cogent import LoadSeqs
from cogent.parse.tree import DndParser, PhyloNode
from cogent.core.alignment import Alignment
from cogent.app.blast import Blastall
from cogent.app.util import get_tmp_filename
from cogent.parse.blast import BlastResult
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.misc import remove_files
from cogent.app.formatdb import build_blast_db_from_fasta_path
from cogent import LoadSeqs
from copy import deepcopy

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
                return parse_otus(StringIO(otu_source))
            except (TypeError, ValueError), e:
                raise OtuMissingError, \
                    "Tried to read OTUs from string starting with # but got "+e
        else:
            try:
                otu_file = open(otu_source, 'U')
            except (TypeError, IOError):
                raise OtuMissingError, \
                    "Couldn't read OTU file at path: %s" % otu_source
            result = parse_otus(otu_file)
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
    
        Although this value is stored in qiime_config, we are currently
         not requiring users to setup a qiime_config. Because developers
         were doing variants of this functionality throughout the code,
         it's safer to define this as a single tested function, and if we
         utimately do require users to set up a qiime_config, we can 
         modify this function to pull the value from there.
    
    """
    # Get the full path of util.py
    current_file_path = abspath(__file__)
    # Get the directory containing util.py
    current_dir_path = dirname(current_file_path)
    # Return the directory containing the directory containing util.py
    return dirname(current_dir_path)

# Begin functions for handling qiime_config file
def parse_qiime_config_file(qiime_config_file):
    """ Parse lines in a qiime_config file
    """
    result = {}
    for line in qiime_config_file:
        line = line.strip()
        # ignore blank lines or lines beginning with '#'
        if not line or line.startswith('#'): continue
        fields = line.split('\t')
        param_id = fields[0]
        param_value = '\t'.join(fields[1:]) or None
        result[param_id] = param_value
    return result

def parse_qiime_config_files(qiime_config_files):
    """ Parse files in (ordered!) list of qiime_config_files
    
        The order of files must be least important to most important.
         Values defined in earlier files will be overwritten if the same 
         values are defined in later files.
    """
    # The qiime_config object is a default dict: if keys are not
    # present, none is returned
    def return_none():
        return None
    results = defaultdict(return_none)
    
    for qiime_config_file in qiime_config_files:
        try:
            results.update(parse_qiime_config_file(qiime_config_file))
        except IOError:
            pass
    
    return results

def load_qiime_config():
    """Return default parameters read in from file"""
    
    qiime_config_filepaths = []
    qiime_project_dir = get_qiime_project_dir()
    qiime_config_filepaths.append(qiime_project_dir + '/qiime_config')
    
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
    
# End functions for handling qiime_config file

# This function is copied from the development branch of PyCogent, and should
# be replaced when we upgrade the requirement to PyCogent1.5.
def build_blast_db_from_fasta_file(fasta_file,is_protein=False,\
    output_dir=None,HALT_EXEC=False):
    """Build blast db from fasta_path; return db name and list of files created
    
        **If using to create temporary blast databases, you can call
        cogent.util.misc.remove_files(db_filepaths) to clean up all the
        files created by formatdb when you're done with the database.
    
        fasta_path: path to fasta file of sequences to build database from
        is_protein: True if working on protein seqs (default: False)
        output_dir: directory where output should be written
         (default: directory containing fasta_path)
        HALT_EXEC: halt just before running the formatdb command and
         print the command -- useful for debugging
    """
    output_dir = output_dir or '.'
    fasta_path = get_tmp_filename(\
     tmp_dir=output_dir, prefix="BLAST_temp_db_", suffix=".fasta")
    
    fasta_f = open(fasta_path,'w')
    for line in fasta_file:
        fasta_f.write('%s\n' % line.strip())
    fasta_f.close()
    
    blast_db, db_filepaths = build_blast_db_from_fasta_path(\
     fasta_path,is_protein=False,output_dir=None,HALT_EXEC=False)
     
    db_filepaths.append(fasta_path)
    
    return blast_db, db_filepaths

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
    sample_ids, otu_ids, otu_table, lineages = parse_otus(otu_f)
    for i in range(otu_table.shape[1]):
        counts.append(sum(otu_table[:,i]))
        
    return min(counts), max(counts), median(counts), mean(counts),\
     dict(zip(sample_ids,counts))

## Begin functions for handling command line interfaces
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
    lines = (line1,
             '', # Blank line
             optional_input_line,
             required_input_line,
             '', # Blank line
             script_description,
             '', # Blank line
             'Example usage:',
             ' Print help message and exit:',
             '  %prog -h',
             script_usage)
    return '\n'.join(lines)

def parse_command_line_parameters(script_description,\
    script_usage,\
    version,\
    required_options=None,
    optional_options=None,
    suppress_verbose=False,
    disallow_positional_arguments=True,
    help_on_no_arguments=True,
    optional_input_line = '[] indicates optional input (order unimportant)',
    required_input_line = '{} indicates required input (order unimportant)'):
    """ Constructs the OptionParser object and parses command line arguments
    """
    # Get the options, or empty lists if none were provided.
    required_options = required_options or []
    optional_options = optional_options or []
    
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
