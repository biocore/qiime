#!/usr/bin/env python

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight", "Daniel McDonald", "Greg Caporaso"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"


"""Contains general utility code in support of the Qiime project.

A lot of this might migrate into cogent at some point.
"""

from StringIO import StringIO
from os import getenv
from os.path import split
from numpy import min, max, median, mean
from qiime.parse import parse_otus
from cogent.parse.tree import DndParser, PhyloNode
from cogent.core.alignment import Alignment
from cogent.app.blast import Blastall
from cogent.parse.blast import BlastResult
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.misc import remove_files
from cogent.app.formatdb import build_blast_db_from_fasta_path,\
    build_blast_db_from_fasta_file
from cogent import LoadSeqs

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

def parse_qiime_config_filepaths(qiime_config_filepaths):
    """ Parse files in (ordered!) list of qiime_config_filepaths
    
        The order of file paths must be least important to most important.
         Values defined in earlier filepaths will be overwritten if the same 
         values are defined in later filepaths.
    """
    files_read = []
    
    results = {}
    for qiime_config_filepath in qiime_config_filepaths:
        try:
            results.update(parse_qiime_config_file(open(qiime_config_filepath)))
            files_read.append(qiime_config_filepath)
        except IOError:
            pass
            
    # If no qiime_config file could be read, raise an error
    if not files_read:
        raise IOError,\
         "No qiime_config were read.\n Tried %s\n Do these exist? Do you have read access?"\
         % ' '.join(qiime_config_filepaths)
         
    return results

def load_qiime_config():
    """Return default parameters read in from file"""
    result = {}
    qiime_config_filepaths = []
    qiime_parallel_dir = split(__file__)[0]
    qiime_config_filepaths.append(qiime_parallel_dir + '/../qiime_config')
    
    qiime_config_env_filepath = getenv('QIIME_CONFIG_FP')
    if qiime_config_env_filepath:
        qiime_config_filepaths.append(qiime_config_env_filepath)
    
    home_dir = getenv('HOME')
    if home_dir:
        qiime_config_home_filepath = home_dir + '/.qiime_config'
        qiime_config_filepaths.append(qiime_config_home_filepath)
    return parse_qiime_config_filepaths(qiime_config_filepaths)
    
qiime_config = load_qiime_config()
# End functions for handling qiime_config file

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