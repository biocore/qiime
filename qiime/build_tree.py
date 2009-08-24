#!/usr/bin/env python

__author__ = "YOUR NAME HERE"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "YOUR NAME HERE"
__email__ = "YOUR EMAIL ADDRESS HERE"
__status__ = "Prototype"


"""Contains code for aligning sequences, using several techniques.

This module has the responsibility for taking a set of sequences and
returning an alignment. Mostly, it will be thin wrappers for code 
already in cogent.app.*, to which wrappers for e.g. NAST need to be
added..
"""
from os import remove
from commands import getoutput
from optparse import OptionParser
from cogent import LoadSeqs, DNA
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.util import get_tmp_filename
from qiime.util import FunctionWithParams
#app controllers that implement align_unaligned_seqs
import cogent.app.muscle
import cogent.app.clustalw
import cogent.app.mafft
import cogent.app.raxml
import cogent.app.fasttree
import cogent.app.fasttree_v1
import cogent.app.clearcut

class TreeBuilder(FunctionWithParams):
    """A TreeBuilder takes a aligned set of sequences and returns a tree.

    This is an abstract class: subclasses should implement the __call__
    method.

    Note: sequence ids should be preserved during this process, i.e. the
    description lines should be saved/restored if the phylogeny app is
    destructive to them.

    Specific wrappers need to be written for nontraditional approaches,
    including:

    (a) RAxML assignment, where sequences are assigned to internal nodes of
        an existing tree
    (b) BLAST-based assignment, where sequences are assigned to existing
        nodes of a tree based on best blast hit (assigned only to terminal
        nodes if using a single hit, but could easily imagine assigning to
        an internal node based on a set of indistinguishably good hits).
    (c) BLAST-like approach using VMATCH or other suffix array library, or 
        using oligonucleotide freqs like the RDP classifier does to assign
        to an arbitrary internal node in an existing tree, etc.
    """
    Name = 'TreeBuilder'

    def __init__(self, params):
        """Return new TreeBuilder object with specified params.
        
        Note: expect params to contain both generic and per-method (e.g. for
        raxml vs. fasttree vs. whatever) params, so leaving it as a dict 
        rather than setting attributes. Some standard entries in params are:

        Application: 3rd-party application used, if any, e.g. raxml
        """
        self.Params = params

    def __call__ (self, aln_path, result_path=None, log_path=None):
        """Returns tree from alignment.
        
        Parameters:
        aln_path: path to file of aligned sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path as fasta, otherwise should
        return cogent.core.alignment.DenseAlignment object.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError, "TreeBuilder is an abstract class"

class CogentTreeBuilder(TreeBuilder):
    """Generic tree builder using Cogent tree methods."""

    Name = 'CogentTreeBuilder'

    def getResult(self, aln_path):
        """Returns alignment from sequences.
        
        Currently does not allow parameter tuning of program and uses
        default parameters -- this is bad and should be fixed.

        #TODO: allow command-line access to important aln params.
        """
        module = self.Params['Module']
        seqs = self.getAlignment(aln_path)
        result = module.build_tree_from_alignment(seqs, moltype=DNA)    
        
        return result

    def __call__(self, result_path=None, log_path=None, *args, **kwargs):
        """Calls superclass method to align seqs"""
        return FunctionWithParams.__call__(self, result_path=result_path,
            log_path=log_path, *args, **kwargs)

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage =\
     'usage: %prog [options] input_alignment_filepath'
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-t','--tree_method',action='store',\
          type='string',dest='tree_method',help='Method for tree building'+\
          ' [default: %default]')
          
    parser.add_option('-o','--result_fp',action='store',\
          type='string',dest='result_fp',help='Path to store '+\
          'result file [default: <input_sequences_filename>_aligned.fasta]')
          
    parser.add_option('-l','--log_fp',action='store',\
          type='string',dest='log_fp',help='Path to store '+\
          'log file [default: No log file created.]')
          
    parser.set_defaults(tree_method='fasttree')

    opts,args = parser.parse_args()
    num_args = 1
    if len(args) != num_args:
       parser.error('Exactly one argument is required.')
       
    if not (opts.tree_method in tree_method_constructors or
            opts.tree_method in tree_module_names):
        parser.error(\
         'Invalid alignment method: %s.\nValid choices are: %s'\
         % (opts.tree_method,\
            ' '.join(tree_method_constructors.keys() +
                tree_module_names.keys())))
            
    return opts,args

tree_method_constructors = {}  #currently unused
tree_module_names = {'muscle':cogent.app.muscle, 
    'clustalw':cogent.app.clustalw, 'mafft':cogent.app.mafft,
    'fasttree':cogent.app.fasttree,'fasttree_v1':cogent.app.fasttree_v1,
    'raxml':cogent.app.raxml, 'clearcut':cogent.app.clearcut,
    }

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    #verbose = opts.verbose
 
    try:
        tree_builder_constructor =\
        tree_method_constructors[opts.tree_method]
        tree_builder_type = 'Constructor'
        params = {}
        tree_builder = tree_builder_constructor(params)
    except KeyError:
        tree_builder = CogentTreeBuilder({
        'Module':tree_module_names[opts.tree_method],
        'Method':opts.tree_method
        })
        tree_builder_type = 'Cogent'
     
    input_seqs_filepath = args[0]
   
    result_path = opts.result_fp or\
     input_seqs_filepath.replace('.fasta','_aligned.fasta')
     
    log_path = opts.log_fp

    if tree_builder_type=='Constructor':
        tree_builder(input_seqs_filepath,\
        result_path=result_path,log_path=log_path,failure_path=failure_path)
    elif tree_builder_type=='Cogent':
        tree_builder(result_path, aln_path=input_seqs_filepath,
            log_path=log_path)

