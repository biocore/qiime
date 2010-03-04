#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski"] 
__license__ = "GPL"
__version__ = "0.92"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"


"""Contains code for making a phylogeny from aligned sequences, using several techniques.

This module has the responsibility for taking a set of sequences and
returning an alignment. Mostly, it will be thin wrappers for code 
already in cogent.app.*, to which wrappers for e.g. NAST need to be
added..
"""

from cogent import LoadSeqs, DNA
from cogent.parse.fasta import MinimalFastaParser
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

    def getResult(self, aln_path, *args, **kwargs):
        """Returns alignment from sequences.
        
        Currently does not allow parameter tuning of program and uses
        default parameters -- this is bad and should be fixed.

        #TODO: allow command-line access to important aln params.
        """
        module = self.Params['Module']
        seqs = self.getAlignment(aln_path)
        result = module.build_tree_from_alignment(seqs, moltype=DNA)

        try: 
            root_method = kwargs['root_method']
            if root_method == 'midpoint':
                result = result.rootAtMidpoint()
            elif root_method == 'tree_method_default':
                pass
        except KeyError:
            pass
        return result

    def __call__(self, result_path=None, log_path=None, *args, **kwargs):
        """Calls superclass method to align seqs"""
        return FunctionWithParams.__call__(self, result_path=result_path,
                                           log_path=log_path, *args, **kwargs)


tree_method_constructors = {}  #currently unused
tree_module_names = {'muscle':cogent.app.muscle, 
    'clustalw':cogent.app.clustalw,
    #'mafft':cogent.app.mafft,   
    # current version of Mafft does not support tree building
    'fasttree':cogent.app.fasttree,'fasttree_v1':cogent.app.fasttree_v1,
    'raxml':cogent.app.raxml, 'clearcut':cogent.app.clearcut,
    }
