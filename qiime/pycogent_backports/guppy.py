#!/usr/bin/env python
"""Application controller for guppy 1.1"""
       
__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Prototype"

from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, FilePath, system, \
       CommandLineAppResult, ResultPath, remove, ApplicationError
from cogent.core.alignment import Alignment
from os.path import splitext,split,join
from os import listdir
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode

class Guppy(CommandLineApplication):
    """guppy Application Controller
    """

    _command = 'guppy'
    _input_handler = '_input_as_multiline_string'
    _parameters = {
    #visualizations
        # makes trees with edges fattened in proportion to the number of reads
        'fat': FlagParameter('', Name='fat'),
        
        # maps an an arbitrary vector of the correct length to the tree
        'heat': FlagParameter('', Name='heat'),
        
        # writes a taxonomically annotated reference tree and an induced 
        # taxonomic tree
        'ref_tree': FlagParameter('', Name='ref_tree'),
        
        # makes one tree for each query sequence, showing uncertainty
        'sing': FlagParameter('', Name='sing'),
        
        # makes a tree with each of the reads represented as a pendant edge
        'tog': FlagParameter('', Name='tog'),
      
      #statistical comparison
        # draws the barycenter of a placement collection on the reference tree
        'bary': FlagParameter('', Name='bary'),
        
        # makes a phyloXML tree showing the bootstrap values
        'bootviz': FlagParameter('', Name='bootviz'),
        
        # calculates the EDPL uncertainty values for a collection of pqueries
        'edpl': FlagParameter('', Name='edpl'),
        
        # calculates the Kantorovich-Rubinstein distance and corresponding 
        # p-values
        'kr': FlagParameter('', Name='kr'),
        
        # makes a heat tree
        'kr_heat': FlagParameter('', Name='kr_heat'),
        
        # performs edge principal components
        'pca': FlagParameter('', Name='pca'),
        
        # writes out differences of masses for the splits of the tree
        'splitify': FlagParameter('', Name='splitify'),
        
        # performs squash clustering
        'squash': FlagParameter('', Name='squash'),
      
      #classification
        # outputs classification information in a tabular or SQLite format
        'classify': FlagParameter('', Name='classify'),
      
      #utilities
        # check a reference package
        'check_refpkg': FlagParameter('', Name='check_refpkg'),
        
        # splits apart placements with multiplicity, undoing a round procedure
        'demulti': FlagParameter('', Name='demulti'),
        
        # prints out a pairwise distance matrix between the edges
        'distmat': FlagParameter('', Name='distmat'),
        
        # filters one or more placefiles by placement name
        'filter': FlagParameter('', Name='filter'),
        
        # writes the number of leaves of the reference tree and the number of 
        # pqueries
        'info': FlagParameter('', Name='info'),
        
        # merges placefiles together
        'merge': FlagParameter('', Name='merge'),
        
        # restores duplicates to deduped placefiles
        'redup': FlagParameter('', Name='redup'),
        
        # clusters the placements by rounding branch lengths
        'round': FlagParameter('', Name='round'),
        
        # makes SQL enabling taxonomic querying of placement results
        'taxtable': FlagParameter('', Name='taxtable'),
        
        # converts old-style .place files to .json placement files
        'to_json': FlagParameter('', Name='to_json'),
        
        # Run the provided batch file of guppy commands
        'batch': FlagParameter('--', Name='batch'),
        
        # Print version and exit
        'version': FlagParameter('--', Name='version'),
        
        # Print a list of the available commands.
        'cmds': FlagParameter('--', Name='cmds'),
        
        # Display this list of options
        '--help': FlagParameter('--', Name='help'),
        
        # Display this list of options
        '-help': FlagParameter('-', Name='help'),
    }
 
    def getTmpFilename(self, tmp_dir='/tmp/',prefix='tmp',suffix='.json',\
           include_class_id=False,result_constructor=FilePath):
        """ Define Tmp filename to contain .json suffix, since guppy requires
            the suffix to be .json """
        
        return super(Guppy,self).getTmpFilename(tmp_dir=tmp_dir,
                                    prefix=prefix,
                                    suffix=suffix,
                                    include_class_id=include_class_id,
                                    result_constructor=result_constructor)
                                    
    def _handle_app_result_build_failure(self,out,err,exit_status,result_paths):
        """ Catch the error when files are not produced """
        raise ApplicationError, \
         'Guppy failed to produce an output file due to the following error: \n\n%s ' \
         % err.read()
    
    def _get_result_paths(self,data):
        basepath,basename=split(splitext(self._input_filename)[0])
        outfile_list=listdir(split(self._input_filename)[0])

        result = {}
        for i in outfile_list:
            if i.startswith(basename) and not i.endswith('.json') and \
                                                    not i.endswith('.txt'):
                result['result'] = ResultPath(Path=join(basepath,i))
                
        return result
    
def build_tree_from_json_using_params(fname,output_dir='/tmp/',params={}):
    """Returns a tree from a json.
    
    fname: filepath to input json 
    
    output_dir: location of output files
    
    params: dict of parameters to pass in to the RAxML app controller.

    The result will be a Tree.
    """

    # convert aln to fasta in case it is not already a fasta file
    
    ih = '_input_as_multiline_string'    

    guppy_app = Guppy(params=params,
                      InputHandler=ih,
                      WorkingDir=output_dir,
                      TmpDir=output_dir,
                      SuppressStderr=True,
                      SuppressStdout=True,
                      HALT_EXEC=False)
                      
    guppy_result = guppy_app(open(fname).read())
    
    try:
        new_tree=guppy_result['result'].read()
    except:
        # catch the error of not producing any results and print the command
        # run so user can check error
        guppy_cmd=Guppy(params=params,
                      InputHandler=ih,
                      WorkingDir=output_dir,
                      TmpDir=output_dir,
                      SuppressStderr=True,
                      SuppressStdout=True,
                      HALT_EXEC=True)
        out_msg=guppy_cmd(open(fname).read())
        
    tree = DndParser(new_tree, constructor=PhyloNode)
    
    guppy_result.cleanUp()

    return tree



    
