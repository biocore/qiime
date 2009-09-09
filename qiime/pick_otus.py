#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight","Greg Caporaso", "Kyle Bittinger"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

"""Contains code for OTU picking, using several techniques.

This module has the responsibility for taking a set of sequences and
grouping those sequences by similarity.
"""

from copy import copy
from optparse import OptionParser
from cogent.app.cd_hit import cdhit_clusters_from_seqs
from cogent.app.dotur import dotur_from_alignment
from cogent.core.sequence import DnaSequence
from cogent import LoadSeqs, DNA, Alignment
from qiime.util import FunctionWithParams


class OtuPicker(FunctionWithParams):
    """An OtuPicker dereplicates a set of sequences at a given similarity.

    This is an abstract class: subclasses should implement the __call__
    method.
    """
    
    Name = 'OtuPicker'

    def __init__(self, params):
        """Return new OtuPicker object with specified params.
        
        Note: expect params to contain both generic and per-method (e.g. for
        cdhit) params, so leaving it as a dict rather than setting
        attributes. Some standard entries in params are:

        Algorithm: algorithm used (e.g. nearest-neighbor, furthest-neighbor)
        Similarity: similarity threshold, e.g. 0.97
        Application: 3rd-party application used, if any, e.g. cdhit
        """
        self.Params = params

    def __call__ (self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError, "OtuPicker is an abstract class"


class CdHitOtuPicker(OtuPicker):
    
    Name = 'CdHitOtuPicker'
    # Application = 'cd-hit-est'
    # Algorithm = 'cdhit: "longest-sequence-first list removal algorithm"'
    # Params = {'Similarity':0.96}
    
    def __init__(self, params):
        """Return new OtuPicker object with specified params.
        
        params contains both generic and per-method (e.g. for
        cdhit application controller) params.
        
        Some generic entries in params are:
    
        Algorithm: algorithm used
        Similarity: similarity threshold, default 0.96, corresponding to
         genus-level OTUs ('Similarity' is a synonym for the '-c' parameter
         to the cd-hit application controllers)
        Application: 3rd-party application used
        """
        _params = {'Similarity':0.96,\
         'Application':'cdhit',\
         'Algorithm':'cdhit: "longest-sequence-first list removal algorithm"'}
        _params.update(params)
        OtuPicker.__init__(self, _params)
    
    def __call__ (self, seq_path, result_path=None, log_path=None, id_len=0):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified,
        dumps the result to the desired path instead of returning it.
        log_path: path to log, which includes dump of params.
        id_len: if set, truncates ids to n chars (you don't want this!)

        """
        moltype = DNA
        
        # Load the seq path. Right now, cdhit_clusters_from_seqs
        # doesn't support being passed a file path even though the 
        # seqs do get written to a fasta file before being passed
        # to cd-hit-est. We may want to change that in the future 
        # to avoid the overhead of loading large sequence collections
        # during this step. 
        seqs = LoadSeqs(seq_path,moltype=moltype,aligned=False)
        
        # create the params dict to pass to cd-hit-est -- IS THERE A
        # BETTER WAY TO MAKE self.Params INTO THE params DICT TO PASS
        # TO THE APPLICATION CONTROLLERS?
        cd_hit_params = copy(self.Params)
        del cd_hit_params['Application']
        del cd_hit_params['Algorithm']
        cd_hit_params['-d'] = id_len  #turn off id truncation
        
        # Get the clusters by running cd-hit-est against the
        # sequence collection
        clusters = cdhit_clusters_from_seqs(\
         seqs=seqs,moltype=moltype,params=cd_hit_params)

        if result_path:
            # if the user provided a result_path, write the 
            # results to file with one tab-separated line per 
            # cluster
            of = open(result_path,'w')
            for i,cluster in enumerate(clusters):
                of.write('%s\t%s\n' % (i,'\t'.join(cluster)))
            of.close()
            result = None
            log_str = 'Result path: %s' % result_path
        else:
            # if the user did not provide a result_path, store
                # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = dict(enumerate(clusters))
            log_str = 'Result path: None, returned as dict.'
 
        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path,'w')
            log_file.write(str(self))
            log_file.write('\n')
            log_file.write('%s\n' % log_str)
    
        # return the result (note this is None if the data was
        # written to file)
        return result


class DoturOtuPicker(OtuPicker):
    Name = 'DoturOtuPicker'

    def __init__(self, params):
        """Return new DoturOtuPicker object with specified params.

        The Behavior of a DoturOtuPicker is controlled by the
        following parameters:

        Similarity: similarity threshold corresponding to genus-level
        OTUs (default 0.96)

        Distance Function: a function which takes two sequences and
        returns a distance ranging between 0 and 1 (default
        DnaSequence.fracDiff())

        Moltype: type of molecule (default: DNA)
        """
        params['Application'] = 'dotur'
        params['Algorithm'] = 'Furthest-neighbor clustering algorithm'

        if 'Similarity' not in params:
            params['Similarity'] = 0.96
        if 'Distance Function' not in params:
            def default_distance_function(first, second):
                first = DnaSequence(first)
                return first.fracDiff(second)
            params['Distance Function'] = default_distance_function 
        if 'Moltype' not in params:
            params['Moltype'] = DNA

        OtuPicker.__init__(self, params)

    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {otu_id:[seq_ids]} for each otu.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified,
        dumps the result to the desired path instead of returning it.
        log_path: path to log, which includes dump of params.
        """
        moltype = self.Params['Moltype']
        distance_function = self.Params['Distance Function']
        similarity_threshold = self.Params['Similarity']

        seqs = LoadSeqs(seq_path, moltype=moltype, aligned=True)
        parsed_results = dotur_from_alignment(seqs, moltype, distance_function)
        clusters = self.__pick_clusters(parsed_results, similarity_threshold)

        # From here down, this is all copied straight from
        # CdHitOtuPicker, and prime for refactoring into a private
        # method of OtuPicker

        if result_path:
            # if the user provided a result_path, write the 
            # results to file with one tab-separated line per 
            # cluster
            of = open(result_path,'w')
            for i,cluster in enumerate(clusters):
                of.write('%s\t%s\n' % (i,'\t'.join(cluster)))
            of.close()
            result = None
            log_str = 'Result path: %s' % result_path
        else:
            # if the user did not provide a result_path, store
            # the clusters in a dict of {otu_id:[seq_ids]}, where
            # otu_id is arbitrary
            result = dict(enumerate(clusters))
            log_str = 'Result path: None, returned as dict.'
 
        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path,'w')
            log_file.write(str(self))
            log_file.write('\n')
            log_file.write('%s\n' % log_str)
    
        # return the result (note this is None if the data was
        # written to file)
        return result

    @staticmethod
    def __pick_clusters(dotur_results, similarity_threshold):
        """Returns OTU's that satisfy the given similarity threshold.

        A higher similarity threshold results in a larger number of
        distinct OTU's.

        Since DOTUR returns all sets of OTU's at various degrees of
        similarity, this process is sufficiently involved to warrant
        its own static method.  This method should probably live with
        the DOTUR application controller in PyCogent, and has been set
        up for easy refactoring.
        """
        # In dotur, a lower similarity score means more otu's.  To
        # find otu's that satisfy a similarity threshold of 0.9, we
        # must find the largest dotur score less than (1 - 0.9 =) 0.1.
        dotur_score_threshold = 1 - similarity_threshold

        # info stored in various positions of the results list
        i_score = 0
        i_otu_list = 2

        # 'unique' in first result set indicates score of 0.0
        # change to numerical value
        dotur_results[0][i_score] = 0.0

        prev_res = dotur_results[0]
        for res in dotur_results:
            score = res[i_score]
            if score <= dotur_score_threshold:
                prev_res = res
            else:
                break
        
        # return only the list of otu's
        retval = prev_res[i_otu_list]

        return retval


def parse_command_line_parameters():
    """ Parses command line arguments """
    usage =\
     'usage: %prog [options] input_sequences_filepath'
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    # 
    # parser.add_option('-v','--verbose',action='store_true',\
    #     dest='verbose',help='Print information during execution -- '+\
    #     'useful for debugging [default: %default]')

    parser.add_option('-m','--otu_picking_method',action='store',\
          type='string',dest='otu_picking_method',help='Method for picking'+\
          ' OTUs [default: %default]')
    
    parser.add_option('-M','--max_cdhit_memory',action='store',\
          type=int,help='max available memory to cdhit (cd-hit\'s -M)'+\
          ' (Mbyte) [default: %default]',default=400)
          
    parser.add_option('-o','--result_fp',action='store',\
          type='string',dest='result_fp',help='Path to store '+\
          'result file [default: <input_sequences_filepath>.otu]')
          
    parser.add_option('-l','--log_fp',action='store',\
          type='string',dest='log_fp',help='Path to store '+\
          'log file [default: No log file created.]')
          
    parser.add_option('-s','--similarity',action='store',\
          type='float',dest='similarity',help='Sequence similarity '+\
          'threshold [default: %default]')

    parser.set_defaults(verbose=False,otu_picking_method='cdhit',\
        similarity=0.96)

    opts,args = parser.parse_args()
    num_args = 1
    if len(args) != num_args:
       parser.error('Exactly one argument is required.')
       
    if opts.otu_picking_method not in otu_picking_method_constructors:
        parser.error(\
         'Invalid assignment method: %s.\nValid choices are: %s'\
         % (opts.otu_picking_method,\
            ' '.join(otu_picking_method_constructors.keys())))

    return opts,args

otu_picking_method_constructors = dict([\
 ('cdhit',CdHitOtuPicker)])#,('dotur',DoturOtuPicker)])

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    #verbose = opts.verbose
 
    otu_picker_constructor =\
     otu_picking_method_constructors[opts.otu_picking_method]
     
    input_seqs_filepath = args[0]
   
    result_path = opts.result_fp or\
     '%s.otu' % input_seqs_filepath
     
    log_path = opts.log_fp
    
    params = {'Similarity':opts.similarity,'-M':opts.max_cdhit_memory}
    
    otu_picker = otu_picker_constructor(params)
    otu_picker(input_seqs_filepath,\
     result_path=result_path,log_path=log_path)
    
    
