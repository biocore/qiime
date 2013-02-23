#!/usr/bin/env python
"""Application controller for DNACLUST

Includes application controllers for DNACLUST and
convenience wrappers for different functions of DNACLUST including
sorting fasta files, finding clusters, and searching and aligning 
against a database. 

Written by Christopher Hill
"""

__author__ = "Christopher Hill"
__copyright__ = "Copyright 2007-2013, The PyCogent Project"
__credits__ = ["Christopher Hill"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Christopher Hill"
__email__ = "cmhill@umd.edu"
__status__ = "Production"

from os.path import split, splitext, basename, isdir, abspath, isfile, join

from cogent.app.util import CommandLineApplication, ResultPath, \
 ApplicationError
from cogent.app.parameters import ValuedParameter, FlagParameter

"""  -h [ --help ]                              produce help message
  -s [ --similarity ] arg (=0.99)            set similarity between cluster 
                                             center and cluster sequences
  -i [ --input-file ] arg                    input file
  -p [ --predetermined-cluster-centers ] arg file containing predetermined 
                                             cluster centers
  -r [ --recruit-only ]                      when used with 'predetermined-clus
                                             ter-centers' option, only clusters
                                             the input sequences that are 
                                             similar to the predetermined 
                                             centers
  -d [ --header ]                            output header line indicating run 
                                             options
  -l [ --left-gaps-allowed ]                 allow for gaps on the left of 
                                             shorter string in semi-global 
                                             alignment
  -k [ --k-mer-length ] arg                  length of k-mer for filtering
  --approximate-filter                       use faster approximate k-mer 
                                             filter
  --no-k-mer-filter                          do not use k-mer filter
  --no-overlap                               cluster some of sequences such 
                                             that the cluster centers are at 
                                             distance at least two times the 
                                             radius of the clusters
  -t [ --threads ] arg                       number of threads
  -u [ --use-full-query-header ]             use the full query header instead 
                                             of the first word
"""

class Dnaclust(CommandLineApplication):
    """ DNACLUST Application controller

    """

    _command = 'dnaclust'
    _input_handler = '_input_as_parameters'
    _parameters = {\

        # Similarity percent for clustering.  By default, 99%.
        'similarity':ValuedParameter('--', Name='similarity', Delimiter=' ',
            IsPath=False, Value='0.97'),
        
        # Fasta input file to be clustered.
        'input-file':ValuedParameter('--', Name='input-file', Delimiter=' ',
            IsPath=True),

        # Don't use the k-mer filter heuristic.
        'no-k-mer-filter':FlagParameter('--', Name='no-k-mer-filter'),

        # Cluster sequences such that the cluster centers are at least twice the
        # radius of the other centers.
        'no-overlap':FlagParameter('--', Name='no-overlap'),

        # Number of threads to use.
        'threads':ValuedParameter('--', Name='threads', Delimiter=' ',
            IsPath=False, Value='1'),

        # Allow for gaps on the left of the shorter string in semi-global alignment
        'left-gaps-allowed':FlagParameter('--', Name='left-gaps-allowed')
    }
    #_synonyms = {'similarity':'s'}

    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self, data):
        """ Set the input path (a fasta filepath)
        """
        #print self.Parameters['--left-gaps-allowed'].isOn()
        #print self.Parameters['--no-overlap'].isOn()
        #print 'param'
        if data:
            self.Parameters['input-file'].on(data)
        
        print 'test'
        return ''

    """def _get_result_paths(self, data):
        result = {}
        name_counter = 0
        seq_counter = 0
        if not isinstance(data,list):
            #means data is file
            data = open(data).readlines()
        for item in data:

        return result
    """

def dnaclust_from_seqs(seq_path,
        no_overlap=False,
        threads=1,
        left_gaps_allowed=False,
        HALT_EXEC=False,
        params=None):

    params = {}
    
    app = Dnaclust(params, HALT_EXEC=HALT_EXEC)
    
    app.Parameters['no-k-mer-filter'].on()

    if left_gaps_allowed:
        app.Parameters['left-gaps-allowed'].on()
        
    if no_overlap:
        app.Parameters['no-overlap'].on()

    try:
        app_result = app(seq_path)
    except ApplicationError:
        raise ValueError, ('No data following filter steps, please check '+\
         'parameter settings for usearch_qf.')
    
    return app_result


