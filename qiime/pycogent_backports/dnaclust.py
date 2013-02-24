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

        # Allow for gaps on the left of the shorter string in semi-global alignment.
        'left-gaps-allowed':FlagParameter('--', Name='left-gaps-allowed')
    }

    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self, data):
        """ Set the input path (a fasta filepath).
        """
        
        if data:
            self.Parameters['input-file'].on(data)
        
        return ''

class DnaclustRef(CommandLineApplication):
    """ DNACLUST-ref Application controller
    """

    _command = 'dnaclust-ref'
    _input_handler = '_input_as_parameters'
    _parameters = {\

        # Similarity percent for clustering.  By default, 99%.
        'similarity':ValuedParameter('-', Name='r', Delimiter=' ',
            IsPath=False, Value='0.97'),
        
        # Fasta input file to be clustered.
        'input-file':ValuedParameter('-', Name='i', Delimiter=' ',
            IsPath=True),

        # Fasta input file of the reference database.
        'centers-file':ValuedParameter('-', Name='c', Delimiter=' ',
            IsPath=True),
        
        # After clustering with a reference perform de novo clustering.
        'de-novo-cluster':FlagParameter('-', Name='d'),

        # Cluster sequences such that the cluster centers are at least twice the
        # radius of the other centers.
        'no-overlap':FlagParameter('-', Name='n'),

        # Number of threads to use.
        'threads':ValuedParameter('-', Name='t', Delimiter=' ',
            IsPath=False, Value='1'),

        # Allow for gaps on the left of the shorter string in semi-global alignment
        'left-gaps-allowed':FlagParameter('-', Name='l')
    }

    _suppress_stdout = False
    _suppress_stderr = False

    def _input_as_parameters(self, data):
        """ Set the input path (a fasta filepath)
        """
        
        if data:
            self.Parameters['input-file'].on(data)
        
        return ''

    def _get_result_paths(self,data):
        """ Set the resultpaths equal to *.db.clusters and *.denovo.clusters.
        """
        result = {}

        extensions = ['db.clusters', 'denovo.clusters']

        for extension in extensions:
            result[extension] = ResultPath(\
             Path=data + '.' + extension,\
             IsWritten=True)
            print str(data + '.' + extension)
        return result

def dnaclust_from_seqs(seq_path,
        similarity=0.97,
        no_overlap=False,
        threads=1,
        left_gaps_allowed=False,
        reference_path=None,
        HALT_EXEC=False):
    """Helper function to run dnaclust/dnaclust-ref.

    seq_path = fasta sequences filepath.
    similarity = similarity threshold for clustering.
    no_overlap = choose cluster centers such that no sequence 
     falls within the radius of another.
    threads = number of threads to use.
    left_gaps_allowed = allow gaps in the left end of alignment.
    reference_path = path of the reference fasta file.
    HALT_EXEC = used for debugging app controller.
    """

    params = {'threads': threads,
        'similarity': similarity}
    
    if reference_path:
        app = DnaclustRef(params, HALT_EXEC=HALT_EXEC)
        app.Parameters['centers-file'].on(reference_path)
        app.Parameters['de-novo-cluster'].on()
    else:
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

def dnaclust_get_clusters(reference_path=None,
        results=None):
    """ Builds a dictionary of the cluster center: [clusters].

    reference_path = path of the reference fasta file.
    results = results dictionary after running the application.
    """

    clusters = {}

    cluster_file = results['StdOut']
    if reference_path:
        cluster_file = results['denovo.clusters']

    # DNACLUST does not give denovo clusterings a name,
    # so we just prepend a count to each cluster.
    otu_count = 0
    for line in cluster_file.readlines():
        clusters['dnaclust_' + str(otu_count)] = line.strip().split()
        otu_count += 1

    # If we have a reference fasta, the database clustered
    # sequences are printed to a special file.
    if reference_path:
        cluster_file = results['db.clusters']

        for line in cluster_file.readlines():
            center_and_seqs = line.strip().split()
            clusters[center_and_seqs[0]] = center_and_seqs[1:]

    return clusters