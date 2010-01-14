#!/usr/bin/env python

__author__ = "Rob Knight, Greg Caporaso"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger"] 
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"


"""Contains code for assigning taxonomy, using several techniques.

This module has the responsibility for taking a set of sequences and
providing a taxon assignment for each sequence.
"""

#import csv
import logging
import re
from os import system, remove, path, mkdir
from os.path import split, splitext
from glob import glob
from itertools import count
from string import strip
from shutil import copy as copy_file
from optparse import OptionParser
from tempfile import NamedTemporaryFile
from cogent import LoadSeqs, DNA
from cogent.app.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from cogent.app.blast import blast_seqs, Blastall, BlastResult
from cogent.app.rdp_classifier import assign_taxonomy, \
    train_rdp_classifier_and_assign_taxonomy
from cogent.parse.fasta import MinimalFastaParser
from qiime.util import FunctionWithParams


class TaxonAssigner(FunctionWithParams):
    """A TaxonAssigner assigns a taxon to each of a set of sequences.

    This is an abstract class: subclasses should implement the __call__
    method.
    """
    Name = 'TaxonAssigner'

    def __init__(self, params):
        """Return new TaxonAssigner object with specified params.
        
        Note: expect params to contain both generic and per-method (e.g. for
        RDP classifier w/ Hugenholtz taxonomy) params, so leaving it as a dict 
        rather than setting attributes. Some standard entries in params are:

        Taxonomy: taxonomy used (e.g. RDP, Hugenholtz)
        Similarity: similarity threshold for assignment, e.g. 0.97
        Bootstrap: bootstrap support for assignment, e.g. 0.80
        Application: 3rd-party application used, if any, e.g. RDP classifier
        """
        self.Params = params

    def __call__ (self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {seq_id:(taxonomy, confidence)} for each seq.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError, "TaxonAssigner is an abstract class"

    @staticmethod
    def _parse_id_to_taxonomy_file(f):
        """ parse the id_to_taxonomy file into a dict mapping id -> taxonomy
        """
        result = {}
        for line in f:
            line = line.strip()
            if line:
                identifier, taxonomy = map(strip, line.split('\t'))
                result[identifier] = taxonomy
        return result 


class BlastTaxonAssigner(TaxonAssigner):
    """ Assign taxon best on best blast hit above a threshold
    """
    Name = 'BlastTaxonAssigner'
    SeqsPerBlastRun = 1000
    def __init__(self, params):
        """ Initialize the object
        """
        _params = {
            'Min percent identity': 0.90,
            'Max E value': 1e-30,
            'Application': 'blastn/megablast'
            }
        _params.update(params)
        TaxonAssigner.__init__(self, _params)
    
    def __call__(self, seq_path=None, seqs=None, result_path=None, log_path=None):
        """Returns dict mapping {seq_id:(taxonomy, confidence)} for each seq.
        """
        assert seq_path or seqs, \
         "Must provide either seqs or seq_path when calling a BlastTaxonAssigner."
         
        # initialize the logger
        logger = self._get_logger(log_path)
        logger.info(str(self))
        
        # assign the blast database, either as a pre-exisiting database
        # specified as self.Params['blast_db'] or by creating a 
        # temporary database from the sequence file specified
        # as self.Params['reference_seqs_filepath']
        try:
            blast_db = self.Params['blast_db']
        except KeyError:
            # build a temporary blast_db
            reference_seqs_path = self.Params['reference_seqs_filepath']
            refseqs_dir, refseqs_name = split(reference_seqs_path)
            blast_db, db_files_to_remove = \
             build_blast_db_from_fasta_path(reference_seqs_path)
        
        # build the mapping of sequence identifier 
        # (wrt to the blast db seqs) to taxonomy
        id_to_taxonomy_map = self._parse_id_to_taxonomy_file(\
         open(self.Params['id_to_taxonomy_filepath'])) 
        
        ## Iterate over the input self.SeqsPerBlastRun seqs at a time. 
        # There are two competing issues here when dealing with very large
        # inputs. If all sequences are read in at once, the containing object
        # can be very large, causing the system to page. On the other hand,
        # in such cases it would be very slow to treat each sequence 
        # individually, since blast requires a filepath. Each call would
        # therefore involve writing a single sequence to file, opening/closing
        # and removing the file. To balance this, sequences are read in and
        # blasted in chunks of self.SeqsPerBlastRun (defualt: 1000) at a time.
        # This appears to solve the problem with the largest sets I've worked
        # with so far. 
        
        if seq_path:
            # Get a seq iterator
            seqs = MinimalFastaParser(open(seq_path))
        # Build object to keep track of the current set of sequence to be
        # blasted, and the results (i.e., seq_id -> (taxonomy,quaility score) 
        # mapping)
        current_seqs = []
        result = {}
        
        # Iterate over the (seq_id, seq) pairs
        for seq_id, seq in seqs:
            # append the current seq_id,seq to list of seqs to be blasted
            current_seqs.append((seq_id,seq))
            
            # When there are 1000 in the list, blast them
            if len(current_seqs) == self.SeqsPerBlastRun:
                # update the result object
                result.update(self._seqs_to_taxonomy(\
                 current_seqs,blast_db,id_to_taxonomy_map))
                # reset the list of seqs to be blasted
                current_seqs = []
        # Assign taxonomy to the remaining sequences
        result.update(self._seqs_to_taxonomy(\
         current_seqs,blast_db,id_to_taxonomy_map))
        ## End iteration over the input self.SeqsPerBlastRun seqs at a time. 
        
        # Write log data if we have a path (while the logger can handle
        # being called if we are not logging, some of these steps are slow).
        if log_path is not None:
            num_inspected = len(result)
            logger.info('Number of sequences inspected: %s' % num_inspected)
            num_null_hits = [r[1] for r in result.values()].count(None)
            logger.info('Number with no blast hits: %s' % num_null_hits)

        if result_path:
            # if the user provided a result_path, write the 
            # results to file
            of = open(result_path,'w')
            for seq_id, (lineage, confidence) in result.items():
                of.write('%s\t%s\t%s\n' % (seq_id, lineage, confidence))
            of.close()
            result = None
            logger.info('Result path: %s' % result_path)
        else:
            # Returning the data as a dict, so no modification to result
            # is necessary.
            pass
                 
            # if no result_path was provided, return the data as a dict
            logger.info('Result path: None, returned as dict.')

        # clean-up temp blastdb files, if a temp blastdb was created
        if 'reference_seqs_filepath' in self.Params:
            map(remove,db_files_to_remove)

        # return the result
        return result
        
    def _seqs_to_taxonomy(self,seqs,blast_db,id_to_taxonomy_map):
        """ Assign taxonomy to (seq_id,seq) pairs
        """
        # Handle the case of no seqs passed in
        if not seqs: 
            return {}
        # blast the seqs
        blast_hits = self._get_blast_hits(blast_db,seqs)

        # select the best blast hit for each query sequence
        best_blast_hit_ids = self._get_first_blast_hit_per_seq(blast_hits)
 
        # map the identifier of the best blast hit to (taxonomy, e-value)
        return self._map_ids_to_taxonomy(\
             best_blast_hit_ids,id_to_taxonomy_map)

    def _get_logger(self, log_path=None):
        if log_path is not None:
            handler = logging.FileHandler(log_path, mode='w')
        else:
            class NullHandler(logging.Handler):
                def emit(self, record): pass
            handler = NullHandler()
        logger = logging.getLogger("BlastTaxonAssigner logger")
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def _map_ids_to_taxonomy(self, hits, id_to_taxonomy_map):
        """ map {query_id:(best_blast_seq_id,e-val)} to {query_id:(tax,None)}
        """
        for query_id, hit in hits.items():
            query_id=query_id.split()[0]
            try:
                hit_id, e_value = hit 
                hits[query_id] = \
                  (id_to_taxonomy_map.get(hit_id, None),e_value)
            except TypeError:
                hits[query_id] = ('No blast hit', None)

        return hits
        
    def _get_blast_hits(self,blast_db,seqs):
        """ blast each seq in seqs against blast_db and retain good hits
        """
        max_evalue = self.Params['Max E value']
        min_percent_identity = self.Params['Min percent identity']
        seq_ids = [s[0] for s in seqs]
        result = {}
        
        blast_result = blast_seqs(\
         seqs,Blastall,blast_db=blast_db,\
         params={'-p':'blastn','-n':'T'},\
         add_seq_names=False)
         
        if blast_result['StdOut']:
            lines = [x for x in blast_result['StdOut']]
            blast_result = BlastResult(lines)
        else:
            return {}.fromkeys(seq_ids,[])
            
        for seq_id in seq_ids:
            blast_result_id = seq_id.split()[0]
            try:
                result[seq_id] = [(e['SUBJECT ID'],float(e['E-VALUE'])) \
                 for e in blast_result[blast_result_id][0]
                 if (float(e['E-VALUE']) <= max_evalue and \
                  float(e['% IDENTITY']) >= min_percent_identity)]
            except KeyError:
                result[seq_id] = []

        return result
        
    def _get_first_blast_hit_per_seq(self,blast_hits):
        """ discard all blast hits except the best for each query sequence
        """
        result = {}
        for k,v in blast_hits.items():
            k = k.split()[0]    #get rid of spaces
            try:
                result[k] = v[0]
            except IndexError:
                # If there is no good blast hit, do we want to 
                # leave the key out, or have it point to None?
                result[k] = None
        
        return result


class RdpTaxonAssigner(TaxonAssigner):
    """Assign taxon using RDP's naive Bayesian classifier

    NOT YET IMPLEMENTED:
    * Training data to be used
    """
    Name = "RdpTaxonAssigner"
    Application = "RDP classfier"
    Citation = "Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7."
    Taxonomy = "RDP"
    _tracked_properties = ['Application','Citation','Taxonomy']

    def __init__(self, params):
        """Return new RdpTaxonAssigner object with specified params.
        
        Standard entries in params are:

        Taxonomy: taxonomy used (e.g. RDP, Hugenholtz)
        """
        _params = {
            'Confidence': 0.80,
            'id_to_taxonomy_fp':None,
            'reference_sequences_fp':None
             
            }
        _params.update(params)
        TaxonAssigner.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {seq_id:(taxonomy, confidence)} for
        each seq.
        
        Parameters:
        seq_path: path to file of sequences 
        result_path: path to file of results. If specified, dumps the 
            result to the desired path instead of returning it.
        log_path: path to log, which should include dump of params.
        """
        
        min_confidence = self.Params['Confidence']
        reference_sequences_fp = self.Params['reference_sequences_fp']
        id_to_taxonomy_fp = self.Params['id_to_taxonomy_fp']
        
        seq_file = open(seq_path, 'r')

        if reference_sequences_fp and id_to_taxonomy_fp:
            # Train and assign taxonomy
            taxonomy_file, training_seqs_file = self._generate_training_files()
            results = train_rdp_classifier_and_assign_taxonomy(
                training_seqs_file, taxonomy_file, seq_file, 
                min_confidence=min_confidence,
                classification_output_fp=result_path)
        else:
            # Just assign taxonomy
            results = assign_taxonomy(
                seq_file, min_confidence=min_confidence, output_fp=result_path)

        if log_path:
            self.writeLog(log_path)

        return results

    def _generate_training_files(self):
        """Returns a tuple of file objects suitable for passing to the
        RdpTrainer application controller.
        """
        id_to_taxonomy_file = open(self.Params['id_to_taxonomy_fp'], 'r')
        reference_seqs_file = open(self.Params['reference_sequences_fp'], 'r')

        # Generate taxonomic tree and write to file 
        tree = self._build_tree(id_to_taxonomy_file)
        id_to_taxonomy_file.seek(0)

        rdp_taxonomy_file = NamedTemporaryFile(
            prefix='RdpTaxonAssigner_taxonomy_', suffix='.txt')
        rdp_taxonomy_file.write(tree.rdp_taxonomy())
        rdp_taxonomy_file.seek(0)

        # Generate a set of training seqs and write to file
        training_seqs = self._generate_training_seqs(
            reference_seqs_file, id_to_taxonomy_file)

        rdp_training_seqs_file = NamedTemporaryFile(
            prefix='RdpTaxonAssigner_training_seqs_', suffix='.fasta')
        for rdp_id, seq in training_seqs:
            rdp_training_seqs_file.write('>%s\n%s\n' % (rdp_id, seq))
        rdp_training_seqs_file.seek(0)

        return rdp_taxonomy_file, rdp_training_seqs_file

    @staticmethod
    def _generate_training_seqs(reference_seqs_file, id_to_taxonomy_file):
        """Returns an iterator of valid training sequences in
        RDP-compatible format

        Each training sequence is represented by a tuple (rdp_id,
        seq).  The rdp_id consists of two items: the original sequence
        ID with whitespace replaced by underscores, and the lineage
        with taxa separated by semicolons.
        """
        # Rdp requires unique sequence IDs without whitespace.  Can't
        # trust user IDs to not have whitespace, so we replace all
        # whitespace with an underscore.  Classification may fail if
        # the replacement method generates a name collision.

        id_to_taxonomy_map = RdpTaxonAssigner._parse_id_to_taxonomy_file(
            id_to_taxonomy_file)

        for id, seq in MinimalFastaParser(reference_seqs_file):
            taxonomy = id_to_taxonomy_map[id]
            lineage = RdpTaxonAssigner._parse_lineage(taxonomy)
            rdp_id = '%s %s' % (re.sub('\s', '_', id), ';'.join(lineage))
            yield rdp_id, seq

    @staticmethod
    def _build_tree(id_to_taxonomy_file):
        """Returns an RdpTree object representing the taxonomic tree
        derived from the id_to_taxonomy file.
        """
        tree = RdpTaxonAssigner.RdpTree()
        for line in id_to_taxonomy_file:
            id, taxonomy = map(strip, line.split('\t'))
            lineage = RdpTaxonAssigner._parse_lineage(taxonomy)
            tree.insert_lineage(lineage)
        return tree

    @staticmethod
    def _parse_lineage(lineage_str):
        """Returns a list of taxa from the semi-colon-separated
        lineage string of an id_to_taxonomy file.
        """
        lineage = lineage_str.strip().split(';')
        # The RDP Classifier can only deal with a lineage that is 6
        # levels deep.  We detect this problem now to avoid an
        # ApplicationError later on.
        if len(lineage) != 6:
            raise ValueError(
                'Each reference assignment must contain 6 items, specifying '
                'domain, phylum, class, order, family, and genus.  '
                'Detected %s items in "%s": %s.' % \
                (len(lineage), lineage_str, lineage))
        return lineage

    # The RdpTree class is defined as a nested class to prevent the
    # implementation details of creating RDP-compatible training files
    # from being exposed at the module level.
    class RdpTree(object):
        """Simple, specialized tree class used to generate a taxonomy
        file for the Rdp Classifier.
        """
        taxonomic_ranks = ['domain', 'phylum', 'class', 
                           'order', 'family', 'genus']

        def __init__(self, name='Root', parent=None):
            self.name = name
            self.parent = parent
            if parent is None:
                self.depth = 0
            else:
                self.depth = parent.depth + 1
            self.children = dict()  # name => subtree

        def insert_lineage(self, lineage):
            """Inserts a lineage into the taxonomic tree.
            
            Lineage must support the iterator interface, or provide an
            __iter__() method that returns an iterator.
            """
            lineage = lineage.__iter__()
            try:
                taxon = lineage.next()
                if taxon not in self.children:
                    self.children[taxon] = RdpTaxonAssigner.RdpTree(name=taxon, parent=self)
                self.children[taxon].insert_lineage(lineage)
            except StopIteration:
                pass
            return self

        def rdp_taxonomy(self, counter=None):
            """Returns a string, in Rdp-compatible format.
            """
            if counter is None:
                # initialize a new counter to assign IDs
                counter = count(0)

            # Assign ID to current node; used by child nodes
            self.id = counter.next()

            if self.parent is None:
                # do not print line for Root node
                retval = ''
            else:
                # Rdp taxonomy file does not count the root node
                # when considering the depth of the taxon.  We
                # therefore specify a taxon depth which is one
                # less than the tree depth.
                taxon_depth = self.depth - 1

                # In this simplest-possible implementation, we
                # also use the taxon depth to retrieve the string
                # value of the taxonomic rank.  Taxa beyond the
                # depth of our master list are given a rank of ''.
                try:
                    taxon_rank = self.taxonomic_ranks[taxon_depth]
                except IndexError:
                    taxon_rank = ''

                fields = [self.id, self.name, self.parent.id, taxon_depth,
                          taxon_rank]
                retval = '*'.join(map(str, fields)) + "\n"

            # Recursively append lines from sorted list of subtrees
            child_names = self.children.keys()
            child_names.sort()
            subtrees = [self.children[name] for name in child_names]
            for subtree in subtrees:
                retval += subtree.rdp_taxonomy(counter)
            return retval

usage_str = """usage: %prog [options] {-i INPUT_SEQUENCES_FILEPATH}

[] indicates optional input (order unimportant) 
{} indicates required input (order unimportant) 

Example usage:
Assign taxonomy of sequences in inseqs.fasta (-i) using the RDP classifier 
 (-m). Output files will be written to rdp_assigned_taxonomy (default).
 python assign_taxonomy.py -i inseqs.fasta -m rdp
 
Assign taxonomy of sequences in inseqs.fasta (-i) using the RDP classifier
 (-m) trained on-the-fly from provided refseqs and taxon assignments (-r, -t) 
 respectively. Output files will be written to custom_rdp.
 python assign_taxonomy.py -i inseqs.fasta -r refseqs.fasta -t id_to_taxonomy.txt -m rdp -o custom_rdp
 
Assign taxonomy of sequences in inseqs.fasta (-i) using BLAST
 (-m) against provided refseqs and taxon assignments (-r, -t) 
 respectively. Output files will be written to blast_assigned_taxonomy (default).
 python assign_taxonomy.py -i inseqs.fasta -r at_refseqs.fasta -t at_id_to_taxonomy.txt -m blast
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i', '--input_seqs_fp',
        help='Path to fasta file of sequences to be assigned [REQUIRED]')

    parser.add_option('-t', '--id_to_taxonomy_fp',
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a comma-separated '
         'list. For assignment with rdp, each assigned taxonomy must be '
         'exactly 6 levels deep. [default: %default; REQUIRED when method is '
         'blast]')

    parser.add_option('-r', '--reference_seqs_fp',
        help='Path to reference sequences.  For assignment with blast, these '
        'are used to generate a blast database. For assignment with rdp, they '
        'are used as training sequences for the classifier.'
        '[default: %default; REQUIRED if -b is not provided when method is blast]')

    assignment_method_choices = assignment_method_constructors.keys()
    parser.add_option('-m','--assignment_method',\
          type='choice',help='Taxon assignment method [default:%default]',\
          choices=assignment_method_choices)
          
    parser.add_option('-b', '--blast_db',
        help='Database to blast against.  Must provide either --blast_db or '
        '--reference_seqs_db for assignment with blast [default: %default]')

    parser.add_option('-c', '--confidence', type='float',
        help='Minimum confidence to record an assignment, only used for rdp '
        'method [default: %default]')

    parser.add_option('-e', '--e_value', type='float',
        help='Maximum e-value to record an assignment, only used for blast '
        'method [default: %default]')

    parser.add_option('-o','--output_dir',\
          help='Path to store result file '+\
          '[default: <ASSIGNMENT_METHOD>_assigned_taxonomy]')

    parser.set_defaults(
        verbose=False,
        assignment_method='rdp',
        confidence=0.80,
        e_value=0.001,
        )

    opts,args = parser.parse_args()

    required_options = ['input_seqs_fp']    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    if opts.assignment_method == 'blast':
        if not opts.id_to_taxonomy_fp:
            parser.error('Option --id_to_taxonomy_fp is required when ' 
                         'assigning with blast.')
        if not (opts.reference_seqs_fp or opts.blast_db):
            parser.error('Either a blast db (via -b) or a collection of '
                         'reference sequences (via -r) must be passed to '
                         'assign taxonomy using blast.')

    if opts.assignment_method == 'rdp':
        if opts.id_to_taxonomy_fp:
            if opts.reference_seqs_fp is None:
                parser.error('A filepath for reference sequences must be '
                             'specified (via -r) along with the id_to_taxonomy '
                             'file to train the Rdp Classifier.')
        elif opts.reference_seqs_fp:
                parser.error('A filepath for an id to taxonomy map must be '
                             'specified (via -t) along with the reference '
                             'sequences fp to train the Rdp Classifier.')
        else:
            pass
    
    return opts, args


assignment_method_constructors = {
    'blast': BlastTaxonAssigner,
    'rdp': RdpTaxonAssigner,
    }


if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    assignment_method = opts.assignment_method
    taxon_assigner_constructor =\
     assignment_method_constructors[assignment_method]
    input_sequences_filepath = opts.input_seqs_fp
    
    try:
        id_to_taxonomy_fp = opts.id_to_taxonomy_fp
        params = {'id_to_taxonomy_filepath':id_to_taxonomy_fp}
    except IndexError:
        params = {}
    
    # Build the output filenames
    output_dir = opts.output_dir or assignment_method + '_assigned_taxonomy'
    try:
        mkdir(output_dir)
    except OSError:
        # output_dir already exists
        pass
        
    fpath, ext = splitext(input_sequences_filepath)
    input_dir, fname = split(fpath)
    result_path = output_dir + '/' + fname + '_tax_assignments.txt'
    log_path = output_dir + '/' + fname + '_tax_assignments.log'
    
    if opts.assignment_method == 'blast':
        # one of these must have a value, otherwise we'd have 
        # an optparse error
        if opts.blast_db:
            params['blast_db'] = opts.blast_db
        else:
            params['reference_seqs_filepath'] = opts.reference_seqs_fp
        params['Max E value'] = opts.e_value

    elif opts.assignment_method == 'rdp':
        params['Confidence'] = opts.confidence
        params['id_to_taxonomy_fp'] = opts.id_to_taxonomy_fp
        params['reference_sequences_fp'] = opts.reference_seqs_fp

    else:
        # should not be able to get here as an unknown classifier would
        # have raised an optparse error
        exit(1)

    taxon_assigner = taxon_assigner_constructor(params)
    taxon_assigner(input_sequences_filepath,\
     result_path=result_path,log_path=log_path)


