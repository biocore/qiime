#!/usr/bin/env python

__author__ = "Rob Knight, Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

#import csv
import logging
import re
from os import system, remove, path, mkdir
from os.path import split, splitext
from glob import glob
from itertools import count
from string import strip
from shutil import copy as copy_file
from tempfile import NamedTemporaryFile
from cogent import LoadSeqs, DNA
from cogent.app.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from cogent.app.blast import blast_seqs, Blastall, BlastResult
from cogent.app.rdp_classifier import assign_taxonomy, \
    train_rdp_classifier_and_assign_taxonomy
from cogent.parse.fasta import MinimalFastaParser
from qiime.util import FunctionWithParams


"""Contains code for assigning taxonomy, using several techniques.

This module has the responsibility for taking a set of sequences and
providing a taxon assignment for each sequence."""

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
            for seq_id, (lineage, confidence, blast_hit_id) in result.items():
                of.write('%s\t%s\t%s\t%s\n' % 
                 (seq_id, lineage, confidence, blast_hit_id))
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
        """ map {query_id:(best_blast_seq_id,e-val)} to {query_id:(tax,e-val,best_blast_seq_id)}
        """
        for query_id, hit in hits.items():
            query_id=query_id.split()[0]
            try:
                hit_id, e_value = hit 
                hits[query_id] = \
                  (id_to_taxonomy_map.get(hit_id, None),e_value,hit_id)
            except TypeError:
                hits[query_id] = ('No blast hit', None, None)

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
