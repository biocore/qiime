#!/usr/bin/env python

__author__ = "Rob Knight, Greg Caporaso"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Rob Knight","Greg Caporaso", "Kyle Bittinger"] 
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"


"""Contains code for assigning taxonomy, using several techniques.

This module has the responsibility for taking a set of sequences and
providing a taxon assignment for each sequence.
"""

import logging
import re
from os import system, remove
from glob import glob
from itertools import imap
from shutil import copy as copy_file
from optparse import OptionParser
from cogent import LoadSeqs, DNA
from cogent.app.util import get_tmp_filename
from cogent.app.blast import blast_seqs, Blastall, BlastResult
from cogent.app.rdp_classifier import assign_taxonomy
from cogent.parse.fasta import MinimalFastaParser
from pipe454.util import FunctionWithParams

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


class BlastTaxonAssigner(TaxonAssigner):
    """ Assign taxon best on best blast hit above a threshold
    """
    Name = 'BlastTaxonAssigner'
    
    def __init__(self,params):
        """ Initialize the object
        
        """
        _params = {'Min percent identity':0.90,\
         'Max E value':1e-30,\
         'Application':'blastn/megablast'}
        _params.update(params)
        TaxonAssigner.__init__(self, _params)
    
    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns dict mapping {seq_id:(taxonomy, confidence)} for each seq.
        """
        # Load the sequence collection containing the files to blast
        seq_coll = LoadSeqs(seq_path,aligned=False,moltype=DNA).degap()
        
        # assign the blast database, either as a pre-exisiting database
        # specified as self.Params['blast_db'] or by creating a 
        # temporary database from the sequence file specified
        # as self.Params['reference_seqs_filepath']
        try:
            blast_db = self.Params['blast_db']
        except KeyError:
            # build the blast_db
            reference_seqs_path = self.Params['reference_seqs_filepath']
            # get a temp file name for the blast database
            blast_db = get_tmp_filename(tmp_dir='',\
             prefix='BlastTaxonAssigner_temp_db_',suffix='')
            # copy the reference seqs file so formatdb can safely
            # create files based on the filename (without the danger of
            # overwriting existing files)
            copy_file(reference_seqs_path,blast_db)    
            # create the blast db
            if system('formatdb -i ' + blast_db + ' -o T -p F') != 0:
                # WHAT TYPE OF ERROR SHOULD BE RAISED IF THE BLAST_DB
                # BUILD FAILS?
                raise RuntimeError,\
                 "Creation of temporary Blast database failed."
            
            # create a list of the files to clean-up
            db_files_to_remove = ['formatdb.log'] + glob(blast_db + '*')
        
        # build the mapping of sequence identifier (wrt to the blast db seqs)
        # to taxonomy
        id_to_taxonomy_map = self._parse_id_to_taxonomy_file(\
         open(self.Params['id_to_taxonomy_filepath'])) 
        # blast the sequence collection against the database
        blast_hits = self._get_blast_hits(blast_db,seq_coll)
        # select the best blast hit for each query sequence
        best_blast_hit_ids = self._get_first_blast_hit_per_seq(blast_hits)
         
        # map the identifier of the best blast hit to (taxonomy, e-value)
        seq_id_to_taxonomy = self._map_ids_to_taxonomy(\
             best_blast_hit_ids,id_to_taxonomy_map)
             
        if result_path:
            # if the user provided a result_path, write the 
            # results to file
            of = open(result_path,'w')
            for seq_id, data in seq_id_to_taxonomy.items():
                if data:
                    of.write('%s\t%s\t%s\n' % (seq_id,data[0],str(data[1])))
                else:
                    # write to failures report file -- need specs on where
                    # this file comes from; for now, just ignore it
                    pass
            of.close()
            result = None
            log_str = 'Result path: %s' % result_path
        else:
            # if no result_path was provided, return the data as a dict
            result = seq_id_to_taxonomy
            log_str = 'Result path: None, returned as dict.'

        # clean-up temp blastdb files, if a temp blastdb was created
        if 'reference_seqs_filepath' in self.Params:
            map(remove,db_files_to_remove)
            
        if log_path:
            # if the user provided a log file path, log the run
            log_file = open(log_path,'w')
            log_file.write(str(self))
            log_file.write('\n')
            log_file.write('%s\n' % log_str)
        
        # return the result
        return result
        
    def _map_ids_to_taxonomy(self,hits,id_to_taxonomy_map): 
        """ map {query_id:(best_blast_seq_id,e-val)} to {query_id:(tax,None)}
        """
        for query_id, hit in hits.items():
            query_id=query_id.split()[0]
            try:
                hit_id, e_value = hit 
                hits[query_id] = \
                  (id_to_taxonomy_map.get(hit_id, None),None)
            except TypeError:
                hits[query_id] = None

        return hits
        
    def _parse_id_to_taxonomy_file(self,f):
        """ parse the id_to_taxonomy file into a dict mapping id -> taxonomy
        """
        result = {}
        for line in f:
            line = line.strip()
            if line:
                identifier, taxonomy = line.split('\t')
                result[identifier] = taxonomy
        return result 

    def _get_blast_hits(self,blast_db,seq_coll):
        """ blast each seq in seq_coll against blast_db and retain good hits
        """
        max_evalue = self.Params['Max E value']
        min_percent_identity = self.Params['Min percent identity']
        result = {}
        
        blast_result = blast_seqs(\
         seq_coll.toFasta(),Blastall,blast_db=blast_db,\
         params={'-p':'blastn','-n':'T'},\
         input_handler='_input_as_multiline_string',
         add_seq_names=False)
         
        if blast_result['StdOut']:
            lines = [x for x in blast_result['StdOut']]
            blast_result = BlastResult(lines)
        else:
            return {}.fromkeys(seq_coll.Names,[])
            
        for seq in seq_coll.iterSeqs():
            seq_id = seq.Name
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
        
        if reference_sequences_fp and id_to_taxonomy_fp:
            raise NotImplementedError,\
             "Training on-the-fly is not yet implemented."
        
        results = assign_taxonomy(\
             open(seq_path),min_confidence=min_confidence,\
             output_fp=result_path)

        if log_path:
            self.writeLog(log_path)

        return results

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage =\
     'usage: %prog [options] input_sequences_filepath [id_to_taxonomy_filepath]'
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-m','--assignment_method',action='store',\
          type='string',dest='assignment_method',help='Method for assigning'+\
          ' taxonomy [default: %default]')
          
    parser.add_option('-b','--blast_db',action='store',\
          type='string',dest='blast_db',help='Database to blast'+\
          ' against [default: %default]')
          
    parser.add_option('-c','--confidence',action='store',\
          type='float',dest='confidence',help='Minimum confidence to'+\
          ' record an assignment [default: %default]')
          
    parser.add_option('-r','--reference_seqs_fp',action='store',\
          type='string',dest='reference_seqs_fp',help='Path to reference '+\
          ' sequences to blast against [default: %default]')
          
    parser.add_option('-o','--result_fp',action='store',\
          type='string',dest='result_fp',help='Path to store '+\
          'result file [default: <input_sequences_filepath>.txt]')
          
    parser.add_option('-l','--log_fp',action='store',\
          type='string',dest='log_fp',help='Path to store '+\
          'log file [default: No log file created.]')

    parser.set_defaults(verbose=False,assignment_method='blast',\
     confidence=0.80)

    opts,args = parser.parse_args()
    if opts.assignment_method == 'blast' and len(args) != 2:
       parser.error(\
        'Exactly two arguments are required when assigning with blast.')
       
    if opts.assignment_method == 'rdp' and len(args) != 1:
       parser.error(\
        'Exactly one argument is required when assigning with RDP.')
       
    if opts.assignment_method not in assignment_method_constructors:
        parser.error(\
         'Invalid assignment method: %s.\nValid choices are: %s'\
         % (opts.assignment_method,\
            ' '.join(assignment_method_constructors.keys())))
            
    if opts.assignment_method == 'blast' and \
     not (opts.reference_seqs_fp or opts.blast_db):
        parser.error('Either a blast db (via -b) or a ' +\
        'collection of reference sequences (via -r) must be ' +\
        'passed to assign taxonomy using blast.')

    return opts,args

assignment_method_constructors = dict([
        ('blast', BlastTaxonAssigner),
        ('rdp', RdpTaxonAssigner),
        ])

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
 
    taxon_assigner_constructor =\
     assignment_method_constructors[opts.assignment_method]
     
    input_seqs_filepath = args[0]
    try:
        id_to_taxonomy_fp = args[1] 
        params = {'id_to_taxonomy_filepath':id_to_taxonomy_fp}
    except IndexError:
        params = {}
   
    result_path = opts.result_fp or \
     input_seqs_filepath.replace('.fasta','_tax_assignments.txt')
     
    log_path = opts.log_fp
    
    if opts.assignment_method == 'blast':
        # one of these must have a value, otherwise we'd have 
        # an optparse error
        if opts.blast_db:
            params['blast_db'] = opts.blast_db
        else:
            params['reference_seqs_filepath'] = opts.reference_seqs_fp
    elif opts.assignment_method == 'rdp':
        params['Confidence'] = opts.confidence
    else:
        # should not be able to get here as an unknown classifier would
        # have raised an optparse error
        exit(1)
            
    
    taxon_assigner = taxon_assigner_constructor(params)
    taxon_assigner(input_seqs_filepath,\
     result_path=result_path,log_path=log_path)


