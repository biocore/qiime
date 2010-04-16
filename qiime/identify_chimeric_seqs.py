#!/usr/bin/env python
# File created on 05 Oct 2009.

from __future__ import division

from cogent.util.misc import remove_files
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.util import FunctionWithParams
from qiime.assign_taxonomy import BlastTaxonAssigner

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

class ChimeraChecker(FunctionWithParams):
    """A ChimeraChecker takes a sequence collection and returns list of chimeric seqs.

    This is an abstract class: subclasses should implement the __call__
    method.

    Note: sequence ids should be preserved during this process, i.e. the
    description lines should be saved/restored if the alignment app is
    destructive to them.
    """
    Name = 'ChimeraChecker'

    def __init__(self, params):
        """Return new ChimeraChecker object with specified params.
        
        """
        self.Params = params
    
    def __call__ (self, seq_path, result_path=None, log_path=None):
        """Returns alignment from sequences.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
         dump the result to the desired path as fasta, otherwise will
         return list of sequence identifiers.
        log_path: path to log, which should include dump of params.
        
        """
        # return list of the chimeric sequences
        return self.getResult(seq_path)
    
    def getResult(self,seq_path):
        """ Must be overwritten in subclasses """
        raise NotImplementedError, "ChimeraChecker is an abstract class"
        
    def cleanUp(self):
        """ Overwrite if any clean-up is necessasry 
        
            Included here so all ChimeraCheckers can call cleanUp() without
             having to know if it is defined for a certain sub-class.
        """
        pass
        
class MetaChimeraChecker(ChimeraChecker):
    """ """
    Name = 'MetaChimeraChecker'
    
    def __init__(self,params):
        """ 
        """ 
        
        try:
            chimera_checker_constructors = \
             params['chimera_checker_constructors']
        except KeyError:
            raise ValueError,\
              "Must provide 'chimera_checker_constructors' to define %s." \
               % self.Name
        
        chimera_checkers = []
        for constructor, constructor_params in chimera_checker_constructors:
            chimera_checkers.append(constructor(constructor_params))
        self.chimera_checkers = chimera_checkers
        
    def getResult(self,seqs_fp):
        raise NotImplementedError, "%s is not yet written" % self.Name
        
class BlastFragmentsChimeraChecker(ChimeraChecker):
    """  """
    Name = 'BlastFragmentsChimeraChecker'
    
    def __init__(self, params):
        """Return new BlastFragmentsChimeraChecker object with specified params.
        
        """
        _params = {'max_e_value':1e-30,\
                   'min_pct_id':0.90,\
                   'num_fragments':3,\
                   'taxonomy_depth':4}
        _params.update(params)
            
        try:
            id_to_taxonomy_fp = params['id_to_taxonomy_fp']
        except KeyError:
            raise ValueError,\
             "id_to_taxonomy_filepath must be provided to %s" % self.Name
        
        # Create the blast database if it hasn't been provided
        if 'blast_db' not in params or params['blast_db'] == None:
            try:
                reference_seqs_fp = params['reference_seqs_fp']
            except KeyError:
                raise ValueError, \
                 "refseqs_fp or blast_db must be provided to  %s" % self.Name
            blast_db, self._db_files_to_remove = \
             build_blast_db_from_fasta_path(reference_seqs_fp)
        else:
            blast_db = params['blast_db']
            self._db_files_to_remove = []
            
        self._taxon_assigner = BlastTaxonAssigner(\
         {'blast_db':blast_db,\
          'id_to_taxonomy_filepath':id_to_taxonomy_fp,
          'Max E value':_params['max_e_value'],
          'Min percent identity':_params['min_pct_id']
         })
        
        ChimeraChecker.__init__(self, _params)
        
    def cleanUp(self):
        """ Remove temporary blast database files, if applicable
        """
        remove_files(self._db_files_to_remove,error_on_missing=False)
        

    def _fragment_seq(self,seq):
        """ Returns list of n roughly equal-sized fragments of seq
        
            Modified from Python Cookbook entry comment at:
             http://code.activestate.com/recipes/425397/
             
            ** This method needs to be parameterized, as there are other ways
            we might want to fragment sequences. For example, we might want to
            specify breakpoints. Currently, this will split seq into 
            self.num_fragments roughly equal-sized fragments.
        """
        num_fragments = self.Params['num_fragments']
        results = []
        start = 0
        for i in range(num_fragments): 
            # My notes:
            # len(seq[i::n]) gives the number of even multiples of
            # num_fragments exist between i (inclusive) and the end of the seq.
            stop = start + len(seq[i::num_fragments])
            results.append(seq[start:stop])
            start = stop
        return results
        
    def getResult(self,seq_path):
        """ """
        # Iterate over seq_id, seq pairs from seq_path
        for seq_id, seq in MinimalFastaParser(open(seq_path)):
            # Generate a list of fragments 
            fragments = self._fragment_seq(seq)
            # Assign the taxonomy for each of the fragments
            taxon_assignments = \
             [self._get_taxonomy(fragment) for fragment in fragments]
            # Test whether the taxon_assignments suggest that the 
            # current seq is a chimera
            if self._is_chimeric(taxon_assignments):
                # If so, yield the seq_id and the list of taxon_assignments
                yield seq_id, taxon_assignments
        
    def _is_chimeric(self,taxon_assignments):
        """From list of taxon assignments, determine if sequence is chimeric 
        
            ** This method needs to be parameterized, as there are other ways
            we might want to test for chimeras based on taxon_assignments.
            
            This method checks that taxon assignments are identical to a 
            depth of depth (the input parameter). Assignments to
            'No blast hit' are ignored, as these cannot provide evidence
            for a sequence being chimeric. This function considers a sequnce
            non-chimeric until it finds evidence that it is chimeric.
            Currently only the first blast hit is considered -- this is 
            something we will want to change. 
            
            When tested against a balanced test set containing 84 chimeric 
             sequences and 84 non-chimeric sequences 
             (precision, recall, f-measure)
             respectively are: 
              @ depth=3 (0.67,0.36,0.47)
              @ depth=4 (0.67,0.60,0.63)
              @ depth=5 (0.56,0.64,0.60)
            
        """
        depth = self.Params['taxonomy_depth']
        # split the assignments, ignoring fragments which had no blast
        # result (as these cannot provide evidence for a sequence being
        # chimeric)
        split_assignments = [ta.split(';') \
         for ta in taxon_assignments if ta != 'No blast hit']
        for i in range(depth):
            try:
                assignments_at_depth_i =\
                 [a[i].strip() for a in split_assignments]
            except IndexError:
                # We return False if there is no information at current
                # depth because there is not enough information to call
                # the sequence chimeric
                return False
            if len(set(assignments_at_depth_i)) > 1:
                return True
                
        return False
        
    def _get_taxonomy(self,fragment):
        """ Return the taxonomy of fragment 
        
            This function is awful right now, because the blast app 
             controller takes a file. This therefore results in each 
             fragment being written to file and subsequently cleaned up.
             This can be _very_ slow. I recently updated BlastTaxonAssigner
             to handle this better, but it will involve some re-writes here
             to get it all going. 
        """
        
        # Pass the temporary file to the taxon assigner, and get the 
        # taxonomy and quality score back 
        r = self._taxon_assigner(seqs=[('fragment',fragment)])['fragment']
        taxonomy, quality_score = r[0], r[1]
        
        # Return the taxonomy
        return taxonomy

def blast_fragments_identify_chimeras(seqs_fp,id_to_taxonomy_fp,\
    reference_seqs_fp=None,blast_db=None,min_pct_id=0.90,\
    max_e_value=1e-30,num_fragments=3,output_fp=None,taxonomy_depth=4):
    """ """
    params = {'id_to_taxonomy_fp': id_to_taxonomy_fp,\
              'reference_seqs_fp':reference_seqs_fp,\
              'blast_db':blast_db,\
              'min_pct_id':min_pct_id,\
              'max_e_value':max_e_value,\
              'num_fragments':num_fragments,\
              'taxonomy_depth':taxonomy_depth}
    bcc = BlastFragmentsChimeraChecker(params)
    
    if output_fp:
        of = open(output_fp,'w')
        for seq_id, tax_assignments in bcc(seqs_fp):
            of.write('\t'.join([seq_id] + tax_assignments))
            of.write('\n')
        of.close()
        result = None
    else:
        result = list(bcc(seqs_fp))
    
    bcc.cleanUp()
    return result
