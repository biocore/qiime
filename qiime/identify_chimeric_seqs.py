#!/usr/bin/env python
# File created on 05 Oct 2009.

from __future__ import division
from optparse import OptionParser
from os.path import split, splitext
from cogent.util.misc import remove_files
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.util import FunctionWithParams
from qiime.assign_taxonomy import BlastTaxonAssigner

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, Qiime"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

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
        try:
            max_e_value = params['max_e_value']
        except KeyError:
            max_e_value = 1e-30
            
        try:
            min_pct_id = params['min_pct_id']
        except KeyError:
            min_pct_id = 0.90
            
        try:
            self.num_fragments = params['num_fragments']
        except KeyError:
            self.num_fragments = 3
            
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
            
        self._taxon_assigner = BlastTaxonAssigner(\
         {'blast_db':blast_db,\
          'id_to_taxonomy_filepath':id_to_taxonomy_fp,
          'Max E value':max_e_value,
          'Min percent identity':min_pct_id
         })
        
        ChimeraChecker.__init__(self, params)
        
    def cleanUp(self):
        """ Remove temporary blast database files, if applicable
        """
        remove_files(self._db_files_to_remove,error_on_missing=False)
        

    def fragmentSeq(self,seq):
        """ Returns list of n roughly equal-sized fragments of seq
        
            Modified from Python Cookbook entry comment at:
             http://code.activestate.com/recipes/425397/
             
            ** This method needs to be parameterized, as there are other ways
            we might want to fragment sequences. For example, we might want to
            specify breakpoints. Currently, this will split seq into 
            self.num_fragments roughly equal-sized fragments.
        """
        num_fragments = self.num_fragments
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
            fragments = self.fragmentSeq(seq)
            # Assign the taxonomy for each of the fragments
            taxon_assignments = \
             [self.getTaxonomy(fragment) for fragment in fragments]
            # Test whether the taxon_assignments suggest that the 
            # current seq is a chimera
            if self.isChimeric(taxon_assignments):
                # If so, yield the seq_id and the list of taxon_assignments
                yield seq_id, taxon_assignments

    def isChimericStrict(self,taxon_assignments):
        """From list of taxon assignments, determine if sequence is chimeric 
        
            ** This method needs to be parameterized, as there are other ways
            we might want to test for chimeras based on taxon_assignments. 
            This method requires fully identical taxonomy assignments,
            but we may want to be more flexible. For example, in cases 
            where assignments differ only in the most specific assignment 
            (species or genus, depending on the taxonomy)  we may be losing
            of interesting sequences representing a new genus.
            
        """
        return len(set(taxon_assignments)) > 1
        
    def isChimeric(self,taxon_assignments,depth=4):
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
        
    def getTaxonomy(self,fragment):
        """ Return the taxonomy of fragment 
        
            This function is awful right now, because the taxon 
             assigners take a file. Will update this to either do
             all fragments at once, or will update the taxonomy 
             assigners to optionally take a sequence.
        
        """
        # Write the fragment to a temporary file
        tmp_filename = get_tmp_filename(prefix='%s_' % self.Name,\
            suffix='.fasta')
        f = open(tmp_filename,'w')
        f.write('>fragment\n%s\n' % fragment)
        f.close()
        
        # Pass the temporary file to the taxon assigner, and get the 
        # taxonomy and quality score back 
        taxonomy, quality_score = \
         self._taxon_assigner(tmp_filename)['fragment']
        
        # Clean up the temp file
        remove_files([tmp_filename])
        
        # Return the taxonomy
        return taxonomy

def blast_fragments_identify_chimeras(seqs_fp,id_to_taxonomy_fp,\
    reference_seqs_fp=None,blast_db=None,min_pct_id=0.90,\
    max_e_value=1e-30,num_fragments=3,output_fp=None):
    """ """
    params = {'id_to_taxonomy_fp': id_to_taxonomy_fp,\
              'reference_seqs_fp':reference_seqs_fp,\
              'blast_db':blast_db,\
              'min_pct_id':min_pct_id,\
              'max_e_value':max_e_value,\
              'num_fragments':num_fragments}
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


usage_str = """usage: %prog [options] {-i INPUT_FASTA_FILE}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:

Test whether each seq in inseqs.fasta (-i) is chimeric by splitting into 
 4 (-n) roughly equal-sized fragments, and assigning taxonomy to each fragment
 using the blast taxon assigner (-t and -r, see qiime.assign_taxonomy.py -h).
 Chimeras are sequences where different fragments are assigned to different
 taxonomies. Seq ids for putative chimeras will be written to 
 inseqs_chimeric.txt (default, derived from -i) along with the taxonomies 
 assigned to each fragment.

 python /Qiime/qiime/identify_chimeric_seqs.py -i inseqs.fasta -t id_to_taxonomy.txt -r refseqs.fasta -n 4
"""

chimera_detection_method_choices = ['blast_fragments']

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i', '--input_seqs_fp',
        help='Path to fasta file of sequences to be assigned [REQUIRED]')

    parser.add_option('-t', '--id_to_taxonomy_fp',
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a comma-separated '
         'list. [default: %default; REQUIRED when method is blast_fragments]')

    parser.add_option('-r', '--reference_seqs_fp',
        help='Path to reference sequences.  For assignment with blast, these '
        'are used to generate a blast database. For assignment with rdp, they '
        'are used as training sequences for the classifier '
        '[default: %default; REQUIRED when method is blast_fragments]')
        
    parser.add_option('-m','--chimera_detection_method',\
          type='choice',help='Chimera detection method [default:%default]',\
          choices=chimera_detection_method_choices)
          
    parser.add_option('-n','--num_fragments',\
          type='int',help='Number of fragments to split sequences into' +\
          ' (i.e., number of expected breakpoints + 1) [default: %default]')
          
    parser.add_option('-o', '--output_fp',
        help='Path to store output [derived from input_seqs_fp]')

    parser.add_option('-v','--verbose',action='store_true',\
        dest='verbose',help='Print information during execution -- '+\
        'useful for debugging [default: %default]')

    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False,chimera_detection_method='blast_fragments',\
        num_fragments=3)

    opts,args = parser.parse_args()
    required_options = ['input_seqs_fp']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    if opts.num_fragments < 2:
        parser.error('Invalid number of fragments (-n %d) Must be >= 2.' \
         % opts.num_fragments)

    return opts,args

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    input_seqs_fp = opts.input_seqs_fp
    id_to_taxonomy_fp = opts.id_to_taxonomy_fp
    reference_seqs_fp = opts.reference_seqs_fp
    chimera_detection_method = opts.chimera_detection_method
    num_fragments = opts.num_fragments
    output_fp = opts.output_fp
    
    if not output_fp:
        input_basename = splitext(split(input_seqs_fp)[1])[0]
        output_fp = '%s_chimeric.txt' % input_basename
    
    if chimera_detection_method == 'blast_fragments':
        blast_fragments_identify_chimeras(input_seqs_fp,id_to_taxonomy_fp,\
            reference_seqs_fp,num_fragments=opts.num_fragments,output_fp=output_fp)
        
    