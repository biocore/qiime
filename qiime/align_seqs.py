#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Jeremy Widmann", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"


"""Contains code for aligning sequences, using several techniques.

This module has the responsibility for taking a set of sequences and
returning an alignment. Mostly, it will be thin wrappers for code 
already in cogent.app.*, to which wrappers for e.g. PyNAST need to be
added..
"""
from os import remove, makedirs
from os.path import exists, splitext, split
from commands import getoutput
from optparse import OptionParser
from cogent import LoadSeqs, DNA
from cogent.core.alignment import DenseAlignment, SequenceCollection, Alignment
from cogent.core.sequence import DnaSequence as Dna
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.util import get_tmp_filename, ApplicationNotFoundError
from qiime.util import FunctionWithParams, qiime_config
from cogent.app.infernal import cmalign_from_alignment
from cogent.parse.rfam import MinimalRfamParser, ChangedSequence
#app controllers that implement align_unaligned_seqs
import cogent.app.muscle
import cogent.app.clustalw
import cogent.app.mafft

# Load PyNAST if it's available. If it's not, skip it if not but set up
# to raise errors if the user tries to use it.
try:
    from pynast.util import pynast_seqs, pair_hmm_align_unaligned_seqs,\
        muscle_align_unaligned_seqs, mafft_align_unaligned_seqs,\
        clustal_align_unaligned_seqs, blast_align_unaligned_seqs
    from pynast.logger import NastLogger

except ImportError:
    def raise_pynast_not_found_error(*args, **kwargs):
        raise ApplicationNotFoundError,\
         "PyNAST cannot be found.\nIs PyNAST installed? Is it in your $PYTHONPATH?"+\
         "\nYou can obtain PyNAST from http://pynast.sourceforge.net/." 
    # set functions which cannot be imported to raise_pynast_not_found_error
    pynast_seqs = NastLogger = pair_hmm_align_unaligned_seqs = \
    muscle_align_unaligned_seqs = mafft_align_unaligned_seqs =\
    clustal_align_unaligned_seqs = blast_align_unaligned_seqs = \
    raise_pynast_not_found_error


class Aligner(FunctionWithParams):
    """An Aligner takes an unaligned set of sequences and returns an alignment.

    This is an abstract class: subclasses should implement the __call__
    method.

    Note: sequence ids should be preserved during this process, i.e. the
    description lines should be saved/restored if the alignment app is
    destructive to them.
    """
    Name = 'Aligner'

    def __init__(self, params):
        """Return new Aligner object with specified params.
        
        Note: expect params to contain both generic and per-method (e.g. for
        infernal vs. PyNAST vs. whatever) params, so leaving it as a dict 
        rather than setting attributes. Some standard entries in params are:

        Application: 3rd-party application used, if any, e.g. infernal
        [can't actually think of any other params that apply to all of
         e.g. PyNAST, infernal, and muscle]
        """
        self.Params = params

    def __call__ (self, seq_path, result_path=None, log_path=None):
        """Returns alignment from sequences.
        
        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path as fasta, otherwise should
        return cogent.core.alignment.DenseAlignment object.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError, "Aligner is an abstract class"

class CogentAligner(Aligner):
    """Generic aligner using Cogent multiple alignment methods."""

    Name = 'CogentAligner'

    def getResult(self, seq_path):
        """Returns alignment from sequences.
        
        Currently does not allow parameter tuning of program and uses
        default parameters -- this is bad and should be fixed.

        #TODO: allow command-line access to important aln params.
        """
        module = self.Params['Module']
        seqs = self.getData(seq_path)
        result = module.align_unaligned_seqs(seqs, moltype=DNA)    
        #TODO: add params
        return result

    def __call__(self, result_path=None, log_path=None, *args, **kwargs):
        """Calls superclass method to align seqs"""
        return FunctionWithParams.__call__(self, result_path=result_path,
            log_path=log_path, *args, **kwargs)
            
class InfernalAligner(Aligner):
    Name = 'InfernalAligner'

    def __init__(self, params):
        """Return new InfernalAligner object with specified params.
        """
        _params = {
            'moltype': DNA,
            'Application': 'Infernal',
            }
        _params.update(params)
        Aligner.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None, \
        failure_path=None, cmbuild_params=None, cmalign_params=None):
        
        log_params = []
        # load candidate sequences
        candidate_sequences = dict(MinimalFastaParser(open(seq_path,'U')))
        
        # load template sequences
        info, template_alignment, struct = list(MinimalRfamParser(open(\
            self.Params['template_filepath'],'U'),\
            seq_constructor=ChangedSequence))[0]
        
        moltype = self.Params['moltype']
        
        #Need to make separate mapping for unaligned sequences
        unaligned = SequenceCollection(candidate_sequences,MolType=moltype)
        int_map, int_keys = unaligned.getIntMap(prefix='unaligned_')
        int_map = SequenceCollection(int_map,MolType=moltype)
        
        #Turn on --gapthresh option in cmbuild to force alignment to full model
        if cmbuild_params is None:
            cmbuild_params = {}
        cmbuild_params.update({'--gapthresh':1.0})
        
        #record cmbuild parameters
        log_params.append('cmbuild parameters:')
        log_params.append(str(cmbuild_params))
        
        #Turn on --sub option in Infernal, since we know the unaligned sequences
        # are fragments.
        #Also turn on --gapthresh to use same gapthresh as was used to build
        # model
        
        if cmalign_params is None:
            cmalign_params = {}
        cmalign_params.update({'--sub':True,'--gapthresh':1.0})
        
        #record cmalign parameters
        log_params.append('cmalign parameters:')
        log_params.append(str(cmalign_params))
        
        #Align sequences to alignment including alignment gaps.
        aligned, struct_string = cmalign_from_alignment(aln=template_alignment,\
            structure_string=struct,\
            seqs=int_map,\
            moltype=moltype,\
            include_aln=True,\
            params=cmalign_params,\
            cmbuild_params=cmbuild_params)
        
        #Pull out original sequences from full alignment.
        infernal_aligned={}
        aligned_dict = aligned.NamedSeqs
        for key in int_map.Names:
            infernal_aligned[int_keys.get(key,key)]=aligned_dict[key]
        
        #Create an Alignment object from alignment dict
        infernal_aligned = Alignment(infernal_aligned,MolType=moltype)
        
        if log_path is not None:
            log_file = open(log_path,'w')
            log_file.write('\n'.join(log_params))
            log_file.close()
        
        if result_path is not None:
            result_file = open(result_path,'w')
            result_file.write(infernal_aligned.toFasta())
            result_file.close()
            return None
        else:
            try:
                return infernal_aligned
            except ValueError:
                return {}


class PyNastAligner(Aligner):
    Name = 'PyNastAligner'

    def __init__(self, params):
        """Return new PyNastAligner object with specified params.
        """
        _params = {
            'min_pct': 75.0,
            'min_len': 1000,
            'blast_db': None,
            'template_filepath': None,
            'pairwise_alignment_method': 'pair_hmm',
            'Application': 'PyNAST',
            'Algorithm': 'NAST',
            }
        _params.update(params)
        Aligner.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None, 
                 failure_path=None):
        # load candidate sequences
        seq_file = open(seq_path, 'r')
        candidate_sequences = MinimalFastaParser(seq_file)

        # load template sequences
        template_alignment = LoadSeqs(
            self.Params['template_filepath'], moltype=DNA, 
            format='fasta', aligned=DenseAlignment)

        # initialize_logger
        logger = NastLogger(log_path)

        # get function for pairwise alignment method
        pairwise_alignment_fcn = pairwise_alignment_methods[
            self.Params['pairwise_alignment_method']]

        pynast_aligned, pynast_failed = pynast_seqs(
            candidate_sequences,
            template_alignment,
            min_pct=self.Params['min_pct'],
            min_len=self.Params['min_len'],
            blast_db=self.Params['blast_db'],
            align_unaligned_seqs_f=pairwise_alignment_fcn,
            logger=logger,
            )

        logger.record(str(self))

        if failure_path is not None:
            fail_file = open(failure_path,'w')
            for seq in pynast_failed:
                fail_file.write(seq.toFasta())
                fail_file.write('\n')
            fail_file.close()

        if result_path is not None:
            result_file = open(result_path,'w')
            for seq in pynast_aligned:
                result_file.write(seq.toFasta())
                result_file.write('\n')
            result_file.close()
            return None
        else:
            try:
                return LoadSeqs(data=pynast_aligned,aligned=DenseAlignment)
            except ValueError:
                return {}

usage_str = """usage: %prog [options] {-i INPUT_SEQUENCES_FILEPATH}

[] indicates optional input (order unimportant) 
{} indicates required input (order unimportant) 

Example usage:

Align 10_seq.fasta (-i) using muscle (-m). Output files will be stored in 
 ./muscle_aligned/
 python align_seqs.py -i 10_seq.fasta -m muscle
 
Align 10_seq.fasta (-i) with pynast (default) against template_aln.fasta (-t).
 Output files will be stored in ./pynast_aligned/
 python align_seqs.py -i 10_seq.fasta -t template_aln.fasta

Align 10_seq.fasta (-i) with infernal against template_aln.sto (-t).
 Output files will be stored in ./infernal_aligned/
 python align_seqs.py -i 10_seq.fasta -t template_aln.sto -m infernal
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i','--input_fasta_fp',\
        help='Input sequences to align [REQUIRED]')
          
    parser.add_option('-t','--template_fp',\
          type='string',dest='template_fp',help='Filepath for '+\
          'template against [default: %default; REQUIRED if -m pynast'+\
          'or -m infernal]')

    alignment_method_choices = \
     alignment_method_constructors.keys() + alignment_module_names.keys()
    parser.add_option('-m','--alignment_method',\
          type='choice',help='Method for aligning'+\
          ' sequences [default: %default]',choices=alignment_method_choices)
          
    pairwise_alignment_method_choices = pairwise_alignment_methods.keys()
    parser.add_option('-a','--pairwise_alignment_method',\
          type='choice',help='method for performing pairwise ' +\
          'alignment in PyNAST [default: %default]',\
          choices=pairwise_alignment_method_choices)

    parser.add_option('-d','--blast_db',\
          dest='blast_db',help='Database to blast against when -m pynast '+\
          '[default: created on-the-fly from template_alignment]')
          
    parser.add_option('-b','--blast_executable',\
        help='Path to blast executable when -m pynast [default: %default]')
          
    parser.add_option('-o','--output_dir',\
          help='Path to store '+\
          'result file [default: <ALIGNMENT_METHOD>_aligned]')
          
    parser.add_option('-e','--min_length',\
          type='int',help='Minimum sequence '+\
          'length to include in alignment [default: %default]')
          
    parser.add_option('-p','--min_percent_id',\
          type='float',help='Minimum percent '+\
          'sequence identity to closest blast hit to include sequence in'+\
          ' alignment [default: %default]')
          

    parser.set_defaults(verbose=False, alignment_method='pynast',\
     pairwise_alignment_method='blast', min_percent_id=75.0,min_length=1000,\
     blast_db=None,blast_executable=qiime_config['blastall_fp'])

    opts,args = parser.parse_args()
    
    required_options = ['input_fasta_fp']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
    
    if opts.min_percent_id <= 1.0:
        opts.min_percent_id *= 100
        
    if not (1.0 <= opts.min_percent_id <= 100.0):
        parser.error('Minimum percent sequence identity must be' +\
        ' between 1.0 and 100.0: %2.2f' % opts.min_percent_id)
        
    if not opts.template_fp and opts.alignment_method == 'pynast':
        parser.error('PyNAST requires a template alignment to be passed via -t')

    return opts,args


alignment_method_constructors ={'pynast':PyNastAligner,\
    'infernal':InfernalAligner}

pairwise_alignment_methods = {
    'muscle':muscle_align_unaligned_seqs,
    'mafft':mafft_align_unaligned_seqs,
    'clustal':clustal_align_unaligned_seqs,
    'blast':blast_align_unaligned_seqs,
    'pair_hmm':pair_hmm_align_unaligned_seqs,
}

alignment_module_names = {'muscle':cogent.app.muscle, 
    'clustalw':cogent.app.clustalw, 'mafft':cogent.app.mafft, \
    'infernal':cogent.app.infernal}


if __name__ == "__main__":
    opts,args = parse_command_line_parameters()

    input_seqs_filepath = opts.input_fasta_fp
    alignment_method = opts.alignment_method
    output_dir = opts.output_dir or alignment_method + '_aligned'
    try:
        makedirs(output_dir)
    except OSError:
        # output directory exists
        pass
    
    fpath, ext = splitext(input_seqs_filepath)
    input_dir, fname = split(fpath)
    
    result_path = output_dir + '/' + fname + "_aligned" + ext
    log_path = output_dir + '/' + fname + "_log.txt"
    failure_path = output_dir + '/' + fname + "_failures.fasta"
 
    try:
        # define the aligner params
        aligner_constructor =\
         alignment_method_constructors[alignment_method]
        aligner_type = alignment_method
        params = {'min_len':opts.min_length,\
                  'min_pct':opts.min_percent_id,\
                  'template_filepath':opts.template_fp,\
                  'blast_db':opts.blast_db,
                  'pairwise_alignment_method': opts.pairwise_alignment_method}
        # build the aligner object
        aligner = aligner_constructor(params)
        # apply the aligner
        aligner(input_seqs_filepath,result_path=result_path,\
         log_path=log_path,failure_path=failure_path)
        
    except KeyError:
        # define the aligner params
        aligner = CogentAligner({
        'Module':alignment_module_names[alignment_method],
        'Method':alignment_method
        })
        # build the aligner
        aligner_type = 'Cogent'
        # apply the aligner
        aligner(result_path, seq_path=input_seqs_filepath,
            log_path=log_path)

