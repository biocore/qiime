#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Rob Knight","Greg Caporaso"]
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
from os import remove
from os.path import exists
from commands import getoutput
from optparse import OptionParser
from cogent import LoadSeqs, DNA
from cogent.core.alignment import DenseAlignment
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.util import get_tmp_filename, ApplicationNotFoundError
from qiime.util import FunctionWithParams
#app controllers that implement align_unaligned_seqs
import cogent.app.muscle
import cogent.app.clustalw
import cogent.app.mafft

# Load PyNAST if it's available. If it's not, skip it if not but set up
# to raise errors if the user tries to use it.
try:
    from pynast.util import pynast_seqs, classic_align_unaligned_seqs,\
        muscle_align_unaligned_seqs, mafft_align_unaligned_seqs,\
        clustal_align_unaligned_seqs, blast_align_unaligned_seqs
    from pynast.logger import NastLogger

except ImportError:
    def raise_pynast_not_found_error(*args, **kwargs):
        raise ApplicationNotFoundError,\
         "PyNAST not installed - pynast aligner currently unavailable." 
    # set functions which cannot be imported to raise_pynast_not_found_error
    pynast_seqs = NastLogger = classic_align_unaligned_seqs = \
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
            'pairwise_alignment_method': 'classic',
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

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage =\
     'usage: %prog [options] input_sequences_filepath'
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-m','--alignment_method',action='store',\
          type='string',dest='alignment_method',help='Method for aligning'+\
          ' sequences [default: %default]')
          
    parser.add_option('-a','--pairwise_alignment_method',action='store',\
          type='string',dest='pairwise_alignment_method',help='method '+\
          'for performing pairwise alignment in PyNAST ' +\
          '[default: %default]')

    parser.add_option('-d','--blast_db',action='store',\
          type='string',dest='blast_db',help='Database to blast'+\
          ' against [default: %default]')
          
    parser.add_option('-b','--blast_executable',action='store',\
        type='string',
        default='/home/caporaso/bin/blastall',
        help='Path to blast executable [default: %default]')
          
    parser.add_option('-t','--template_fp',action='store',\
          type='string',dest='template_fp',help='Filepath for '+\
          'template against [default: %default]')
          
    parser.add_option('-o','--result_fp',action='store',\
          type='string',dest='result_fp',help='Path to store '+\
          'result file [default: <input_sequences_filename>_aligned.fasta]')
          
    parser.add_option('-l','--log_fp',action='store',\
          type='string',dest='log_fp',help='Path to store '+\
          'log file [default: No log file created.]')
          
    parser.add_option('-e','--min_length',action='store',\
          type='int',dest='min_length',help='Minimum sequence '+\
          'length to include in alignment [default: %default]')
          
    parser.add_option('-p','--min_percent_id',action='store',\
          type='float',dest='min_percent_id',help='Minimum percent '+\
          'sequence identity to closest blast hit to include sequence in'+\
          ' alignment [default: %default]')
          
    parser.add_option('-f','--fail_fp',action='store',\
          type='string',dest='fail_fp',help='Path to store '+\
          'seqs which fail to align [default: No fail file created.]')

    parser.set_defaults(verbose=False, alignment_method='muscle',\
     pairwise_alignment_method='classic', min_percent_id=75.0,min_length=1000,\
     blast_db=None)

    opts,args = parser.parse_args()
    
    num_args = 1
    if len(args) != num_args:
       parser.error('Exactly one argument is required.')
       
    if not (opts.alignment_method in alignment_method_constructors or
            opts.alignment_method in alignment_module_names):
        parser.error(\
         'Invalid alignment method: %s.\nValid choices are: %s'\
         % (opts.alignment_method,\
            ' '.join(alignment_method_constructors.keys() +
                alignment_module_names.keys())))

    if opts.pairwise_alignment_method not in pairwise_alignment_methods:
        parser.error(
            'Invalid pairwise alignment method: %s.\nValid choices are: %s'\
                % (opts.pairwise_alignment_method,
                   ' '.join(pairwise_alignment_methods.keys())))
            
    if opts.min_percent_id <= 1.0:
        opts.min_percent_id *= 100
        
    if not (1.0 <= opts.min_percent_id <= 100.0):
        parser.error('Minimum percent sequence identity must be' +\
        ' between 1.0 and 100.0: %2.2f' % opts.min_percent_id)
        
    if not opts.template_fp and opts.alignment_method == 'pynast':
        parser.error('PyNAST requires a template alignment to be passed via -t')

    return opts,args


alignment_method_constructors = dict([\
        ('pynast',PyNastAligner)])

pairwise_alignment_methods = {
    'muscle':muscle_align_unaligned_seqs,
    'mafft':mafft_align_unaligned_seqs,
    'clustal':clustal_align_unaligned_seqs,
    'blast':blast_align_unaligned_seqs,
    'classic':classic_align_unaligned_seqs,
}

alignment_module_names = {'muscle':cogent.app.muscle, 
    'clustalw':cogent.app.clustalw, 'mafft':cogent.app.mafft}


if __name__ == "__main__":
    opts,args = parse_command_line_parameters()

    input_seqs_filepath = args[0]
   
    result_path = opts.result_fp or\
     input_seqs_filepath.replace('.fasta','_aligned.fasta')
     
    log_path = opts.log_fp
 
    try:
        # define the aligner params
        aligner_constructor =\
         alignment_method_constructors[opts.alignment_method]
        aligner_type = opts.alignment_method
        params = {'min_len':opts.min_length,\
                  'min_pct':opts.min_percent_id,\
                  'template_filepath':opts.template_fp,\
                  'blast_db':opts.blast_db,
                  'pairwise_alignment_method': opts.pairwise_alignment_method}
        failure_path = opts.fail_fp
        # build the aligner object
        aligner = aligner_constructor(params)
        # apply the aligner
        aligner(input_seqs_filepath,result_path=result_path,\
         log_path=log_path,failure_path=failure_path)
        
    except KeyError:
        # define the aligner params
        aligner = CogentAligner({
        'Module':alignment_module_names[opts.alignment_method],
        'Method':opts.alignment_method
        })
        # build the aligner
        aligner_type = 'Cogent'
        # apply the aligner
        aligner(result_path, seq_path=input_seqs_filepath,
            log_path=log_path)

