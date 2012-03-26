#!/usr/bin/env python
"""Application controller for pplacer 1.1"""

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Kyle Bittinger","Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.6.0.dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Prototype"

from cogent.app.parameters import ValuedParameter, FlagParameter
from cogent.app.util import CommandLineApplication, FilePath, system, \
       CommandLineAppResult, ResultPath, remove, ApplicationError, \
       get_tmp_filename
from cogent.core.alignment import Alignment
from cogent.app.guppy import build_tree_from_json_using_params
from os.path import splitext,abspath,join,split
from StringIO import StringIO
from cogent.parse.phylip import get_align_for_phylip
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode

class Pplacer(CommandLineApplication):
    """pplacer Application Controller
    """

    _command = 'pplacer'
    _input_handler = '_input_as_multiline_string'
    _parameters = {
        # -c Specify the path to the reference package.
        '-c': ValuedParameter('-', Name='c', Delimiter=' ', IsPath=True),

        # -t Specify the reference tree filename.
        '-t': ValuedParameter('-', Name='t', Delimiter=' ', IsPath=True),

        # -r Specify the reference alignment filename.
        '-r': ValuedParameter('-', Name='r', Delimiter=' ', IsPath=True),

        # -s Supply a phyml stats.txt or a RAxML info file giving the model parameters.
        '-s': ValuedParameter('-', Name='s', Delimiter=' ', IsPath=True),

        # -d Specify the directory containing the reference information.
        '-d': ValuedParameter('-', Name='d', Delimiter=' ', IsPath=True),

        # -p Calculate posterior probabilities.
        '-p': FlagParameter('-', Name='p'),

        # -m Substitution model. Protein: are LG, WAG, or JTT. Nucleotides: GTR.
        '-m': ValuedParameter('-', Name='m', Delimiter=' '),

        # --model-freqs Use model frequencies instead of reference alignment frequencies.
        '--model-freqs': FlagParameter('--', Name='model-freqs'),

        # --gamma-cats Number of categories for discrete gamma model.
        '--gamma-cats': ValuedParameter('--', Name='gamma-cats', Delimiter=' '),

        # --gamma-alpha Specify the shape parameter for a discrete gamma model.
        '--gamma-alpha': ValuedParameter('--', Name='gamma-alpha', Delimiter=' '),

        # --ml-tolerance 1st stage branch len optimization tolerance (2nd stage to 1e-5). Default: 0.01.
        '--ml-tolerance': ValuedParameter('--', Name='ml-tolerance', Delimiter=' '),

        # --pp-rel-err Relative error for the posterior probability calculation. Default is 0.01.
        '--pp-rel-err': ValuedParameter('--', Name='pp-rel-err', Delimiter=' '),

        # --unif-prior Use a uniform prior rather than exponential.
        '--unif-prior': FlagParameter('--', Name='unif-prior'),

        # --start-pend Starting pendant branch length. Default is 0.1.
        '--start-pend': ValuedParameter('--', Name='start-pend', Delimiter=' '),
        
        # --max-pend Set the maximum ML pendant branch length. Default is 2.
        '--max-pend': ValuedParameter('--', Name='max-pend', Delimiter=' '),
        
        # --max-strikes Maximum number of strikes for baseball. 0 -> no ball playing. Default is 6.
        '--max-strikes': ValuedParameter('--', Name='max-strikes', Delimiter=' '),
        
        # --strike-box Set the size of the strike box in log likelihood units. Default is 3.
        '--strike-box': ValuedParameter('--', Name='strike-box', Delimiter=' '),
        
        # --max-pitches Set the maximum number of pitches for baseball. Default is 40.
        '--max-pitches': ValuedParameter('--', Name='max-pitches', Delimiter=' '),
        
        # --fantasy Desired likelihood cutoff for fantasy baseball mode. 0 -> no fantasy.
        '--fantasy': ValuedParameter('--', Name='fantasy', Delimiter=' '),
        
        # --fantasy-frac Fraction of fragments to use when running fantasy baseball. Default is 0.1.
        '--fantasy-frac': ValuedParameter('--', Name='fantasy-frac', Delimiter=' '),
        
        # --write-masked Write alignment masked to the region without gaps in the query.
        '--write-masked': FlagParameter('--', Name='write-masked'),
        
        # --verbosity Set verbosity level. 0 is silent, and 2 is quite a lot. Default is 1.
        '--verbosity': ValuedParameter('--', Name='verbosity', Delimiter=' '),
        
        # --unfriendly Do not run friend finder pre-analysis.
        '--unfriendly': FlagParameter('--', Name='unfriendly'),
        
        # --out-dir Specify the directory to write place files to.
        '--out-dir': ValuedParameter('--', Name='out-dir', Delimiter=' ', IsPath=True),
        
        # --pretend Only check out the files then report. Do not run the analysis.
        '--pretend': FlagParameter('--', Name='pretend'),

        # --csv Make a CSV file with the results.
        '--csv': FlagParameter('--', Name='csv'),

        # --old-format Make an old-format placefile with the resuls.
        '--old-format': FlagParameter('--', Name='old-format'),

        # --diagnostic Write file describing the 'diagnostic' mutations for various clades.
        '--diagnostic': FlagParameter('--', Name='diagnostic'),

        # --check-like Write out the likelihood of the reference tree, calculated two ways.
        '--check-like': FlagParameter('--', Name='check-like'),

        # --version Write out the version number and exit.
        '--version': FlagParameter('--', Name='version'),

        # --help  Display this list of options
        '--help': FlagParameter('--', Name='help'),
    }
 
    def getTmpFilename(self, tmp_dir="/tmp",prefix='tmp',suffix='.fasta',\
           include_class_id=False,result_constructor=FilePath):
        """ Define Tmp filename to contain .fasta suffix, since pplacer requires
            the suffix to be .fasta """
            
        return super(Pplacer,self).getTmpFilename(tmp_dir=tmp_dir,
                                    prefix=prefix,
                                    suffix=suffix,
                                    include_class_id=include_class_id,
                                    result_constructor=result_constructor)
    
    def _handle_app_result_build_failure(self,out,err,exit_status,result_paths):
        """ Catch the error when files are not produced """
        raise ApplicationError, \
         'Pplacer failed to produce an output file due to the following error: \n\n%s ' \
         % out.read()

    def _get_result_paths(self,data):
        """ Define the output filepaths """
        output_dir = self.Parameters['--out-dir'].Value
        result = {}
        result['json'] = ResultPath(Path=join(output_dir,
                                splitext(split(self._input_filename)[-1])[0] + \
                                '.jplace'))
        return result

def insert_sequences_into_tree(aln, moltype, params={},
                                           write_log=True):
    """Returns a tree from Alignment object aln.

    aln: an xxx.Alignment object, or data that can be used to build one.

    moltype: cogent.core.moltype.MolType object

    params: dict of parameters to pass in to the RAxML app controller.

    The result will be an xxx.Alignment object, or None if tree fails.
    """

    # convert aln to phy since seq_names need fixed to run through pplacer
    
    new_aln=get_align_for_phylip(StringIO(aln))

    # convert aln to fasta in case it is not already a fasta file
    aln2 = Alignment(new_aln)
    seqs = aln2.toFasta()

    ih = '_input_as_multiline_string'    

    pplacer_app = Pplacer(params=params,
                      InputHandler=ih,
                      WorkingDir=None,
                      SuppressStderr=False,
                      SuppressStdout=False)
    
    pplacer_result = pplacer_app(seqs)

    # write a log file
    if write_log:
        log_fp = join(params["--out-dir"],'log_pplacer_' + \
                      split(get_tmp_filename())[-1])
        log_file=open(log_fp,'w')
        log_file.write(pplacer_result['StdOut'].read())
        log_file.close()
        
    # use guppy to convert json file into a placement tree
    guppy_params={'tog':None}
    
    new_tree=build_tree_from_json_using_params(pplacer_result['json'].name, \
                                               output_dir=params['--out-dir'], \
                                               params=guppy_params)

    pplacer_result.cleanUp()
    
    return new_tree



    
