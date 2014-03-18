#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = [
    "Rob Knight",
    "Greg Caporaso",
    "Jeremy Widmann",
    "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


"""Contains code for aligning sequences, using several techniques.

This module has the responsibility for taking a set of sequences and
returning an alignment. Mostly, it will be thin wrappers for code
already in cogent.app.*, to which wrappers for e.g. PyNAST need to be
added..
"""
import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
from os import remove
from numpy import median
from cogent import LoadSeqs, DNA
from cogent.core.alignment import DenseAlignment, SequenceCollection, Alignment
from cogent.core.sequence import DnaSequence as Dna
from skbio.parse.sequences import parse_fasta
from skbio.core.exception import RecordError
from skbio.app.util import ApplicationNotFoundError
from cogent.app.infernal import cmalign_from_alignment
from cogent.parse.rfam import MinimalRfamParser, ChangedSequence
import cogent.app.clustalw
import cogent.app.mafft

from qiime.util import (get_tmp_filename,
                        FunctionWithParams,
                        get_qiime_temp_dir)
import cogent.app.muscle_v38


# Load PyNAST if it's available. If it's not, skip it if not but set up
# to raise errors if the user tries to use it.
try:
    from pynast.util import pynast_seqs, pairwise_alignment_methods
    from pynast.logger import NastLogger

except ImportError:
    def raise_pynast_not_found_error(*args, **kwargs):
        raise ApplicationNotFoundError("PyNAST cannot be found.\nIs PyNAST installed? Is it in your $PYTHONPATH?" +
                                       "\nYou can obtain PyNAST from http://qiime.org/pynast/.")
    # set functions which cannot be imported to raise_pynast_not_found_error
    pynast_seqs = NastLogger = raise_pynast_not_found_error
    pairwise_alignment_methods = {}


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

    def __call__(self, seq_path, result_path=None, log_path=None):
        """Returns alignment from sequences.

        Parameters:
        seq_path: path to file of sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path as fasta, otherwise should
        return cogent.core.alignment.DenseAlignment object.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError("Aligner is an abstract class")


class CogentAligner(Aligner):

    """Generic aligner using Cogent multiple alignment methods."""

    Name = 'CogentAligner'

    def getResult(self, seq_path):
        """Returns alignment from sequences.

        By convention, app parameters begin with a '-'.  Key-value
        pairs in self.Params following this convention will be passed
        as parameters to the module's alignment function.
        """
        module = self.Params['Module']
        seqs = self.getData(seq_path)
        params = dict(
            [(k, v) for (k, v) in self.Params.items() if k.startswith('-')])
        result = module.align_unaligned_seqs(seqs, moltype=DNA, params=params)
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

    def __call__(self, seq_path, result_path=None, log_path=None,
                 failure_path=None, cmbuild_params=None, cmalign_params=None):

        log_params = []
        # load candidate sequences
        candidate_sequences = dict(parse_fasta(open(seq_path, 'U')))

        # load template sequences
        try:
            info, template_alignment, struct = list(MinimalRfamParser(open(
                self.Params['template_filepath'], 'U'),
                seq_constructor=ChangedSequence))[0]
        except RecordError:
            raise ValueError(
                "Template alignment must be in Stockholm format with corresponding secondary structure annotation when using InfernalAligner.")

        moltype = self.Params['moltype']

        # Need to make separate mapping for unaligned sequences
        unaligned = SequenceCollection(candidate_sequences, MolType=moltype)
        int_map, int_keys = unaligned.getIntMap(prefix='unaligned_')
        int_map = SequenceCollection(int_map, MolType=moltype)

        # Turn on --gapthresh option in cmbuild to force alignment to full
        # model
        if cmbuild_params is None:
            cmbuild_params = {}
        cmbuild_params.update({'--gapthresh': 1.0})

        # record cmbuild parameters
        log_params.append('cmbuild parameters:')
        log_params.append(str(cmbuild_params))

        # Turn on --sub option in Infernal, since we know the unaligned sequences
        # are fragments.
        # Also turn on --gapthresh to use same gapthresh as was used to build
        # model

        if cmalign_params is None:
            cmalign_params = {}
        cmalign_params.update({'--sub': True, '--gapthresh': 1.0})

        # record cmalign parameters
        log_params.append('cmalign parameters:')
        log_params.append(str(cmalign_params))

        # Align sequences to alignment including alignment gaps.
        aligned, struct_string = cmalign_from_alignment(aln=template_alignment,
                                                        structure_string=struct,
                                                        seqs=int_map,
                                                        moltype=moltype,
                                                        include_aln=True,
                                                        params=cmalign_params,
                                                        cmbuild_params=cmbuild_params)

        # Pull out original sequences from full alignment.
        infernal_aligned = {}
        aligned_dict = aligned.NamedSeqs
        for key in int_map.Names:
            infernal_aligned[int_keys.get(key, key)] = aligned_dict[key]

        # Create an Alignment object from alignment dict
        infernal_aligned = Alignment(infernal_aligned, MolType=moltype)

        if log_path is not None:
            log_file = open(log_path, 'w')
            log_file.write('\n'.join(log_params))
            log_file.close()

        if result_path is not None:
            result_file = open(result_path, 'w')
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
            'min_len': 150,
            'blast_db': None,
            'template_filepath': None,
            'pairwise_alignment_method': 'blast',
            'Application': 'PyNAST',
            'Algorithm': 'NAST',
        }
        _params.update(params)
        Aligner.__init__(self, _params)

    def __call__(self, seq_path, result_path=None, log_path=None,
                 failure_path=None):
        # load candidate sequences
        seq_file = open(seq_path, 'U')
        candidate_sequences = parse_fasta(seq_file)

        # load template sequences
        template_alignment = []
        template_alignment_fp = self.Params['template_filepath']
        for seq_id, seq in parse_fasta(open(template_alignment_fp)):
            # replace '.' characters with '-' characters
            template_alignment.append((seq_id, seq.replace('.', '-').upper()))
        try:
            template_alignment = LoadSeqs(data=template_alignment, moltype=DNA,
                                          aligned=DenseAlignment)
        except KeyError as e:
            raise KeyError('Only ACGT-. characters can be contained in template alignments.' +
                           ' The offending character was: %s' % e)

        # initialize_logger
        logger = NastLogger(log_path)

        # get function for pairwise alignment method
        pairwise_alignment_f = pairwise_alignment_methods[
            self.Params['pairwise_alignment_method']]

        pynast_aligned, pynast_failed = pynast_seqs(
            candidate_sequences,
            template_alignment,
            min_pct=self.Params['min_pct'],
            min_len=self.Params['min_len'],
            align_unaligned_seqs_f=pairwise_alignment_f,
            logger=logger,
            temp_dir=get_qiime_temp_dir())

        logger.record(str(self))

        if failure_path is not None:
            fail_file = open(failure_path, 'w')
            for seq in pynast_failed:
                fail_file.write(seq.toFasta())
                fail_file.write('\n')
            fail_file.close()

        if result_path is not None:
            result_file = open(result_path, 'w')
            for seq in pynast_aligned:
                result_file.write(seq.toFasta())
                result_file.write('\n')
            result_file.close()
            return None
        else:
            try:
                return LoadSeqs(data=pynast_aligned, aligned=DenseAlignment)
            except ValueError:
                return {}


def compute_min_alignment_length(seqs_f, fraction=0.75):
    """ compute the min alignment length as n standard deviations below the mean """
    med_length = median([len(s) for _, s in parse_fasta(seqs_f)])
    return int(med_length * fraction)


alignment_method_constructors = {'pynast': PyNastAligner,
                                 'infernal': InfernalAligner}

alignment_module_names = {
    'muscle': cogent.app.muscle_v38,
    'clustalw': cogent.app.clustalw,
    'mafft': cogent.app.mafft,
    'infernal': cogent.app.infernal,
}
