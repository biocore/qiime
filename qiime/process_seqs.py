#!/usr/bin/env python

from qiime.workflow.core import Workflow, requires, priority, _continuous
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import MinimalQualParser
from itertools import chain, izip
from qiime.util import MetadataMap

from qiime.hamming import decode as decode_hamming_8
from qiime.golay import decode as decode_golay_12

def _fasta_qual_strict(fasta_gen, qual_gen):
    """Yield fasta and qual together

    Raises ValueError if the sequence IDs and quality IDs are not in the same
    order. Raises ValueError if the sequence length does not match the length 
    of the quality score. 
    """
    for (seq_id, seq), (qual_id, qual) in izip(fasta_gen, qual_gen):
        if seq_id != qual_id:
            raise ValueError("%s is not equal to %s!" % (seq_id, qual_id))
        if len(seq) != len(qual):
            raise ValueError("%s is not equal length to %s!" % (seq_id,qual_id))

        yield (seq_id, seq, qual)

def fasta_qual_iterator(fasta_fps, qual_fps=None):
    """Yield fasta and qual data

    Expects file-like objects. If qual_fps is not None, quality scores are 
    yielded, otherwise None is yielded for the quality. Specifically, the 
    tuple yielded is always of the form:

    (seq_id, seq, qual)
    """
    fasta_gens = chain(*map(MinimalFastaParser, fasta_fps))
    
    if qual_fps is not None:
        qual_gens = chain(*map(MinimalQualParser, qual_fps))
        gen = _fasta_qual_strict(fasta_gens, qual_gens)
    else:
        qual_gens = None
        gen = ((seq_id, seq, None) for seq_id, seq in fasta_gens)

    return gen

def _count_mismatches(seq1, seq2):
    """Counts mismatches between two sequences"""
    return sum([a == b for a,b in zip(seq1, seq2)])    

SEQ_ID_INDEX = 0
SEQ_INDEX = 1
QUAL_INDEX = 2

class SequenceWorkflow(Workflow):
    FinalState = {'fwd_primer':None,
                  'rev_primer':None,
                  'seq':None,
                  'qual':None,
                  'sample':None,
                  'original_barcode':None,
                  'corrected_barcode':None,
                  'final_barcode':None,
                  'corrected_barcode_errors':None}
    
    def _stage_state(self):
        """Fish out barcodes from the mapping data"""
        # set all the barcodes
        bcs = {}
        for sample in self.Mapping.SampleIds:
            sample_bc = self.Mapping.getCategoryValue(sample, 'barcode')
            if sample_bc in bcs:
                raise ValueError("Duplicate barcode found for sample %s" \
                                 % sample)
            else:
                bcs[sample_bc] = sample
        self.Barcodes = frozenset(bcs)

    def _sanity_check(self):
        name = self.__name__
        if not hasattr(self, 'Mapping'):
            raise AttributeError("%s is missing Mapping!" % name)
        
        if not isinstance(self.Mapping, MetadataMap):
            raise AttributeError("self.Mapping is not of type MetadataMap")

    ### Start Workflow methods

    @priority(1000)
    @no_requirements
    def wf_init(self, item):
        self._init_final_state(item)

    @priority(900)
    @requires(Option='max_bc_errors')
    @requires(Option='barcode_type', Values=['hamming_8','golay_12'])
    def wf_demultiplex_fixed(self, item):
        self._correct_golay12(item)
        self._correct_hamming8(item)

        bc_errors = self.Options['max_bc_errors']
        if self.FinalState['corrected_barcode_errors'] > bc_errors:
            self.Failed = True
            self.Stats['exceeds_bc_errors'] += 1
    
    @priority(900)
    @requires(Option='barcode_type', Values='variable')
    def wf_demultiplex_variable(self, item):
        raise NotImplementedError("variable length barcodes not supported yet")
        self._correct_variable(item)

    @priority(90)
    @requires(Option='min_seq_len')
    def wf_length_check(self, item):
        """Checks minimum sequence length"""
        seq_id, seq, qual_id, qual = item

        if len(seq) < self.Options['min_seq_len']:
            self.Failed = True
            self.Stats['min_seq_len'] += 1

    @priority(89)
    @requires(Option='instrument-type', Values='454')
    @requires(Option='disable_primer_check', Values=False)
    def wf_check_primer(self, item):
        """ """
        self._count_mismatches(item)
        self._local_align_forward_primer(item)

    ### End Workflow methods

    def _check_exact_barcode(self):
        """Check for a match"""
        return self.FinalState['original_barcode'] in self.Barcodes

    @requires(Option='barcode_type', Values='golay_12')
    def _correct_golay12(self, item):
        """ """
        self._correct_encoded_barcode(item, decode_golay_12, 12)

    @requires(Option='barcode_type', Values='hamming_8')
    def _correct_hamming8(self, item):
        """ """
        self._correct_encoded_barcode(item, decode_hamming_8, 8)

    def _correct_encoded_barcode(self, item, method, bc_length):
        putative_bc = item[SEQ_INDEX][:bc_length]
        self.FinalState['original_barcode'] = putative_bc
        
        if self._check_exact_barcode():
            self.FinalState['corrected_barcode_errors'] = 0
            final_bc = putative_bc
            sample = self.Barcodes.get(putative_bc, None)
        else:
            corrected, num_errors = method(putative_bc)
            final_bc = corrected

            self.FinalState['corrected_barcode'] = corrected
            self.FinalState['corrected_barcode_errors'] = num_errors
            self.Stats['barcodes_corrected'] += 1
            sample = self.Barcodes.get(corrected, None)

        self.FinalState['final_barcode'] = final_bc

        if sample is None:
            self.Failed = True
        else:
            self.FinalState['sample'] = sample

    def _init_final_state(self, item):
        """Reset final state"""
        for k in self.FinalState:
            self.FinalState[k] = None
        
    ##### the requires are likely wrong here
    @requires(Option='max_primer_mismatch')
    def _count_mismatches(self, item):
        """ """
        seq = item[SEQ_INDEX]
        qual = item[QUAL_INDEX]
        
        exp_primer = self.Mapping.getCategoryValue(self.FinalState['sample'],
                                                   'LinkerPrimerSequence'))
        len_primer = len(exp_primer)
        obs_primer = seq[:len_primer]
        
        mismatches = _count_mismatches(obs_primer, exp_primer)
        
        if not self.Options['retain_primer']:
            seq = seq[len_primer:]
            qual = qual[len_primer:]

        if mismatches > self.Options['max_primer_mismatch']:
            self.Failed = True
            self.Stats['max_primer_mismatch'] += 1
            self.Stats['exceeds_max_primer_mismatch'] = 1
        
        self.FinalState['fwd_primer'] = obs_primer
        self.FinalState['seq'] = seq

    ##### for truncating i believe, but isn't clear why we need to attempt to 
    ##### align against all possible primers instead of just the one we expect

    ### THIS IS STILL IN PROGRESS
    @requires(Option='local_align_forward_primer', Values=True)
    @requires(Option='max_primer_mismatch')
    def _local_align_forward_primer(self, item):
        """ """
        seq = item[SEQ_INDEX]
        qual = item[QUAL_INDEX]
        
        failed = True
        max_primer_mismatch = self.Options['max_primer_mismatch']
        for primer in self._primers:
            mismatches, hit_start = local_align_primer_seq(primer, fasta_seq)
            if mismatches <= max_primer_mismatch:
                seq = seq[hit_start + len(primer):]
                qual = seq[hit_start + len(primer):]
                failed = False
                break

        if failed:
            self.Stats['max_primer_mismatch'] += 1
            self.Stats['exceeds_max_primer_mismatch'] = 1
        else:
            self.FinalState['fwd_primer'] = primer
            self.FinalState['seq'] = seq
            self.FinalState['qual'] = qual


