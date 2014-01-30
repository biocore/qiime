#!/usr/bin/env python

"""Filter poor quality reads, trim barcodes/primers and assign to samples"""

from qiime.workflow.core import Workflow, requires, priority, no_requirements
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser
from qiime.parse import MinimalQualParser
from itertools import chain, izip
from qiime.util import MetadataMap
from qiime.parse import is_casava_v180_or_later
from qiime.split_libraries import expand_degeneracies
from qiime.hamming import decode_barcode_8 as decode_hamming_8
from qiime.golay import decode as decode_golay_12
from qiime.quality import ascii_to_phred33, ascii_to_phred64
from numpy import array

# pre allocate the iterable return. This is done for performance reasons to
# avoid frequent reallocations and to ensure a consistent object type
_iter_prealloc = {'SequenceID':None,
                  'Sequence':None,
                  'Qual':None,
                  'Barcode':None}

def _reset_iter_prealloc():
    """Reset the buffer"""
    for k in _iter_prealloc:
        _iter_prealloc[k] = None

def _fasta_qual_gen(fasta_gen, qual_gen):
    """Yield fasta and qual together

    Raises ValueError if the sequence IDs and quality IDs are not in the same
    order. Raises ValueError if the sequence length does not match the length
    of the quality score.

    Note: object yielded is updated on each iteration. A new object is _not_
    created on each iteration. This is done for performance reasons, where
    quick testing showed a 50% reduction in runtime.
    """
    for (seq_id, seq), (qual_id, qual) in izip(fasta_gen, qual_gen):
        if seq_id != qual_id:
            raise ValueError("%s is not equal to %s!" % (seq_id, qual_id))
        if len(seq) != len(qual):
            raise ValueError("%s is not equal length to %s!" % (seq_id,
                                                                qual_id))

        _iter_prealloc['SequenceID'] = seq_id
        _iter_prealloc['Sequence'] = seq
        _iter_prealloc['Qual'] = qual

        yield _iter_prealloc

def _fasta_gen(fasta_gens):
    """Yield fasta data

    Note: object yielded is updated on each iteration. A new object is _not_
    created on each iteration. This is done for performance reasons, where
    quick testing showed a 50% reduction in runtime.
    """
    for id_, seq in fasta_gens:
        _iter_prealloc['SequenceID'] = id_
        _iter_prealloc['Sequence'] = seq
        yield _iter_prealloc

def fasta_iterator(fasta_fps, qual_fps=None):
    """Yield fasta and qual data

    Expects file-like objects. If qual_fps is not None, quality scores are
    yielded. The return will either be:

    {'SequenceID':foo, 'Sequence':bar, 'Qual':array([])}

    or

    {'SequenceID':foo, 'Sequence':bar, 'Qual':None}

    Note: object yielded is updated on each iteration. A new object is _not_
    created on each iteration. This is done for performance reasons, where
    quick testing showed a 50% reduction in runtime.
    """
    _reset_iter_prealloc()

    fasta_gens = chain(*[MinimalFastaParser(f) for f in fasta_fps])

    if qual_fps is not None:
        qual_gens = chain(*[MinimalQualParser(f) for f in qual_fps])
        gen = _fasta_qual_gen(fasta_gens, qual_gens)
    else:
        qual_gens = None
        gen = _fasta_gen(fasta_gens)

    return gen

def _fastq_barcode_gen(fastq_gens, barcode_gens, phred_f):
    """Yield fastq and barcode data

    Note: object yielded is updated on each iteration. A new object is _not_
    created on each iteration. This is done for performance reasons, where
    quick testing showed a 50% reduction in runtime.
    """
    _gen = izip(fastq_gens, barcode_gens)
    for (seqid, seq, qual), (bc_seqid, bc_seq, bc_qual) in _gen:
        if seqid != bc_seqid:
            raise ValueError("%s is not equal to %s!" % (seqid, bc_seqid))
        _iter_prealloc['SequenceID'] = seqid
        _iter_prealloc['Sequence'] = seq
        _iter_prealloc['Qual'] = array([phred_f(q) for q in qual])
        _iter_prealloc['Barcode'] = bc_seq

        yield _iter_prealloc

def _fastq_gen(fastq_gens, phred_f):
    """Yield fastq data

    Note: object yielded is updated on each iteration. A new object is _not_
    created on each iteration. This is done for performance reasons, where
    quick testing showed a 50% reduction in runtime.
    """
    for (seqid, seq, qual) in fastq_gens:
        _iter_prealloc['SequenceID'] = seqid
        _iter_prealloc['Sequence'] = seq
        _iter_prealloc['Qual'] = array([phred_f(q) for q in qual])

        yield _iter_prealloc

def fastq_iterator(fastq_fps, barcode_fps=None):
    """Yield fastq data

    Expects file-like objects. If barcode_fps is not None, barcodes are also
    yielded. The return will either be:

    {'SequenceID':foo, 'Sequence':bar, 'Qual':array([]), 'Barcode':foobar}

    or

    {'SequenceID':foo, 'Sequence':bar, 'Qual':array([]), 'Barcode':None}

    Note: object yielded is updated on each iteration. A new object is _not_
    created on each iteration. This is done for performance reasons, where
    quick testing showed a 50% reduction in runtime.
    """
    _reset_iter_prealloc()

    fastq_gens = chain(*[MinimalFastqParser(f) for f in fastq_fps])

    # peek
    first_item = fastq_gens.next()
    seqid, seq, qual = first_item
    fastq_gens = chain([first_item], fastq_gens)

    # from qiime.parse.parse_fastq_qual_score (v1.8.0)
    if is_casava_v180_or_later('@%s' % seqid):
        ascii_to_phred_f = ascii_to_phred33
    else:
        ascii_to_phred_f = ascii_to_phred64

    if barcode_fps:
        barcode_gens = chain(*[MinimalFastqParser(f) for f in barcode_fps])
        gen = _fastq_barcode_gen(fastq_gens, barcode_gens, ascii_to_phred_f)
    else:
        gen = _fastq_gen(fastq_gens, ascii_to_phred_f)

    return gen

### can cythonize
def _count_mismatches(seq1, seq2):
    """Counts mismatches between two sequences"""
    return sum([a != b for a, b in zip(seq1, seq2)])

def _has_qual(item):
    """Check if an item has Qual"""
    return item['Qual'] is not None

class SequenceWorkflow(Workflow):
    """Implement the sequence processing workflow
    
    All workflow methods expect an item that is dict-like with the following
    keys and value types:
        SequenceID : str
        Sequence   : str 
        Qual       : np.array or None
        Barcode    : str or None
    """
    FinalState = {'Forward primer':None,
                  'Reverse primer':None,
                  'Sequence':None,
                  'Qual':None,
                  'Sample':None,
                  'Original barcode':None,
                  'Corrected barcode':None,
                  'Final barcode':None,
                  'Corrected barcode errors':None}

    def _stage_state(self):
        """Fish out barcodes and primers from the mapping data"""
        bcs = {}
        primers = {}
        for sample in self.Mapping.SampleIds:
            sample_bc = self.Mapping.getCategoryValue(sample, 'BarcodeSequence')
            if sample_bc in bcs:
                raise ValueError("Duplicate barcode found for sample %s" \
                                 % sample)
            bcs[sample_bc] = sample

            sample_primers = self.Mapping.getCategoryValue(sample, 
                                                        'LinkerPrimerSequence')
            all_sample_primers = sample_primers.split(',')
            primers[sample_bc] = expand_degeneracies(all_sample_primers)

        self.Barcodes = bcs
        self.Primers = primers

    def _sanity_check(self):
        name = self.__class__.__name__
        if not hasattr(self, 'Mapping'):
            raise AttributeError("%s is missing Mapping!" % name)

        if not isinstance(self.Mapping, MetadataMap):
            raise AttributeError("self.Mapping is not of type MetadataMap")

    ### Start Workflow methods
    ### NEED TO ADD STRONG DEFINITIONS OF EXPECTED STATE CHANGES

    @priority(1000)
    @no_requirements
    def wf_init(self, item):
        """Perform per sequence state initialization
        
        This workflow group will reset FinalState. 
        """
        self._init_final_state(item)

    @priority(200)
    @requires(ValidData=_has_qual)
    def wf_quality(self, item):
        """Check sequence quality"""
        pass

    @priority(150)
    @requires(Option='demultiplex', Values=True)
    def wf_demultiplex(self, item):
        """Demultiplex a sequence

        If the sequence has not Failed, the following fields in FinalState will
        be set:

            Sample
            Original barcode
            Final barcode
        
        In addition, the following field may be set:

            Corrected barcode
            Corrected barcode errors
        
        This workflow group can trigger Failed and update Stats
        """
        self._demultiplex_golay12(item)
        self._demultiplex_hamming8(item)
        self._demultiplex_other(item)

        bc_errors = self.Options['max_bc_errors']
        if self.FinalState['Corrected barcode errors'] > bc_errors:
            self.Failed = True
            self.Stats['exceeds_bc_errors'] += 1

    @priority(100)
    @requires(Option='check_primer', Values=True)
    def wf_primer(self, item):
        """Perform primer validation

        Primer validation may update the following keys in FinalState:

            Forward primer
            Reverse primer

        This workflow group can trigger Failed and update Stats
        """
        self._primer_instrument_454(item)

    @priority(50)
    @no_requirements
    def wf_sequence(self, item):
        """Final sequence level checks
        
        Sequence level checks will not alter FinalState but may trigger Failed
        and update Stats
        """
        self._sequence_ambiguous_count(item)
        self._sequence_length_check(item)

    ### End Workflow methods

    @requires(Option='barcode_type', Values='golay_12')
    def _demultiplex_golay12(self, item):
        """Correct and decode a Golay 12nt barcode"""
        self._demultiplex_encoded_barcode(item, decode_golay_12, 12)

    @requires(Option='barcode_type', Values='hamming_8')
    def _demultiplex_hamming8(self, item):
        """Correct and decode a Hamming 8nt barcode"""
        self._demultiplex_encoded_barcode(item, decode_hamming_8, 8)

    @requires(Option='barcode_type', Values='variable')
    def _demultiplex_other(self, item):
        """Decode a variable length barcode"""
        raise NotImplementedError

    def _demultiplex_encoded_barcode(self, item, method, bc_length):
        """Correct and decode an encoded barcode"""
        if item['Barcode'] is not None:
            putative_bc = item['Barcode']
        else:
            putative_bc = item['Sequence'][:bc_length]

        self.FinalState['Original barcode'] = putative_bc

        if putative_bc in self.Barcodes:
            self.FinalState['Corrected barcode errors'] = 0
            final_bc = putative_bc
            sample = self.Barcodes[putative_bc]
        else:
            corrected, num_errors = method(putative_bc)
            final_bc = corrected
            self.FinalState['Corrected barcode'] = corrected
            self.FinalState['Corrected barcode errors'] = num_errors
            self.Stats['Barcodes corrected'] += 1
            sample = self.Barcodes.get(corrected, None)

        self.FinalState['Final barcode'] = final_bc

        if sample is None:
            self.Failed = True
        else:
            self.FinalState['Sample'] = sample

    def _init_final_state(self, item):
        """Reset per sequence state"""
        for k in self.FinalState:
            self.FinalState[k] = None

    @requires(Option='instrument_type', Values='454')
    def wf_check_primer(self, item):
        """Check for a valid primer"""
        self._count_primer_mismatches(item)
        #self._local_align_forward_primer(item)

    @requires(Option='max_primer_mismatch')
    def _primer_count_mismatches(self, item):
        """Assess primer mismatches"""
        seq = item['Sequence']
        qual = item['Qual']

        obs_barcode = self.FinalState['Final barcode']
        len_barcode = len(obs_barcode)

        exp_primers = self.Primers[obs_barcode]
        len_primer = len(exp_primers[0])

        obs_primer = seq[len_barcode:len_barcode + len_primer]

        mm = array([_count_mismatches(obs_primer, p) for p in exp_primers])

        if (mm > self.Options['max_primer_mismatch']).all():
            self.Failed = True
            self.Stats['max_primer_mismatch'] += 1
            self.Stats['exceeds_max_primer_mismatch'] += 1

        ### should decompose
        if not self.Options['retain_primer']:
            seq = seq[len_primer:]
            if qual is not None:
                qual = qual[len_primer:]

        self.FinalState['Forward primer'] = obs_primer
        self.FinalState['Sequence'] = seq

    ##### for truncating i believe, but isn't clear why we need to attempt to 
    ##### align against all possible primers instead of just the one we expect

    ### THIS IS STILL IN PROGRESS
    @requires(Option='local_align_forward_primer', Values=True)
    @requires(Option='max_primer_mismatch')
    def _primer_local_align_forward(self, item):
        """ """
        seq = item['Sequence']
        qual = item['Qual']
        
        failed = True
        max_primer_mismatch = self.Options['max_primer_mismatch']
        for primer in self._primers:
            mismatches, hit_start = local_align_primer_seq(primer, seq)
            if mismatches <= max_primer_mismatch:
                seq = seq[hit_start + len(primer):]
                qual = seq[hit_start + len(primer):]
                failed = False
                break

        if failed:
            self.Stats['max_primer_mismatch'] += 1
            self.Stats['exceeds_max_primer_mismatch'] = 1
        else:
            self.FinalState['Forward primer'] = primer
            self.FinalState['Sequence'] = seq
            self.FinalState['Qual'] = qual
    
    @requires(Option='min_seq_len')
    def _sequence_length_check(self, item):
        """Checks minimum sequence length"""
        if len(item['Sequence']) < self.Options['min_seq_len']:
            self.Failed = True
            self.Stats['min_seq_len'] += 1

