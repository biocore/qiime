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

### notes on splitlib fastq options:
# barcode_read_fps: via Command
# store_qual_scores: via Command
# sample_ids: via Command
# store_demultiplexed_fastq: via Command
# retain_unassigned_reads: via Command (Failed == False, Sample == None)
# max_bad_run_length: via wf_quality, 
#       UNTESTED
# min_per_read_length_fraction: via wf_quality, if truncation happens, do
#       in place update on item
#       STUBBED OUT
# sequence_max_n: via wf_sequence
#        STUBBED OUT (ambiguous_count), UNTESTED
# start_seq_id: via Command, but also hopefully deprecated in favor of 
#       HDF5 format
# rev_comp_barcode: via Command and iterators? only if the barcodes are separate
#       then it is possible to do at the iterator level...
# rev_comp_mapping_barcodes: via Command
# rev_comp: via Command and iterators
# phred_quality_threshold: via wf_quality
#       STUBBED OUT, basically implemented? split_libraries_fastq is difficult to read...
# barcode_type: via wf_demultiplex
#       DONE
# max_barcode_error: via wf_demultiplex
# phred_offset: via Command and iterators

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

    @priority(1000)
    @no_requirements
    def wf_init(self, item):
        """Perform per sequence state initialization

        This workflow group will reset FinalState and will set the following in
        FinalState:

            Sequence
        """
        self._init_final_state(item)

    @priority(200)
    @requires(ValidData=_has_qual)
    def wf_quality(self, item):
        """Check sequence quality

        This workflow group may update _item_ in the event of a sequence
        truncation due to quality!

        """
        self._quality_max_bad_run_length(item)
        self._quality_min_per_read_length_fraction(item)

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
        self._demultiplex_max_barcode_error(item)

    ### should this be wf_instrument for instriument specific checks?
    @priority(100)
    @requires(Option='check_primer', Values=True)
    def wf_primer(self, item):
        """Perform primer validation

        Primer validation may update the following keys in FinalState:

            Sequence
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
        self._sequence_length_check(item)
        self._sequence_ambiguous_count(item)

    ### End Workflow methods

    ### Start quality methods

    @requires(Option='phred_quality_threshold')
    @requires(Option='max_bad_run_length')
    def _quality_max_bad_run_length(self, item):
        """Fail sequence if there is a poor quality run

        Warning: this method can modify item in place
        """
        max_bad_run_length = self.Options['max_bad_run_length']
        phred_quality_threshold = self.Options['phred_quality_threshold']

        # can cythonize
        run_length = 0
        max_run_length = 0
        run_start_idx = 0
        max_run_start_idx = 0

        for idx, v in enumerate(item['Qual']):
            if v <= phred_quality_threshold:
                max_run_length += 1
            else:
                if run_length > max_run_length:
                    max_run_length = run_length
                    max_run_start_idx = run_start_idx

                run_length = 0
                run_start_idx = idx

                if max_run_length == 0:
                    max_run_start_idx = run_start_idx

        if max_run_length > max_bad_run_length:
            item['Qual'] = item['Qual'][:max_run_start_idx+1]
            item['Sequence'] = item['Sequence'][:max_run_start_idx+1]
            self.Stats['_quality_max_bad_run_length'] += 1

    @requires(Option='phred_quality_threshold')
    @requires(Option='min_per_read_length_fraction')
    def _quality_min_per_read_length_fraction(self, item):
        """Fail a sequence if a percentage of bad quality calls exist"""
        bad_bases = item['Qual'] < self.Options['phred_quality_threshold']
        bad_bases_count = bad_bases.sum(dtype=float)
        threshold = 1 - self.Options['min_per_read_length_fraction']
        
        if (bad_bases_count / len(item['Sequence'])) > threshold:
            self.Failed = True
            self.Stats['min_per_read_length_fraction'] += 1

    ### End quality methods

    ### Start demultiplex methods
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

    #### use kwargs for method and bc_length
    def _demultiplex_encoded_barcode(self, item, method=decode_golay_12,
                                     bc_length=12):
        """Correct and decode an encoded barcode"""
        if item['Barcode'] is not None:
            putative_bc = item['Barcode']
        else:
            putative_bc = item['Sequence'][:bc_length]
            ### if this case happens, need to update item['Sequence'] to
            ### trim off the barcode!

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
            self.Stats['Unknown barcode'] += 1
        else:
            self.FinalState['Sample'] = sample
   
    ### really need the requires to be a nonnone value:
    # @requires(Option='max_barcode_error', Values=_not_none)
    @requires(Option='max_barcode_error')
    def _demultiplex_max_barcode_error(self, item):
        """ """
        bc_errors = self.Options['max_barcode_error']
        if self.FinalState['Corrected barcode errors'] > bc_errors:
            self.Failed = True
            self.Stats['exceeds_bc_errors'] += 1

    ### End demultiplex methods

    ### Start init methods

    def _init_final_state(self, item):
        """Reset per sequence state"""
        for k in self.FinalState:
            self.FinalState[k] = None
        self.FinalState['Sequence'] = item['Sequence']

    ### End init methods

    ### Start primer methods

    @requires(Option='instrument_type', Values='454')
    def _primer_instrument_454(self, item):
        """Check for a valid primer"""
        self._primer_check_forward(item)

    @requires(Option='retain_primer')
    @requires(Option='max_primer_mismatch')
    def _primer_check_forward(self, item):
        """Attempt to determine if the forward primer exists and trim if there
        
        Warning: this method may do an in place update on item if retain primer
        False.
        """
        seq = item['Sequence']
        qual = item['Qual']

        obs_barcode = self.FinalState['Final barcode']
        exp_primers = self.Primers.get(obs_barcode, None)

        if exp_primers is None:
            self.Stats['unknown_primer_barcode_pair'] += 1
            self.Failed = True
            return

        len_primer = len(exp_primers[0])

        obs_primer = seq[:len_primer]

        mm = array([_count_mismatches(obs_primer, p) for p in exp_primers])

        if (mm > self.Options['max_primer_mismatch']).all():
            self.Failed = True
            self.Stats['max_primer_mismatch'] += 1
            self.Stats['exceeds_max_primer_mismatch'] += 1
            return

        ### should decompose
        if not self.Options['retain_primer']:
            seq = seq[len_primer:]
            item['Sequence'] = seq
            if qual is not None:
                qual = qual[len_primer:]
                item['Qual'] = qual

        self.FinalState['Forward primer'] = obs_primer
        self.FinalState['Sequence'] = seq
        self.FinalState['Qual'] = qual

    ### End primer methods

    ### Start sequence methods

    @requires(Option='min_seq_len')
    def _sequence_length_check(self, item):
        """Checks minimum sequence length"""
        if len(item['Sequence']) < self.Options['min_seq_len']:
            self.Failed = True
            self.Stats['min_seq_len'] += 1

    @requires(Option='ambiguous_count')
    def _sequence_ambiguous_count(self, item):
        """Fail if the number of N characters is greater than threshold"""
        count = item['Sequence'].count('N')
        if count > self.Options['ambiguous_count']:
            self.Failed = True
            self.Stats['ambiguous_count'] += 1

    ### End sequence methods
