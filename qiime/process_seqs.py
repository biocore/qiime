#!/usr/bin/env python

"""Filter poor quality reads, trim barcodes/primers and assign to samples"""

import numpy as np

from collections import Counter
from itertools import izip

from skbio.factory.sequence import factory
from skbio.core.workflow import Workflow, requires, method, not_none
from qiime.hamming import decode_barcode_8 as decode_hamming_8
from qiime.golay import decode as decode_golay_12


def count_mismatches(seq1, seq2):
    """Counts mismatches between two sequences"""
    return sum(a != b for a, b in izip(seq1, seq2))


def has_sequence_qual(state):
    """Check if state has Qual"""
    return state['Qual'] is not None


def has_barcode_qual(state):
    """Check if state has Barcode Qual"""
    return state['Barcode Qual'] is not None


### notes on splitlib fastq options:
# barcode_read_fps: via Command
# store_qual_scores: via Command
# sample_ids: via Command
# store_demultiplexed_fastq: via Command
# retain_unassigned_reads: via Command (Failed == False, Sample == None)
# start_seq_id: via Command, but also hopefully deprecated in favor of
#       HDF5 format
# rev_comp_barcode: via Command and iterators? only if the barcodes are separate
#       then it is possible to do at the iterator level...
# rev_comp_mapping_barcodes: via Command
# rev_comp: via Command and iterators
# phred_offset: via Command and iterators


class Seqs(object):
    """Augmented sequence iterators

    This sequence iterator allows for optionally combining sequence reads with
    barcode data, as well as performing transforms independently on the reads
    or the barcode data. Barcode quality, if available, is also yielded.

    Attributes
    ----------
    reads
    barcodes

    Examples
    --------
    >>> import os
    >>> out = open('test_barcodes.fna', 'w')
    >>> out.write(">s1\nAT\n>s2\nGC\n")
    >>> out.close()
    >>> out = open('test_seqs.fq', 'w')
    >>> out.write("@s1\nAAAT\n+\nghgh\n@s2\nTTGG\n+\nfggh\n")
    >>> outgz.close()

    >>> from qiime.process_seqs import Seqs
    >>> it = Seqs(seq='test_seqs.fq', barcode='test_barcodes.fna')
    >>> for rec in it:
        ...     print rec['SequenceID']
        ...     print rec['Sequence']
        ...     print rec['Qual']
        ...     print rec['BarcodeID']
        ...     print rec['Barcode']
        ...     print rec['BarcodeQual']
    s1
    AAAT
    [39 40 39 40]
    s1
    AT
    None
    s2
    TTGG
    [38 39 39 40]
    s2
    GC
    None
    >>> os.remove('test_seqs.fq')
    >>> os.remove('test_barcodes.fna')

    """

    def __init__(self, seq, qual=None, barcode=None, barcode_qual=None,
                 seq_kwargs=None, barcode_kwargs=None):

        seq_kwargs = {} if seq_kwargs is None else seq_kwargs
        self.reads = factory(seq=seq, qual=qual, **seq_kwargs)

        if barcode is None:
            self.barcodes = None
        else:
            barcode_kwargs = {} if barcode_kwargs is None else barcode_kwargs
            self.barcodes = factory(seq=barcode, qual=barcode_qual,
                                    **barcode_kwargs)

    def __iter__(self):
        remap = (('SequenceID', 'BarcodeID'),
                 ('Sequence', 'Barcode'),
                 ('QualID', 'BarcodeQualID'),
                 ('Qual', 'BarcodeQual'))

        rec = {'SequenceID': None,
               'Sequence': None,
               'QualID': None,
               'Qual': None,
               'BarcodeID': None,
               'Barcode': None,
               'BarcodeQualID': None,
               'BarcodeQual': None}

        if self.barcodes is None:
            for seq in self.reads:
                rec.update(seq)
                yield rec
        else:
            for seq, barcode in izip(self.reads, self.barcodes):
                rec.update(seq)
                rec.update({new_k: barcode[old_k] for new_k, old_k in remap})

                if rec['SequenceID'] != rec['BarcodeID']:
                    raise ValueError("ID mismatch. SequenceID: %s, "
                                     "BarcodeID: %s" % (rec['SequenceID'],
                                                        rec['BarcodeID']))

                yield rec


class SequenceWorkflow(Workflow):
    """Implement the sequence processing workflow

    The sequence processing workflow manages the following tasks, executed in
    the following order::

        1. Quality filtering and trimming of primary sequence data
        2. Demultiplexing and assigning reads to samples
        3. Validating primers
        4. Sequence level quality checks (e.g., ambiguous bases)

    Execution of a task will only happen if it is relevant for the data. For
    instance, quality checks are only performed if the data being operated on
    has quality scores associated. Runtime control through options are also
    supported, such that, for instance, the Golay decoder is only executed if
    indicated by the options passed to the `SequenceWorkflow`.

    Any task can trigger `failed` and update `stats`.

    Parameters
    ----------
    options : dict
        Runtime options. See ``Options`` for more details
    barcodes : dict
        Mapping of barcode nucleotide sequence to a sample ID
    primers : dict
        Mapping of nucleotide sequence to enumerated possible primers

    Options
    -------
    ## DESCRIBE EACH OPTION THAT CAN AFFECT WHAT METHODS ARE EXECUTED

    Attributes
    ----------
    state
    stats
    options
    barcodes
    primers

   """

    def __init__(self, *args, **kwargs):
        if 'barcodes' not in kwargs:
            kwargs['barcodes'] = {}
        if 'primers' not in kwargs:
            kwargs['primers'] = {}

        kwargs['state'] = {'Forward primer': None,
                           'Reverse primer': None,
                           'Sequence': None,
                           'Qual': None,
                           'Barcode': None,
                           'Barcode Qual': None,
                           'Sample': None,
                           'Original barcode': None,
                           'Corrected barcode': None,
                           'Final barcode': None,
                           'Corrected barcode errors': None}

        kwargs['stats'] = Counter()

        super(SequenceWorkflow, self).__init__(self, *args, **kwargs)

    def initialize_state(self, item):
        """Reset `state` and update with the current `item`

        Parameters
        ----------
        item : dict
            An item from the `Workflow` generator
        """
        for k in self.state:
            self.state[k] = None
        self.state.update(item)

    ### Start Workflow methods

    @method(priority=200)
    @requires(state=has_qual)
    def wf_quality(self):
        """Check sequence quality

        Changes to `state`
        ------------------
        This workflow group may trim `state['Sequence']` and `state['Qual']` if
        quality trimming is enabled.

        Triggers for `failed`
        ---------------------
        - If to many nucleotides in `Sequence` are of poor quality.

        Impacted `stats`
        ----------------
        quality_max_bad_run_length
            Incremented if the read contained a run of poor quality bases
        min_per_read_length_fraction
            Incrememted if to many positions in `Sequence` are of poor quality
        """
        self._quality_max_bad_run_length()
        self._quality_min_per_read_length_fraction()

    @method(priority=150)
    @requires(option='demultiplex', values=True)
    def wf_demultiplex(self):
        """Demultiplex a sequence

        Changes to `state`
        ------------------
        Sample
        Original barcode
        Final barcode
        Barcode errors

        Triggers for `failed`
        ---------------------
        - If a sequence could not be associated to a sample
        - If the number of errors observed in the barcode exceed tolerance

        Impacted `stats`
        ----------------
        barcode_corrected
            Incremented if a barcode was corrected
        unknown_barcode
            Incremented if an unknown barcode was observed
        exceed_barcode_error
            Incremented if the number of observed barcode
            errors exceeded tolerance
        """
        self._demultiplex_golay12()
        self._demultiplex_hamming8()
        self._demultiplex_other()
        self._demultiplex_max_barcode_error()

    ### should this be wf_instrument for instriument specific checks?
    @method(priority=100)
    @requires(option='check_primer', values=True)
    def wf_primer(self):
        """Perform primer validation

        Changes to `state`
        ------------------
        Sequence
        Qual
        Forward primer
        Reverse primer

        Triggers for `failed`
        ---------------------
        - If the `primer` mapping does not contain primers associated with the
            nucleotide barcode
        Impacted `stats`
        ----------------
        unknown_primer_barcode_pair
        """
        self._primer_instrument_454()

    @method(priority=50)
    def wf_sequence(self):
        """Final sequence level checks

        Sequence level checks will not alter `state` but may trigger Failed
        and update `stats`

        Changes to `state`
        ------------------

        Triggers for `failed`
        ---------------------

        Impacted `stats`
        ----------------
        """
        self._sequence_length_check()
        self._sequence_ambiguous_count()

    ### End Workflow methods

    ### Start quality methods

    @requires(option='phred_quality_threshold')
    @requires(option='max_bad_run_length')
    def _quality_max_bad_run_length(self):
        """Fail sequence if there is a poor quality run"""
        max_bad_run_length = self.options['max_bad_run_length']
        phred_quality_threshold = self.options['phred_quality_threshold']

        # can cythonize
        run_length = 0
        max_run_length = 0
        run_start_idx = 0
        max_run_start_idx = 0

        for idx, v in enumerate(self.state['Qual']):
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
            self.state['Qual'] = self.state['Qual'][:max_run_start_idx+1]
            self.state['Sequence'] = self.state['Sequence'][:max_run_start_idx+1]
            self.stats['_quality_max_bad_run_length'] += 1

    @requires(option='phred_quality_threshold')
    @requires(option='min_per_read_length_fraction')
    def _quality_min_per_read_length_fraction(self):
        """Fail a sequence if a percentage of bad quality calls exist"""
        bad_bases = self.state['Qual'] < self.options['phred_quality_threshold']
        bad_bases_count = bad_bases.sum(dtype=float)
        threshold = 1 - self.options['min_per_read_length_fraction']

        if (bad_bases_count / len(self.state['Sequence'])) > threshold:
            self.failed = True
            self.stats['min_per_read_length_fraction'] += 1

    ### End quality methods

    ### Start demultiplex methods
    @requires(option='barcode_type', values='golay_12')
    def _demultiplex_golay12(self):
        """Correct and decode a Golay 12nt barcode"""
        self._demultiplex_encoded_barcode(decode_golay_12, 12)

    @requires(option='barcode_type', values='hamming_8')
    def _demultiplex_hamming8(self):
        """Correct and decode a Hamming 8nt barcode"""
        self._demultiplex_encoded_barcode(decode_hamming_8, 8)

    @requires(option='barcode_type', values='variable')
    def _demultiplex_other(self):
        """Decode a variable length barcode"""
        raise NotImplementedError

    #### use kwargs for method and bc_length
    def _demultiplex_encoded_barcode(self, method, bc_length):
        """Correct and decode an encoded barcode"""
        if self.state['Barcode'] is not None:
            from_sequence = False
            putative_bc = self.state['Barcode']
        else:
            from_sequence = True
            putative_bc = self.state['Sequence'][:bc_length]

        self.state['Original barcode'] = putative_bc

        if putative_bc in self.barcodes:
            self.state['Barcode errors'] = 0
            final_bc = putative_bc
            sample = self.barcodes[putative_bc]
        else:
            corrected, num_errors = method(putative_bc)
            final_bc = corrected
            self.state['Barcode errors'] = num_errors
            self.stats['barcode_corrected'] += 1
            sample = self.barcodes.get(corrected, None)

        self.state['Final barcode'] = final_bc

        if from_sequence:
            self.state['Sequence'] = self.state['Sequence'][bc_length:]

        if sample is None:
            self.failed = True
            self.stats['unknown_barcode'] += 1
        else:
            self.state['Sample'] = sample

    @requires(option='max_barcode_error', values=not_none)
    def _demultiplex_max_barcode_error(self):
        """Fail a sequence if it exceeds a max number of barcode errors"""
        bc_errors = self.options['max_barcode_error']
        if self.state['Barcode errors'] > bc_errors:
            self.failed = True
            self.stats['exceed_barcode_error'] += 1

    ### End demultiplex methods

    ### Start primer methods

    @requires(option='instrument_type', values='454')
    def _primer_instrument_454(self):
        """Check for a valid primer"""
        self._primer_check_forward()

    @requires(option='retain_primer')
    @requires(option='max_primer_mismatch', values=not_none)
    def _primer_check_forward(self):
        """Attempt to determine if the forward primer exists and trim if there

        Warning: this method may do an in place update on item if retain primer
        False.
        """
        seq = self.state['Sequence']
        qual = self.state['Qual']

        obs_barcode = self.state['Final barcode']
        exp_primers = self.primers.get(obs_barcode, None)

        if exp_primers is None:
            self.stats['unknown_primer_barcode_pair'] += 1
            self.failed = True
            return

        len_primer = len(exp_primers[0])

        obs_primer = seq[:len_primer]

        mm = np.array([count_mismatches(obs_primer, p) for p in exp_primers])

        if (mm > self.options['max_primer_mismatch']).all():
            self.failed = True
            self.stats['exceeds_max_primer_mismatch'] += 1
            return

        if not self.options['retain_primer']:
            seq = seq[len_primer:]
            if qual is not None:
                qual = qual[len_primer:]

        self.state['Forward primer'] = obs_primer
        self.state['Sequence'] = seq
        self.state['Qual'] = qual

    ### End primer methods

    ### Start sequence methods

    @requires(option='min_seq_len')
    def _sequence_length_check(self):
        """Checks minimum sequence length"""
        if len(self.state['Sequence']) < self.options['min_seq_len']:
            self.failed = True
            self.stats['min_seq_len'] += 1

    @requires(option='ambiguous_count')
    def _sequence_ambiguous_count(self):
        """Fail if the number of N characters is greater than threshold"""
        count = self.state['Sequence'].count('N')
        if count > self.options['ambiguous_count']:
            self.failed = True
            self.stats['ambiguous_count'] += 1

    ### End sequence methods
