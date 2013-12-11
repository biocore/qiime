#!/usr/bin/env python

from qiime.workflow.core import Workflow, requires, priority, _continuous
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import MinimalQualParser
from itertools import chain, izip

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

SEQ_ID_INDEX = 0
SEQ_INDEX = 1
QUAL_INDEX = 2

class QualFilterFastaWorkflow(Workflow):
    FinalState = {'fwd_primer':None,
                  'rev_primer':None,
                  'seq':None,
                  'qual':None,
                  'original_barcode':None,
                  'corrected_barcode':None}
    
    @priority(1000)
    @no_requirements
    def wf_init(self, item):
        # reset final state
        for k in self.FinalState:
            self.FinalState[k] = None

    @priority(90)
    @requires(Option='min_seq_len')
    def wf_length_check(self, item):
        """Checks minimum sequence length"""
        seq_id, seq, qual_id, qual = item

        if len(seq) < self.Options['min_seq_len']:
            self.Failed = True
            self.Stats['min_seq_len'] += 1

    @priority(89)
    @requires(IsValid=True)
    def wf_check_primer(self, item):
        """ """
        self._set_primers(item)

        self._local_align_forward_primer(item)
        self._
    @requires(IsValid=False, Option='ids_primers')
    def _set_primers(self, item):
        """ """
        seq_id = item[SEQ_ID_INDEX]

        if self.Options['suppress_sample_id_check']:
            primers = self.Options['ids_primers']['all_primers']
        else:
            seq_label = seq_id.split('_')[0]
            if seq_label not in self.Options['ids_primers']:
                self.Stats['seq_id_not_in_mapping'] += 1
                self.Failed = True
            else:
                primers = self.Options['ids_primers'][seq_label]
        else:
            primers = ids_primers['all_primers']

        self._primers = primers

    @requires(Option='local_align_forward_primer', Values=False)
    @requires(Option='max_primer_mismatch')
    @requires(Option='retain_primer', Values=False)
    def _count_mismatches(self, item):
        """ """
        seq = item[SEQ_INDEX]
        for primer in self._primers:
            exceeds_mismatch = count_mismatches(seq, primer,
                                       self.Options['max_primer_mismatch'])
            if not exceeds_mismatch:
                self.Stats['exceeds_max_primer_mismatch'] += 1
                if not retain_primer:
                    fasta_seq = fasta_seq[len(primer):]
                    qual_seq = qual_seq[len(primer):]
                failed = False
                break
    @requires(Option='local_align_forward_primer', Values=True)
    @requires(Option='max_primer_mismatch')
    def _local_align_forward_primer(self, item):
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
