#!/usr/bin/env python

from qiime.workflow.core import Workflow, requires, priority, _continuous
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import MinimalQualParser
from itertools import chain, izip

def _fasta_qual_strict(fasta_gen, qual_gen):
    for (seq_id, seq), (qual_id, qual) in izip(fasta_gen, qual_gen):
        if seq_id != qual_id:
            raise ValueError("%s is not equal to %s!" % (seq_id, qual_id))
        if len(seq) != len(qual):
            raise ValueError("%s is not equal length to %s!" % (seq_id,qual_id))

        yield (seq_id, seq, qual_id, qual)

def fasta_qual_iterator(fasta_fps, qual_fps=None):
    fasta_gens = chain(*map(MinimalFastaParser, fasta_fps))
    
    if qual_fps is not None:
        qual_gens = chain(*map(MinimalQualParser, qual_fps))
        gen = _fasta_qual_strict(fasta_gens, qual_gens)
    else:
        qual_gens = None
        gen = ((seq_id, seq, None, None) for seq_id, seq in fasta_gens)

    return gen

class QualFilterFastaWorkflow(Workflow):
    @priority(90)
    @requires(Option='min_seq_len', Values=_continuous)
    def wf_length_check(self, item):
        seq_id, seq, qual_id, qual = item

        if len(seq) < self.Options['min_seq_len']:
            self.Failed = True
            self.Stats['min_seq_len'] += 1
