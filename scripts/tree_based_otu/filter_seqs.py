#!/usr/bin/env python

"""Filter sequences. This code is safe to use on very large sequence files"""

from sys import argv
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.transform import exclude_chars, keep_chars
from collections import defaultdict

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"

def get_good_seqs(open_f, min_length, max_bad_base):
    toss_len_count = 0
    toss_base_count = 0
    keep_count = 0

    keep_f = keep_chars('ATGCU', case_sens=False)
    exclude_f = exclude_chars('ATGCU-.', case_sens=False)

    for seq_id, seq in MinimalFastaParser(open_f):
        seq_len = len(keep_f(seq))

        if seq_len < min_length:
            toss_len_count += 1
            continue

        bad_base_count = len(exclude_f(seq))
        if bad_base_count >= max_bad_base:
            toss_base_count += 1
            continue

        keep_count += 1
        yield (seq_id, seq)

    print 'tossed because of length: %d' % toss_len_count
    print 'tossed because of bad bases: %d' % toss_base_count
    print 'kept: %d' % keep_count

def test_self():
    seqs = ['>seq1\n', 'aaatttgggc..---gtcgtAAATGTGCCATTGT\n',
            '>seq2\n', 'aatgtggnnnntttagtgaagtnnattgaatgatgatgaa\n',
            '>seq3\n', ' aaaa------------a---aaaaaaaa------aaaa\n',
            '>seq4\n', 'aaaatgtgtga--agattgatnaagatgtatagaggtagtg\n']
    min_length = 20
    max_bad_base = 2
    exp = [('seq1', 'aaatttgggc..---gtcgtAAATGTGCCATTGT'),
           ('seq4', 'aaaatgtgtga--agattgatnaagatgtatagaggtagtg')]
    obs = list(get_good_seqs(seqs, min_length, max_bad_base))
    assert obs == exp

def main(args):
    open_f = open(args[1])
    min_length = int(args[2])
    min_bad_base = int(args[3])
    outfile = open(args[4],'w')

    my_buffer = []
    for seq_id, seq in get_good_seqs(open_f, min_length, min_bad_base):
        my_buffer.extend(['>' + seq_id, seq])
        if len(my_buffer) >= 500:
            outfile.write('\n'.join(my_buffer))
            outfile.write('\n')
            my_buffer = []
    outfile.write('\n'.join(my_buffer))
    outfile.close()

if __name__ == '__main__':
    main(argv)
