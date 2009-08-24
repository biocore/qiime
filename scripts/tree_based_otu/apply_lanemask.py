#!/usr/bin/env python
from numpy import nonzero, array, fromstring
from cogent.parse.fasta import MinimalFastaParser

def lanemask_to_positions(s):
    """Converts lanemask to valid positions."""
    return nonzero(array(map(int, s)))

def lanemasked_positions(s, p):
    """Gets valid positons in p from s."""
    return fromstring(s,'b')[p].tostring()

def lanemasked_fasta(lanemask, lines, min_index=None, max_index=None):
    """Applies lanemask to fasta-formatted file, yielding masked seqs.
    
    Note: range returned includes min_index but not max_index, as is standard.
    """
    p = lanemask_to_positions(lanemask)[0]
    if min_index:
        i = p.searchsorted(min_index)
        p = p[i:]
    if max_index:
        i = p.searchsorted(max_index)
        p = p[:i]

    for k, v in MinimalFastaParser(lines):
        yield '>%s\n%s' % (k, lanemasked_positions(v,p))

def test():
    lanemask = '1001011000001'
    s = '>a\nabcdefghijklmnopqrs\n>b\nabcdefghijklmnopqrs\n'.splitlines()
    for i in lanemasked_fasta(lanemask, s, 3, 7): print i

if __name__ == '__main__':
    from sys import argv
    #default values: 188, 1900
    lanemask = open(argv[1]).read().strip()
    infile = open(argv[2],'U')
    if len(argv) > 3:
        min_index, max_index = int(argv[3]), int(argv[4])
    else:
        min_index, max_index = None, None
    for i in lanemasked_fasta(lanemask, infile, min_index, max_index):
        print i



