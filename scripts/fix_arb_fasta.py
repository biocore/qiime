#!/usr/bin/env python
#simple script to fix arb fasta format by repairing incorrect line break chars
#and stripping spaces and replacing . with - chars.

from cogent.parse.fasta import MinimalFastaParser
from sys import argv
for label, seq in MinimalFastaParser(open(argv[1], 'U')):
    print '>%s\n%s' % (label, seq.replace(' ','').replace('.','-'))

