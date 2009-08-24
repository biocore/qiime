#!/usr/bin/env python
#quick and dirty script to grab lines from greengenes file matching fasta

from cogent.parse.fasta import MinimalFastaParser
from sys import argv

ids = []

for k, v in MinimalFastaParser(open(argv[1], 'U')):
    ids.append(k.split()[0])

ids = set(ids) #in case of duplicates, who knows...?

for line in open(argv[2], 'U'):
    if not line.strip():
        continue
    id_, fields = line.split('\t', 1)
    if id_ in ids:
        print line.rstrip()


