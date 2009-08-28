#!/usr/bin/env python
#filter_otu_table

from qiime.parse import parse_otus
from sys import argv, stdout
from numpy import array

min_otu_count = 5
min_otu_samples = 2

otus = open(argv[1])
if len(argv) > 2:
    min_otu_count = int(argv[2])
if len(argv) < 3:
    min_otu_samples = int(argv[3])
for line in otus:
    if line.startswith('#'):
        stdout.write(line)
    else:
        fields = line.split('\t')
        try:
            vals = array(map(int, fields[1:]), dtype=int)
        except ValueError:
            vals = array(map(int, fields[1:-1]), dtype=int)
        if vals.sum() >= min_otu_count and (vals > 0).sum() > min_otu_samples:
            stdout.write(line)


