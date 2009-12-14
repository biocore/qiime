#!/usr/bin/env python
#file prep_sff.py
from os import listdir, system
from os.path import splitext, join
"""Converts directory of sff files into fasta and qual files.

Requires that 454's off-instrument apps are on your path.
"""
def make_fna(filename):
    """Makes fna file from sff file."""
    system('sffinfo -s %s > %s.fna' % (filename, splitext(filename)[0]))

def make_qual(filename):
    """Makes qual file from sff file."""
    system('sffinfo -q %s > %s.qual' % (filename, splitext(filename)[0]))

def prep_sffs_in_dir(pathname):
    """Converts all sffs in dir to fasta/qual."""
    for name in listdir(pathname):
        if name.endswith('.sff'):
            make_fna(join(pathname,name))
            make_qual(join(pathname,name))


if __name__ == '__main__':
    from sys import argv
    pathname = argv[1]
    prep_sffs_in_dir(pathname)
