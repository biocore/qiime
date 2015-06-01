#!/usr/bin/env python
# file make_library_id_lists.py: make id list for each lib from fasta file
from optparse import OptionParser
from string import strip
from os.path import exists, join
from os import makedirs
from collections import defaultdict

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"


def get_ids(lines, field, bad_ids=None, debug=False):
    """Make dict of lib:ids"""
    result = defaultdict(list)
    for line in lines:
        if line.startswith('>'):
            fields = map(strip, line[1:].split())
            label = fields[0]
            if not '_' in label:  # no lib specified
                continue
            lib, id_ = label.rsplit('_', 1)
            if bad_ids and label in bad_ids:
                if debug:
                    print "Excluded bad id: %s" % label
            else:
                result[lib].append(fields[field])
    return result


def get_first_id(lines):
    """Gets first fasta id from each line in lines"""
    result = set()
    for line in lines:
        if line.startswith('>'):
            result.add(line[1:].split()[0])
    return result
