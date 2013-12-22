#!/usr/bin/env python

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"

import os
import tempfile

from cogent.util.unit_test import TestCase, main
from cogent.parse.binary_sff import parse_binary_sff
from qiime.make_per_library_sff import (
    filter_sff_reads, parse_id_list, combine_sff_data, make_per_library_sff,
    make_per_library_sff_with_sfffile,
)


class FunctionTests(TestCase):

    def setUp(self):
        test_dir = os.path.dirname(os.path.abspath(__file__))
        sff_filenames = ['F6AVWTA01.sff', 'F6AVWTA02.sff']
        self.sff_fps = [
            os.path.join(test_dir, 'test_support_files', 'F6AVWTA', fn)
            for fn in sff_filenames
        ]

    def test_filter_sff_reads(self):
        header = {}
        reads = [{'Name': 1}, {'Name': 2}, {'Name': 3}]

        header1, reads1 = filter_sff_reads(
            (header, reads), ids_to_keep=[1])
        self.assertEqual(header1, {'number_of_reads': 1})
        self.assertEqual(list(reads1), [{'Name': 1}])

        header23, reads23 = filter_sff_reads(
            (header, reads), ids_to_remove=[1])
        self.assertEqual(header23, {'number_of_reads': 2})
        self.assertEqual(list(reads23), [{'Name': 2}, {'Name': 3}])

        header2, reads2 = filter_sff_reads(
            (header, reads), ids_to_keep=[1, 2], ids_to_remove=[1])
        self.assertEqual(header2, {'number_of_reads': 1})
        self.assertEqual(list(reads2), [{'Name': 2}])

    def test_parse_id_list(self):
        ids = [
            '>asdfg fldk;s asdf\n',
            'mndmdf\t----',
        ]
        self.assertEqual(parse_id_list(ids), set(['asdfg', 'mndmdf']))

    def test_combine_sff_data(self):
        sff_datasets = [parse_binary_sff(open(fp)) for fp in self.sff_fps]
        observed_header, observed_reads = combine_sff_data(*sff_datasets)
        self.assertEqual(observed_header, combined_header)

        observed_reads = list(observed_reads)
        self.assertEqual(len(observed_reads), 40)
        observed_ids = [r['Name'] for r in observed_reads]
        self.assertEqual(observed_ids, combined_ids)

    def test_make_per_library_sff(self):
        id_list_file = tempfile.NamedTemporaryFile()
        id_list_file.write('GA202I001ER3QL\nGA202I001DBRNC\nGA202I001DJLC5\n')
        id_list_file.seek(0)

        make_per_library_sff(self.sff_fps, id_list_file.name)

        header, reads = parse_binary_sff(open(id_list_file.name + '.sff'))
        self.assertEquals(header, per_library_header)

        self.assertEqual(reads.next()['Name'], 'GA202I001ER3QL')
        self.assertEqual(reads.next()['Name'], 'GA202I001DBRNC')
        self.assertEqual(reads.next()['Name'], 'GA202I001DJLC5')
        self.assertRaises(StopIteration, reads.next)

    def test_make_per_library_sff_with_sfffile(self):
        id_list_file = tempfile.NamedTemporaryFile()
        id_list_file.write('GA202I001ER3QL\nGA202I001DBRNC\nGA202I001DJLC5\n')
        id_list_file.seek(0)

        make_per_library_sff_with_sfffile(self.sff_fps, id_list_file.name)

        header, reads = parse_binary_sff(open(id_list_file.name + '.sff'))
        # The index length varies between versions of sfftools
        del header['index_length']
        self.assertEquals(header, per_library_header_sfffile)

        self.assertEqual(reads.next()['Name'], 'GA202I001ER3QL')
        self.assertEqual(reads.next()['Name'], 'GA202I001DBRNC')
        self.assertEqual(reads.next()['Name'], 'GA202I001DJLC5')
        self.assertRaises(StopIteration, reads.next)


per_library_header_sfffile = {
    'header_length': 440,
    'flowgram_format_code': 1,
    'magic_number': 779314790,
    'number_of_flows_per_read': 400,
    'version': 1,
    'flow_chars':
    'TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG',
    'key_length': 4,
    'key_sequence': 'TCAG',
    'number_of_reads': 3,
    'index_offset': 5392,
}

per_library_header = {
    'header_length': 440,
    'flowgram_format_code': 1,
    'index_length': 0,
    'magic_number': 779314790,
    'number_of_flows_per_read': 400,
    'version': 1,
    'flow_chars':
    'TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG',
    'key_length': 4,
    'key_sequence': 'TCAG',
    'number_of_reads': 3,
    'index_offset': 0,
}

combined_header = {
    'flow_chars':
    'TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG',
    'flowgram_format_code': 1,
    'header_length': 440,
    'index_length': 0,
    'index_offset': 0,
    'key_length': 4,
    'key_sequence': 'TCAG',
    'magic_number': 779314790,
    'number_of_flows_per_read': 400,
    'number_of_reads': 40,
    'version': 1,
}

combined_ids = [
    'GA202I001ER3QL', 'GA202I001DBRNC', 'GA202I001DJLC5', 'GA202I001D4EKM',
    'GA202I001D7W79', 'GA202I001ECKP7', 'GA202I001DUJYF', 'GA202I001EEVFW',
    'GA202I001EZXR9', 'GA202I001DA5ZC', 'GA202I001B9RKV', 'GA202I001E1SL6',
    'GA202I001EGORF', 'GA202I001DRQUB', 'GA202I001COCJJ', 'GA202I001EK5WE',
    'GA202I001C3FPU', 'GA202I001B35KA', 'GA202I001C31UY', 'GA202I001CXAFL',
    'GA202I001DPW7B', 'GA202I001DIPMR', 'GA202I001EM8FL', 'GA202I001E4DQG',
    'GA202I001C11P2', 'GA202I001DWXK5', 'GA202I001DKJWO', 'GA202I001ES8TM',
    'GA202I001EHG2L', 'GA202I001D3VEC', 'GA202I001DMR8D', 'GA202I001D73WL',
    'GA202I001CRNW0', 'GA202I001B7WRX', 'GA202I001DPZYC', 'GA202I001CMX39',
    'GA202I001DGL8H', 'GA202I001DY7CG', 'GA202I001C76B0', 'GA202I001DJ1TL',
]


if __name__ == "__main__":
    main()
