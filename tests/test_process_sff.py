#!/usr/bin/env python
import os
import shutil
import tempfile
from cogent.util.unit_test import TestCase, main
from skbio.app.util import ApplicationNotFoundError
from qiime.process_sff import (
    make_flow_txt, make_fna, make_qual, prep_sffs_in_dir, convert_Ti_to_FLX,
    adjust_sff_cycles, check_sffinfo)
from cogent.parse.binary_sff import parse_binary_sff
from qiime.util import get_qiime_project_dir, qiime_open

"""Tests of the process_sff.py file.
"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = [
    "Rob Knight",
    "Greg Caporaso",
    "Kyle Bittinger",
    "Jesse Stombaugh",
    "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"


class TopLevelTests(TestCase):

    """Top-level tests of functions in process_sff"""

    def setUp(self):
        """Create temporary directory of SFF files"""

        # Cannot use get_qiime_project_dir() due to test errors in virtual box
        test_dir = os.path.dirname(os.path.abspath(__file__))
        sff_original_fp = os.path.join(test_dir, 'test_support_files',
                                       'test.sff')
        sff_original_gz_fp = os.path.join(test_dir, 'test_support_files',
                                          'test_gz.sff.gz')

        # copy sff file to working directory
        self.sff_dir = tempfile.mkdtemp()
        self.gz_sff_dir = tempfile.mkdtemp()
        self.sff_fp = os.path.join(self.sff_dir, 'test.sff')
        self.sff_gz_fp = os.path.join(self.gz_sff_dir, 'test_gz.sff.gz')
        shutil.copy(sff_original_fp, self.sff_fp)
        shutil.copy(sff_original_gz_fp, self.sff_gz_fp)

    def tearDown(self):
        shutil.rmtree(self.sff_dir)
        shutil.rmtree(self.gz_sff_dir)

    def test_adjust_sff_cycles(self):
        sff_data = parse_binary_sff(open(self.sff_fp))
        sff_gz_data = parse_binary_sff(qiime_open(self.sff_gz_fp))
        header, reads = adjust_sff_cycles(sff_data, 2)
        header_gz, reads_gz = adjust_sff_cycles(sff_gz_data, 2)
        expected_header = {
            'header_length': 48,
            'version': 1,
            'index_length': 0,
            'magic_number': 779314790,
            'number_of_flows_per_read': 8,
            'flowgram_format_code': 1,
            'flow_chars': 'TACGTACG',
            'index_offset': 0,
            'key_sequence': 'TCAG',
            'number_of_reads': 1,
            'key_length': 4,
        }
        self.assertEqual(header, expected_header)
        self.assertEqual(header_gz, expected_header)

        expected_read = {
            'name_length': 14,
            'Name': 'FA6P1OK01CGMHQ',
            'flowgram_values':
            [1.04, 0.0, 1.01, 0.0, 0.0, 0.95999999999999996, 0.0, 1.02],
            'clip_adapter_left': 0,
            'read_header_length': 32,
            'Bases': 'TCAG',
            'number_of_bases': 4,
            'flow_index_per_base': (1, 2, 3, 2),
            'clip_qual_left': 4,
            'clip_adapter_right': 0,
            'clip_qual_right': 4,
            'quality_scores': (32, 32, 32, 32),
        }
        reads = list(reads)
        reads_gz = list(reads_gz)
        self.assertEqual(len(reads), 1)
        self.assertEqual(len(reads_gz), 1)
        self.assertEqual(reads[0], expected_read)
        self.assertEqual(reads_gz[0], expected_read)

    def test_convert_Ti_to_FLX(self):
        """test_convert_Ti_to_FLX should do proper conversion from Ti to FLX"""
        sff_flx_fp = os.path.join(self.sff_dir, 'test_FLX.sff')
        sff_flx_gz_fp = os.path.join(self.gz_sff_dir, 'test_FLX_gz.sff')
        convert_Ti_to_FLX(self.sff_fp, sff_flx_fp)
        convert_Ti_to_FLX(self.sff_gz_fp, sff_flx_gz_fp)
        self.assertNotEqual(os.path.getsize(sff_flx_fp), 0)
        self.assertNotEqual(os.path.getsize(sff_flx_gz_fp), 0)

    def test_make_flow_txt(self):
        """test_make_flow_txt should make flowgram file as expected"""
        flow_fp = os.path.join(self.sff_dir, 'test.txt')
        flow_gz_fp = os.path.join(self.gz_sff_dir, 'test_gz.txt')
        make_flow_txt(self.sff_fp, flow_fp)
        make_flow_txt(self.sff_gz_fp, flow_gz_fp)
        self.assertEqual(open(flow_fp).read(), flow_txt)
        self.assertEqual(open(flow_gz_fp).read(), flow_txt)

    def test_make_fna(self):
        """test_make_fna should make fasta file as expected"""
        fna_fp = os.path.join(self.sff_dir, 'test.fna')
        fna_gz_fp = os.path.join(self.gz_sff_dir, 'test_gz.fna')
        make_fna(self.sff_fp, fna_fp)
        make_fna(self.sff_gz_fp, fna_gz_fp)
        self.assertEqual(open(fna_fp).read(), fna_txt)
        self.assertEqual(open(fna_gz_fp).read(), fna_txt)

    def test_make_qual(self):
        """test_make_qual should make qual file as expected"""
        qual_fp = os.path.join(self.sff_dir, 'test.qual')
        qual_gz_fp = os.path.join(self.gz_sff_dir, 'test_gz.qual')
        make_qual(self.sff_fp, qual_fp)
        make_qual(self.sff_gz_fp, qual_gz_fp)
        self.assertEqual(open(qual_fp).read(), qual_txt)
        self.assertEqual(open(qual_gz_fp).read(), qual_txt)

    def test_prep_sffs_in_dir(self):
        """test_prep_sffs_in_dir should make fasta/qual from sffs."""
        prep_sffs_in_dir(self.sff_dir, self.sff_dir, make_flowgram=True)
        prep_sffs_in_dir(self.gz_sff_dir, self.gz_sff_dir, make_flowgram=True)

        fna_fp = os.path.join(self.sff_dir, 'test.fna')
        fna_gz_fp = os.path.join(self.gz_sff_dir, 'test_gz.fna')
        self.assertEqual(open(fna_fp).read(), fna_txt)
        self.assertEqual(open(fna_gz_fp).read(), fna_txt)

        qual_fp = os.path.join(self.sff_dir, 'test.qual')
        qual_gz_fp = os.path.join(self.gz_sff_dir, 'test_gz.qual')
        self.assertEqual(open(qual_fp).read(), qual_txt)
        self.assertEqual(open(qual_gz_fp).read(), qual_txt)

        flow_fp = os.path.join(self.sff_dir, 'test.txt')
        flow_gz_fp = os.path.join(self.gz_sff_dir, 'test_gz.txt')
        self.assertEqual(open(flow_fp).read(), flow_txt)
        self.assertEqual(open(flow_gz_fp).read(), flow_txt)

    def test_prep_sffs_in_dir_FLX(self):
        """test_prep_sffs_in_dir should convert to FLX read lengths."""
        output_dir = tempfile.mkdtemp()
        gz_output_dir = tempfile.mkdtemp()

        prep_sffs_in_dir(
            self.sff_dir, output_dir, make_flowgram=True, convert_to_flx=True)
        prep_sffs_in_dir(
            self.gz_sff_dir, gz_output_dir, make_flowgram=True, convert_to_flx=True)

        fna_fp = os.path.join(output_dir, 'test_FLX.fna')
        fna_gz_fp = os.path.join(gz_output_dir, 'test_gz_FLX.fna')
        self.assertEqual(open(fna_fp).read(), fna_txt)
        self.assertEqual(open(fna_gz_fp).read(), fna_txt)

        qual_fp = os.path.join(output_dir, 'test_FLX.qual')
        qual_gz_fp = os.path.join(gz_output_dir, 'test_gz_FLX.qual')
        self.assertEqual(open(qual_fp).read(), qual_txt)
        self.assertEqual(open(qual_gz_fp).read(), qual_txt)

        flow_fp = os.path.join(output_dir, 'test_FLX.txt')
        flow_gz_fp = os.path.join(gz_output_dir, 'test_gz_FLX.txt')
        self.assertEqual(open(flow_fp).read(), flx_flow_txt)
        self.assertEqual(open(flow_gz_fp).read(), flx_flow_txt)

        shutil.rmtree(output_dir)
        shutil.rmtree(gz_output_dir)

    def test_prep_sffs_in_dir_no_trim(self):
        """test_prep_sffs_in_dir should use the no_trim option only if sffinfo exists."""
        output_dir = tempfile.mkdtemp()
        gz_output_dir = tempfile.mkdtemp()

        try:
            check_sffinfo()
            perform_test = True
        except:
            perform_test = False

        if perform_test:
            prep_sffs_in_dir(self.sff_dir, output_dir, make_flowgram=False,
                             convert_to_flx=False, use_sfftools=True,
                             no_trim=True)

            fna_fp = os.path.join(output_dir, 'test.fna')

            self.assertEqual(open(fna_fp).read(), fna_notrim_txt)

            qual_fp = os.path.join(output_dir, 'test.qual')
            self.assertEqual(open(qual_fp).read(), qual_notrim_txt)

            self.assertRaises(TypeError, "gzipped SFF", prep_sffs_in_dir,
                              self.gz_sff_dir, gz_output_dir, make_flowgram=False,
                              convert_to_flx=False, use_sfftools=True,
                              no_trim=True)

            shutil.rmtree(output_dir)
            shutil.rmtree(gz_output_dir)


flow_txt = """\
Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  1504
  Index Length:  706
  # of Reads:    1
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG

>FA6P1OK01CGMHQ
  Run Prefix:   R_2008_05_28_17_11_38_
  Region #:     1
  XY Location:  0892_1356

  Read Header Len:  32
  Name Length:      14
  # of Bases:       77
  Clip Qual Left:   5
  Clip Qual Right:  52
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.00	1.01	0.00	0.00	0.96	0.00	1.02	0.00	1.02	0.00	0.00	0.99	0.00	1.00	0.00	1.00	0.00	0.00	1.00	0.00	1.10	0.00	1.08	0.00	0.00	1.46	0.00	0.88	0.18	0.00	2.69	1.01	0.08	0.96	0.00	0.02	0.92	0.08	0.00	0.98	0.68	0.00	0.89	0.00	0.00	1.15	0.00	1.13	0.00	0.02	1.12	0.05	0.15	1.84	0.00	1.10	0.00	2.47	0.96	0.86	1.06	0.00	1.96	0.12	0.93	0.13	1.65	1.06	0.06	0.00	0.99	0.00	0.00	1.87	0.44	1.08	0.00	3.25	0.09	0.97	0.50	1.00	1.72	0.07	0.00	0.92	0.58	0.00	0.00	0.59	0.06	0.11	0.09	0.07	0.06	0.16	0.00	0.24	0.03	0.00	0.12	0.06	0.16	0.00	0.18	0.00	0.00	0.14	0.00	0.15	0.00	0.18	0.00	0.03	0.14	0.03	0.13	0.01	0.19	0.00	0.02	0.33	0.05	0.00	0.16	0.10	0.35	0.01	0.21	0.04	0.09	0.18	0.13	0.19	0.00	0.10	0.51	0.26	0.00	0.23	0.19	0.27	0.01	0.29	0.05	0.14	0.17	0.16	0.18	0.27	0.09	0.26	0.10	0.18	0.23	0.15	0.22	0.13	0.37	0.11	0.11	0.26	0.59	0.14	0.06	0.33	0.34	0.26	0.05	0.27	0.44	0.19	0.10	0.35	0.27	0.15	0.34	0.28	0.45	0.14	0.16	0.34	0.27	0.12	0.07	0.25	0.18	0.12	0.04	0.23	0.16	0.12	0.05	0.20	0.16	0.11	0.03	0.21	0.16	0.10	0.02	0.21	0.16	0.12	0.02	0.20	0.15	0.10	0.02	0.23	0.15	0.11	0.02	0.22	0.14	0.09	0.02	0.20	0.13	0.09	0.01	0.19	0.13	0.08	0.02	0.17	0.12	0.08	0.03	0.17	0.09	0.08	0.01	0.14	0.09	0.07	0.01	0.15	0.09	0.06	0.01	0.13	0.08	0.06	0.00	0.13	0.08	0.05	0.02	0.12	0.07	0.05	0.01	0.11	0.07	0.05	0.00	0.10	0.07	0.05	0.01	0.11	0.08	0.04	0.00	0.10	0.06	0.05	0.01	0.09	0.06	0.04	0.01	0.08	0.07	0.05	0.00	0.08	0.06	0.05	0.00	0.09	0.06	0.04	0.00	0.09	0.06	0.04	0.01	0.08	0.06	0.04	0.00	0.09	0.06	0.03	0.00	0.09	0.06	0.02	0.00	0.09	0.06	0.04	0.00	0.08	0.05	0.03	0.00	0.07	0.05	0.02	0.00	0.08	0.04	0.03	0.00	0.07	0.04	0.03	0.00	0.07	0.05	0.02	0.00	0.07	0.05	0.02	0.00	0.06	0.04	0.02	0.00	0.06	0.03	0.03	0.00	0.08	0.02	0.00	0.00	0.07	0.03	0.01	0.00	0.06	0.03	0.02	0.00	0.05	0.03	0.02	0.00	0.05	0.03	0.01	0.00	0.06	0.02	0.00	0.00	0.05	0.01	0.01	0.00	0.04	0.01	0.01	0.00	0.04	0.01	0.01	0.00	0.05	0.01	0.00	0.00	0.04	0.02	0.01	0.00	0.03	0.02	0.01	0.00	0.03	0.01	0.00	0.00	0.03	0.00	0.00	0.00	0.03	0.00	0.00	0.00	0.02	0.00
Flow Indexes:	1	3	6	8	10	13	15	17	20	22	24	27	29	32	32	32	33	35	38	41	42	44	47	49	52	55	55	57	59	59	60	61	62	64	64	66	68	68	69	72	75	75	77	79	79	79	81	82	83	84	84	87	88	91	102	126	130	138	140	145	153	157	161	164	166	171	175	179	183	187	191	195	199	203	211	215	219
Bases:	tcagATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGCgcnnnannnnngnnnnnnnnnnnnn
Quality Scores:	32	32	32	32	32	32	32	32	32	32	32	25	25	21	21	21	28	32	32	31	30	30	32	32	32	33	31	25	18	18	20	18	32	30	28	23	22	22	24	28	18	19	18	16	16	16	17	18	13	17	27	21	20	21	0	0	0	17	0	0	0	0	0	17	0	0	0	0	0	0	0	0	0	0	0	0	0
"""

fna_txt = """\
>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_
ATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGC
"""

qual_txt = """\
>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_
32 32 32 32 32 32 32 25 25 21 21 21 28 32 32 31 30 30 32 32 32 33 31 25 18 18 20 18 32 30 28 23 22 22 24 28 18 19 18 16 16 16 17 18 13 17 27 21
"""

fna_notrim_txt = """\
>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_
tcagATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGCgcnnnann
nnngnnnnnnnnnnnnn
"""

qual_notrim_txt = """\
>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_
32 32 32 32 32 32 32 32 32 32 32 25 25 21 21 21 28 32 32 31 30 30 32 32 32 33 31 25 18 18 20 18 32 30 28 23 22 22 24 28 18 19 18 16 16 16 17 18 13 17 27 21 20 21 0 0 0 17 0 0
0 0 0 17 0 0 0 0 0 0 0 0 0 0 0 0 0
"""

# same as other flow_txt, but index_offset and index_length are now 0,
# since we don't yet support the SFF index section.
flx_flow_txt = """\
Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  0
  Index Length:  0
  # of Reads:    1
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG

>FA6P1OK01CGMHQ
  Run Prefix:   R_2008_05_28_17_11_38_
  Region #:     1
  XY Location:  0892_1356

  Read Header Len:  32
  Name Length:      14
  # of Bases:       77
  Clip Qual Left:   5
  Clip Qual Right:  52
  Clip Adap Left:   0
  Clip Adap Right:  0

Flowgram:	1.04	0.00	1.01	0.00	0.00	0.96	0.00	1.02	0.00	1.02	0.00	0.00	0.99	0.00	1.00	0.00	1.00	0.00	0.00	1.00	0.00	1.10	0.00	1.08	0.00	0.00	1.46	0.00	0.88	0.18	0.00	2.69	1.01	0.08	0.96	0.00	0.02	0.92	0.08	0.00	0.98	0.68	0.00	0.89	0.00	0.00	1.15	0.00	1.13	0.00	0.02	1.12	0.05	0.15	1.84	0.00	1.10	0.00	2.47	0.96	0.86	1.06	0.00	1.96	0.12	0.93	0.13	1.65	1.06	0.06	0.00	0.99	0.00	0.00	1.87	0.44	1.08	0.00	3.25	0.09	0.97	0.50	1.00	1.72	0.07	0.00	0.92	0.58	0.00	0.00	0.59	0.06	0.11	0.09	0.07	0.06	0.16	0.00	0.24	0.03	0.00	0.12	0.06	0.16	0.00	0.18	0.00	0.00	0.14	0.00	0.15	0.00	0.18	0.00	0.03	0.14	0.03	0.13	0.01	0.19	0.00	0.02	0.33	0.05	0.00	0.16	0.10	0.35	0.01	0.21	0.04	0.09	0.18	0.13	0.19	0.00	0.10	0.51	0.26	0.00	0.23	0.19	0.27	0.01	0.29	0.05	0.14	0.17	0.16	0.18	0.27	0.09	0.26	0.10	0.18	0.23	0.15	0.22	0.13	0.37	0.11	0.11	0.26	0.59	0.14	0.06	0.33	0.34	0.26	0.05	0.27	0.44	0.19	0.10	0.35	0.27	0.15	0.34	0.28	0.45	0.14	0.16	0.34	0.27	0.12	0.07	0.25	0.18	0.12	0.04	0.23	0.16	0.12	0.05	0.20	0.16	0.11	0.03	0.21	0.16	0.10	0.02	0.21	0.16	0.12	0.02	0.20	0.15	0.10	0.02	0.23	0.15	0.11	0.02	0.22	0.14	0.09	0.02	0.20	0.13	0.09	0.01	0.19	0.13	0.08	0.02	0.17	0.12	0.08	0.03	0.17	0.09	0.08	0.01	0.14	0.09	0.07	0.01	0.15	0.09	0.06	0.01	0.13	0.08	0.06	0.00	0.13	0.08	0.05	0.02	0.12	0.07	0.05	0.01	0.11	0.07	0.05	0.00	0.10	0.07	0.05	0.01	0.11	0.08	0.04	0.00	0.10	0.06	0.05	0.01	0.09	0.06	0.04	0.01	0.08	0.07	0.05	0.00	0.08	0.06	0.05	0.00	0.09	0.06	0.04	0.00	0.09	0.06	0.04	0.01	0.08	0.06	0.04	0.00	0.09	0.06	0.03	0.00	0.09	0.06	0.02	0.00	0.09	0.06	0.04	0.00	0.08	0.05	0.03	0.00	0.07	0.05	0.02	0.00	0.08	0.04	0.03	0.00	0.07	0.04	0.03	0.00	0.07	0.05	0.02	0.00	0.07	0.05	0.02	0.00	0.06	0.04	0.02	0.00	0.06	0.03	0.03	0.00	0.08	0.02	0.00	0.00	0.07	0.03	0.01	0.00	0.06	0.03	0.02	0.00	0.05	0.03	0.02	0.00	0.05	0.03	0.01	0.00	0.06	0.02	0.00	0.00	0.05	0.01	0.01	0.00	0.04	0.01	0.01	0.00	0.04	0.01	0.01	0.00	0.05	0.01	0.00	0.00	0.04	0.02	0.01	0.00	0.03	0.02	0.01	0.00	0.03	0.01	0.00	0.00	0.03	0.00	0.00	0.00	0.03	0.00	0.00	0.00	0.02	0.00
Flow Indexes:	1	3	6	8	10	13	15	17	20	22	24	27	29	32	32	32	33	35	38	41	42	44	47	49	52	55	55	57	59	59	60	61	62	64	64	66	68	68	69	72	75	75	77	79	79	79	81	82	83	84	84	87	88	91	102	126	130	138	140	145	153	157	161	164	166	171	175	179	183	187	191	195	199	203	211	215	219
Bases:	tcagATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGCgcnnnannnnngnnnnnnnnnnnnn
Quality Scores:	32	32	32	32	32	32	32	32	32	32	32	25	25	21	21	21	28	32	32	31	30	30	32	32	32	33	31	25	18	18	20	18	32	30	28	23	22	22	24	28	18	19	18	16	16	16	17	18	13	17	27	21	20	21	0	0	0	17	0	0	0	0	0	17	0	0	0	0	0	0	0	0	0	0	0	0	0
"""

if __name__ == '__main__':
    main()
