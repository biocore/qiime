#!/usr/bin/env python
# File created on 16 Feb 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"


from tempfile import mkdtemp
from unittest import TestCase, main
from qiime.submit_to_mgrast import parse_and_submit_params, post_multipart,\
    encode_multipart_formdata, get_content_type
from os import mkdir, remove, removedirs, path, listdir
from qiime.util import get_qiime_project_dir
from glob import glob


class TopLevelTests(TestCase):

    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""

        qiime_dir = get_qiime_project_dir()

        self.key = 'qiime_test'
        self.project_id = 'qiime_test'
        self.sample_id = 'qiime_sample1'
        self.params = [('key', self.key), ('sample', self.sample_id),
                       ('project', self.project_id)]
        test_dir = path.dirname(path.abspath(__file__))
        self.seq_file = path.join(test_dir, 'test_support_files',
                                  'qiime_tutorial_split_lib_seqs_subset.fna')
        self.output_dir = mkdtemp()
        self.sample_file = [('file', 'qiime_test.fna', fasta_example)]
        self._paths_to_clean_up = []
        self._dirs_to_clean_up = []

        # make the webfile directory
        try:
            mkdir(self.output_dir)
        except OSError:
            pass

        # define directory to clean up
        self._dirs_to_clean_up = [self.output_dir]

    def tearDown(self):
        map(remove, self._paths_to_clean_up)
        map(removedirs, self._dirs_to_clean_up)

    def test_parse_and_submit_params(self):
        """parse_and_submit_params: separates the split-library sequences \
into separate fasta files and posts the sequence data to mg-rast"""

        obs = parse_and_submit_params(self.key, self.project_id, self.seq_file,
                                      self.output_dir, submit_to_server=False)

        self.assertEqual(obs, mg_rast_exp_log)

        for i in glob('%s/*' % (self.output_dir)):
            self._paths_to_clean_up.append(path.join(self.output_dir, i))

    def test_post_multipart(self):
        """post_multipart: sends the data as a multipart form"""

        obs = post_multipart('metagenomics.nmpdr.org', self.params,
                             self.sample_file, submit_to_server=False)

        self.assertEqual(obs, post_body_exp)

    def test_encode_multipart_formdata(self):
        """encode_multipart_formdata: encodes the params and file to multipart
        posting format"""

        obs1, obs2 = encode_multipart_formdata(self.params, self.sample_file)

        self.assertEqual(obs1, 'multipart/form-data')
        self.assertEqual(obs2, post_body_exp)

    def test_get_content_type(self):
        """get_content_type: this tries to determine the file type based on \
the filename"""

        # the seq_file is .fna so should return text/plain
        obs1 = get_content_type(self.seq_file)
        self.assertEqual(obs1, 'text/plain')

        # fasta file should return text/plain
        obs2 = get_content_type('seqs.fasta')
        self.assertEqual(obs2, 'text/plain')

        # text file should return text/plain
        obs3 = get_content_type('seqs.txt')
        self.assertEqual(obs3, 'text/plain')

        # html file should return text/html
        obs4 = get_content_type('seqs.html')
        self.assertEqual(obs4, 'text/html')

fasta_example = """\
>PC.354_202 FLP3FBN01CUQNU orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTT\n\
>PC.355_85 FLP3FBN01AIELL orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n\
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG"\n\
"""

post_body_exp = """\
--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
qiime_sample1\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="qiime_test.fna"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.354_202 FLP3FBN01CUQNU orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTT\n\
>PC.355_85 FLP3FBN01AIELL orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n\
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG"\n\
\r\n\
--xYzZY--\r\n\
"""

mg_rast_exp_log = """\
<h3>The following jobs were submitted to MG-RAST.</h3>\n<table border=1><tr><th>Fasta File</th><th>Job ID</th>\n<th>md5</th></tr>\n</table><br><h3 style="color:red">--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
PC.354\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="PC.354.fasta"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.354_202 FLP3FBN01CUQNU orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTATCCACGCGTT\n\
>PC.354_224 FLP3FBN01BZOEE orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTCGGTGGGCCGTTACCCCGCCGACTAGCTAATGCGCCGCATGGCCATCCGCAGCCGATAAATCTTTAAACATCGGGAGATGCCTCCCAACGTTGTTACGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGCTGCGGGCAGGTTCCATACGTGTTACTCACCCGTGCG\n\r\n\
--xYzZY--\r\n\
</h3>\n</table><br><h3 style="color:red">--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
PC.355\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="PC.355.fasta"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.355_85 FLP3FBN01AIELL orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n\
TTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGCGGACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCG\n\
>PC.355_66 FLP3FBN01BYAA8 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCTCACCAACCAGCTAATCAGACGCGGGTCCATCTTACACCACCGGAGTTTTTCACACCGAACCATGCGGTTCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTGCAAGGCAGGTTACCCACGCGTTACTCACCCGTCCG\n\r\n\
--xYzZY--\r\n\
</h3>\n</table><br><h3 style="color:red">--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
PC.356\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="PC.356.fasta"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.356_249 FLP3FBN01BXX60 orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0\n\
CTGGGCCGTGTCTCAGTCCCAGTGTGGCCGTCCACCCTCTCAGGCCGGCTGCTGATCGTCGCTTTGGTGGGCCGTTACCCCGCCAACTGGCTAATCGGACGCGGATCCATCGTATGCCGATAAATCTTTTCACACCAGACCATGCGATCCTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTATACGGCAGGTTACCCACGCGTTACT\n\
>PC.356_243 FLP3FBN01DRT2O orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACT\n\r\n\
--xYzZY--\r\n\
</h3>\n</table><br><h3 style="color:red">--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
PC.481\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="PC.481.fasta"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.481_219 FLP3FBN01EFY7W orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCTTTACCCCGCCAACAAGCTAATCAGACGCGGATCCATCGTATACCACCAAAAGCTTTAGCTTTTTGTTTTCCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCT\n\
>PC.481_228 FLP3FBN01ED5UR orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCTGTTACCCCGCCAACCAGCTAATCAGACGCGGATCCATCGTATACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTT\n\r\n\
--xYzZY--\r\n\
</h3>\n</table><br><h3 style="color:red">--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
PC.593\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="PC.593.fasta"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.593_314 FLP3FBN01AGF7L orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\n\
TTGGACCGTGTCTCAGTTCCAATGTGGCCGATCACCCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTAAGCCGTTACCCTACCAACTAGCTAATCAGACGCGGGTCCATCCTGTACCGCAAAAGCTTTGATACTTCTACCATGCGATAAAGTATTTTATCTCGTATTAGCATACCTTTCGGTATGTTATCCGTGTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCG\n\
>PC.593_342 FLP3FBN01D9HWD orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nTTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACGCATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATGCGCCATAAGTCCATCCTCTACCAGTGCCTTGCAGCACTTTTAATACGGTCACCATGCAGTGTCCCTACCTATGCGGTCTTAGCTGCCGTTTCCAGCAGTTATCCCCCTGTAAAGGCCAGGTT\n\r\n\
--xYzZY--\r\n\
</h3>\n</table><br><h3 style="color:red">--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
PC.607\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="PC.607.fasta"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.607_209 FLP3FBN01BAT5U orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGCTGCCCCGCCAACAAGCTAATCAGACGCGGGCCCCTCCCATACCGCCGGAACTTTCCCCAGAAAGGCATGCGCCTCCCTGGTTTATGCGGTATTAGCAGCCGTTTCCGGCTGTTATCCCCCTGTATGGGGCAGGTTGCCCACGCGTTACTCACCCGTCCGCC\n\
>PC.607_197 FLP3FBN01B5DKJ orig_bc=AACTGTGCGTAC new_bc=AACTGTGCGTAC bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGTCCATCCTATACCACCGGAGTTTTTCACACCGGAGCATGCGCTCCTGTGCGCTTATGCGGTATTAACAGTCGTTTCCAACTGTTATCCCCCTGTATAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACT\n\r\n\
--xYzZY--\r\n\
</h3>\n</table><br><h3 style="color:red">--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
PC.634\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="PC.634.fasta"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.634_11 FLP3FBN01CPD70 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n\
CTGGACCGTGTCTCAGTTCCAATGTGGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCG\n\
>PC.634_14 FLP3FBN01AM0P3 orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0\n\
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGGCCGCCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCGCCAACCAGCTAATCAGACGCGGGCCCATCCCGTACCACCGGAGTTTTTCACACTGCTTCATGCGAAGCTGTGCGCTTATGCGGTATTAGCACCTATTTCTAAGTGTTATCCCCCAGTATACGGCAGGTTCTCCACGCGTTACTCACCCG\n\r\n\
--xYzZY--\r\n\
</h3>\n</table><br><h3 style="color:red">--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
PC.635\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="PC.635.fasta"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.635_130 FLP3FBN01EA91G orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\n\
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCCTTGGTGGGCCGTTACCCCTCCAACCAGCTAATCAGACGCGGGTCCATCCTGTACCACCGGAGTTTTTCACACTGTACCATGCGGTACTGTGCGCTTATGCGGTTTTAGCACCTATTTCTAAGTGTTATCCCCCTGTACAGGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCACTAAG\n\
>PC.635_131 FLP3FBN01B06QS orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\n\
TTGGTCCGTGTCTCAGTACCAATGTGGGGGGTTAACCTCTCAGTCCCCCTATGTATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATACAACGCATGCCCATCTTCCACCACAAAAAGCTTTCAACCCAGAGAGATGCCTCTCCGAATTATATGGGGTATTAGTACCAATTTCTCAGTGTTATCCCCCTGTGAAAGGTAGGTTGCATACGCGTTACGCACCCGTCCGCCGGTCG\n\r\n\
--xYzZY--\r\n\
</h3>\n</table><br><h3 style="color:red">--xYzZY\r\n\
Content-Disposition: form-data; name="key"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="sample"\r\n\
\r\n\
PC.636\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="project"\r\n\
\r\n\
qiime_test\r\n\
--xYzZY\r\n\
Content-Disposition: form-data; name="file"; filename="PC.636.fasta"\r\n\
Content-Type: text/plain\r\n\
\r\n\
>PC.636_100 FLP3FBN01D7AVV orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0\n\
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGGCTTGGTGGTCCGTTACACCGCCAACTACCTAGTGCGACGCATGCCCATCCGCTACCGGATCGCTCCTTTGGAATCCCGGGGATGTCCCCGGAACTCGTTATGCGGTATTAGACGGAATTTCTTCCGCTTATCCCCCTGTAGCGGGCAGGTTGCATACGTGTTACTCACCCGTGCGCCGGTCGCCGG\n\
>PC.636_113 FLP3FBN01BOOHE orig_bc=ACGGTGAGTGTC new_bc=ACGGTGAGTGTC bc_diffs=0\n\
CTGGACCGTGTCTCAGTTCCAGTGTGGCCGATCACCCTCTCAGGTCGGCTACGTATCGTTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATACGGCGCGGGTCCATCTATAAGTGACAGCCGAAACCGTCTTTCAACATTGAACCATGCGGTTCAATATATTATCCGGTATTAGCCCCGGTTTCCCGGAGTTATCCCAGTCTTATAGGTAGGTTACCCACGTGTTACTCACCCGTGCGCCGGTCGCCGG\n\r\n\
--xYzZY--\r\n\
</h3>\n</table>\
"""


if __name__ == "__main__":
    main()
