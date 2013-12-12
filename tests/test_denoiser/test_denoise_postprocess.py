#!/usr/bin/env python

"""tests for denoiser_postprocess functions
"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Jens Reeder", "Rob Knight"]#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"

from os import remove, rmdir, mkdir
from cogent.util.unit_test import TestCase, main
from qiime.util import get_tmp_filename

#import as _main to not interfere with TestCase.main
from qiime.denoiser.denoise_postprocess import extract_read_to_sample_mapping,\
    combine_mappings

class DenoiserTests(TestCase):

    def setUp(self):
        self.tmp_dir= None

    def tearDown(self):
        """remove tmp files"""        
      
        if self.tmp_dir:
            remove(self.tmp_dir+"denoised_otu_map.txt")
            remove(self.tmp_dir+"denoised_all.fasta")
            rmdir(self.tmp_dir)

    def test_extract_read_to_sample_mapping(self):
        """extract_read_to_sample_mapping pulls info from label line"""
        
        labels = ['S160_1 E86FECS01DW5V4 orig_bc=CAGTACGATCTT new_bc=CAGTACGATCTT bc_diffs=0',
                  'S160_2 E86FECS01DW5V5 orig_bc=CAGTACGATCTT new_bc=CAGTACGATCTT bc_diffs=0']
         
        expected = {'E86FECS01DW5V4': 'S160_1',
                    'E86FECS01DW5V5': 'S160_2'}
                     
        self.assertEqual(extract_read_to_sample_mapping(labels),
                         expected)
        

    def test_combine_mappings(self):
        """combine_mappings works as expected"""

        self.tmp_dir = get_tmp_filename(tmp_dir = "./", suffix="/")
        mkdir(self.tmp_dir)

        combine_mappings(fasta, denoiser_mapping, denoised_seqs, otu_picker_map, self.tmp_dir)
        
        observed_otu_map = "".join(list(open(self.tmp_dir +"/denoised_otu_map.txt")))
        
        expected_otu_map ="""1:\tS1_1\tS1_2\tS2_4\tS2_5
2:\tS2_3\tS1_6
"""
        self.assertEqual(observed_otu_map, expected_otu_map)

        observed_fasta = "".join(list(open(self.tmp_dir+"/denoised_all.fasta")))
        expected_fasta = """>S1_1 Read1
AAA
>S1_2 Read2
TTT
>S2_3 Read3
GGG
"""
        self.assertEqual(observed_fasta, expected_fasta)


fasta = """>S1_1 Read1
AAA
>S1_2 Read2
TTT
>S2_3 Read3
GGG
>S2_4 Read4
GG
>S2_5 Read5
TTT
>S1_6 Read6
AA""".split("\n")

denoised_seqs = """>Read1
AAA
>Read2
TTT
>Read3
GGG""".split("\n")
        
denoiser_mapping = """Read1:\tRead4\tRead5
Read2:
Read3:\tRead6""".split("\n")

otu_picker_map = """1: Read1 Read2
2: Read3""".split("\n")


if __name__ == "__main__":
    main()
